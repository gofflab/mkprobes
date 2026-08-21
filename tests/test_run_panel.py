"""Tests for the panel batch driver (mkprobes run-panel)."""

import json
from concurrent.futures import Future
from pathlib import Path
from unittest.mock import patch

import pytest
from click.testing import CliRunner

from mkprobes.run_panel import (
    DEFAULT_RESTRICTION,
    final_parquet,
    find_missing_final,
    load_acceptable,
    load_worklist,
    run_panel,
    run_panel_cli,
)

CODEBOOK = {"Och.687.1": [1, 2, 3], "Och.958.1": [4, 5, 6], "Blank-1": [7, 8, 9]}


@pytest.fixture
def codebook_path(tmp_path: Path) -> Path:
    p = tmp_path / "codebook.json"
    p.write_text(json.dumps(CODEBOOK))
    return p


class ImmediateExecutor:
    """Stands in for ProcessPoolExecutor: runs submissions synchronously in-process."""

    def __init__(self, *args, **kwargs): ...

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def submit(self, fn, *args, **kwargs):
        fut: Future = Future()
        try:
            fut.set_result(fn(*args, **kwargs))
        except Exception as e:  # noqa: BLE001
            fut.set_exception(e)
        return fut


class TestWorklist:
    def test_blanks_excluded(self, codebook_path: Path):
        wl = load_worklist(codebook_path)
        assert set(wl) == {"Och.687.1", "Och.958.1"}

    def test_duplicates_rejected(self, tmp_path: Path):
        # json.loads keeps the last value for a repeated key, so a duplicated
        # target would silently design against the wrong bits. Detection has to
        # happen while the document is still a list of pairs.
        p = tmp_path / "dup.json"
        p.write_text('{"A": [1,2,3], "A": [4,5,6], "B": [7,8,9]}')

        with pytest.raises(ValueError, match="more than once: A"):
            load_worklist(p)

    def test_distinct_targets_accepted(self, tmp_path: Path):
        p = tmp_path / "ok.json"
        p.write_text('{"A": [1,2,3], "B": [4,5,6]}')

        assert load_worklist(p) == {"A": [1, 2, 3], "B": [4, 5, 6]}

    def test_final_parquet_naming(self, tmp_path: Path):
        p = final_parquet(tmp_path, "Och.687.1", [3, 1, 2], ("BamHI", "KpnI"))
        assert p.name == "Och.687.1_final_BamHIKpnI_1,2,3.parquet"

    def test_find_missing_final(self, tmp_path: Path):
        wl = {"A": [1, 2, 3], "B": [4, 5, 6]}
        final_parquet(tmp_path, "A", [1, 2, 3], DEFAULT_RESTRICTION).touch()
        assert find_missing_final(wl, tmp_path, DEFAULT_RESTRICTION) == ["B"]

    def test_acceptable_convention_and_override(self, tmp_path: Path, codebook_path: Path):
        assert load_acceptable(codebook_path, None) == {}
        codebook_path.with_suffix(".acceptable.json").write_text('{"Och.687.1": ["Och.1.1"]}')
        assert load_acceptable(codebook_path, None) == {"Och.687.1": ["Och.1.1"]}
        override = tmp_path / "other.json"
        override.write_text('{"Och.958.1": ["Och.2.1"]}')
        assert load_acceptable(codebook_path, override) == {"Och.958.1": ["Och.2.1"]}


class TestRunPanel:
    def _run(self, tmp_path: Path, codebook_path: Path, worker, **kwargs):
        out = tmp_path / "output"
        with (
            patch("mkprobes.run_panel.ProcessPoolExecutor", ImmediateExecutor),
            patch("mkprobes.run_panel.run_gene", side_effect=worker) as mock_worker,
        ):
            summary = run_panel(tmp_path, codebook_path, out, **kwargs)
        return summary, mock_worker, out

    def test_all_genes_submitted(self, tmp_path: Path, codebook_path: Path):
        summary, worker, _ = self._run(tmp_path, codebook_path, lambda *a, **k: None)
        assert sorted(summary["done"]) == ["Och.687.1", "Och.958.1"]
        assert worker.call_count == 2

    def test_finished_genes_skipped(self, tmp_path: Path, codebook_path: Path):
        out = tmp_path / "output"
        out.mkdir()
        final_parquet(out, "Och.687.1", [1, 2, 3], DEFAULT_RESTRICTION).touch()
        summary, worker, _ = self._run(tmp_path, codebook_path, lambda *a, **k: None)
        assert summary["skipped"] == ["Och.687.1"]
        assert summary["done"] == ["Och.958.1"]
        assert worker.call_count == 1

    def test_acceptable_forces_rerun(self, tmp_path: Path, codebook_path: Path):
        out = tmp_path / "output"
        out.mkdir()
        final_parquet(out, "Och.687.1", [1, 2, 3], DEFAULT_RESTRICTION).touch()
        codebook_path.with_suffix(".acceptable.json").write_text('{"Och.687.1": ["Och.1.1"]}')
        summary, worker, _ = self._run(tmp_path, codebook_path, lambda *a, **k: None)
        assert "Och.687.1" in summary["done"]  # finished output, but allow-list forces redo
        forced = {c.kwargs["gene"]: c.kwargs["overwrite"] for c in worker.call_args_list}
        assert forced["Och.687.1"] is True
        assert forced["Och.958.1"] is False

    def test_failure_isolated_and_recorded(self, tmp_path: Path, codebook_path: Path):
        def worker(*a, **k):
            if k["gene"] == "Och.687.1":
                raise Exception("Och.687.1")

        summary, _, _ = self._run(tmp_path, codebook_path, worker)
        assert summary["failed"] == ["Och.687.1"]
        assert summary["done"] == ["Och.958.1"]
        failed_file = codebook_path.parent / "codebook.failed.txt"
        assert failed_file.read_text().strip() == "Och.687.1"

    def test_single_gene_mode_forces_overwrite(self, tmp_path: Path, codebook_path: Path):
        out = tmp_path / "output"
        out.mkdir()
        final_parquet(out, "Och.687.1", [1, 2, 3], DEFAULT_RESTRICTION).touch()
        summary, worker, _ = self._run(tmp_path, codebook_path, lambda *a, **k: None, gene="Och.687.1")
        assert summary["done"] == ["Och.687.1"]
        assert worker.call_args.kwargs["overwrite"] is True

    def test_unknown_single_gene_rejected(self, tmp_path: Path, codebook_path: Path):
        with pytest.raises(ValueError, match="not in the codebook"):
            self._run(tmp_path, codebook_path, lambda *a, **k: None, gene="Nope")


class TestRunPanelCli:
    def test_list_failed(self, tmp_path: Path, codebook_path: Path):
        out = tmp_path / "output"
        out.mkdir()
        final_parquet(out, "Och.958.1", [4, 5, 6], DEFAULT_RESTRICTION).touch()
        res = CliRunner().invoke(
            run_panel_cli, [str(tmp_path), str(codebook_path), "--list-failed", "-o", str(out)]
        )
        assert res.exit_code == 0, res.output
        assert res.output.split() == ["Och.687.1"]

    def test_exit_code_on_failure(self, tmp_path: Path, codebook_path: Path):
        def worker(*a, **k):
            raise Exception(k["gene"])

        with (
            patch("mkprobes.run_panel.ProcessPoolExecutor", ImmediateExecutor),
            patch("mkprobes.run_panel.run_gene", side_effect=worker),
        ):
            res = CliRunner().invoke(run_panel_cli, [str(tmp_path), str(codebook_path)])
        assert res.exit_code == 1
