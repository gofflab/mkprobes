"""Provenance stamping and read-back."""

import json
from pathlib import Path

import polars as pl
import pytest
from click.testing import CliRunner

from mkprobes import cli
from mkprobes.utils.provenance import (
    PROVENANCE_KEY,
    encode,
    provenance_metadata,
    provenance_record,
    read_provenance,
)


class TestProvenanceRecord:
    def test_carries_version_time_and_command(self):
        record = provenance_record()

        assert record["mkprobes_version"]
        assert record["command"]
        # Timestamps must be UTC-qualified, or they are useless across machines.
        assert record["generated_at"].endswith("+00:00")

    def test_dataset_path_is_absolute(self, tmp_path: Path):
        record = provenance_record(tmp_path)

        assert Path(record["dataset"]).is_absolute()

    def test_parameters_are_recorded(self):
        record = provenance_record(stage="screen", overlap=-2, restriction=["BamHI", "KpnI"])

        assert record["stage"] == "screen"
        assert record["overlap"] == -2
        assert record["restriction"] == ["BamHI", "KpnI"]

    def test_none_parameters_are_dropped(self):
        # Callers pass optional arguments straight through.
        assert "fpkm_path" not in provenance_record(stage="screen", fpkm_path=None)

    def test_encodes_to_string_valued_metadata(self):
        metadata = provenance_metadata(stage="construct")

        assert set(metadata) == {PROVENANCE_KEY}
        assert isinstance(metadata[PROVENANCE_KEY], str)
        assert json.loads(metadata[PROVENANCE_KEY])["stage"] == "construct"


class TestRoundTrip:
    def test_survives_a_parquet_write(self, tmp_path: Path):
        path = tmp_path / "out.parquet"
        record = provenance_record(tmp_path, stage="screen", overlap=-2)
        pl.DataFrame({"seq": ["ACGT"]}).write_parquet(path, metadata=encode(record))

        assert read_provenance(path) == record

    def test_missing_provenance_reads_as_none(self, tmp_path: Path):
        path = tmp_path / "bare.parquet"
        pl.DataFrame({"seq": ["ACGT"]}).write_parquet(path)

        assert read_provenance(path) is None

    def test_unreadable_file_reads_as_none(self, tmp_path: Path):
        path = tmp_path / "not.parquet"
        path.write_text("this is not a parquet file")

        assert read_provenance(path) is None


class TestProvenanceCommand:
    def test_prints_the_record(self, tmp_path: Path):
        path = tmp_path / "out.parquet"
        pl.DataFrame({"seq": ["ACGT"]}).write_parquet(
            path, metadata=provenance_metadata(stage="screen")
        )

        result = CliRunner().invoke(cli.main, ["provenance", str(path)])

        assert result.exit_code == 0
        assert json.loads(result.output)["stage"] == "screen"

    def test_explains_when_absent(self, tmp_path: Path):
        path = tmp_path / "bare.parquet"
        pl.DataFrame({"seq": ["ACGT"]}).write_parquet(path)

        result = CliRunner().invoke(cli.main, ["provenance", str(path)])

        assert result.exit_code != 0
        assert "carries no provenance" in result.output


@pytest.mark.parametrize(
    "module", ["screen", "candidates", "codebook.finalconstruct", "assembly"]
)
def test_every_pipeline_parquet_write_is_stamped(module: str):
    """
    Guards the invariant rather than the current call sites: a new
    `write_parquet` in a pipeline module must pass `metadata=`, or a future
    output silently loses its provenance.
    """
    import ast
    import importlib

    source = importlib.import_module(f"mkprobes.{module}").__loader__.get_source(f"mkprobes.{module}")
    unstamped = [
        node.lineno
        for node in ast.walk(ast.parse(source))
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == "write_parquet"
        and not any(keyword.arg == "metadata" for keyword in node.keywords)
    ]
    assert not unstamped, f"mkprobes/{module}.py write_parquet without metadata= at lines {unstamped}"
