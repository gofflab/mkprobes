"""
Failure modes a bench biologist will actually hit.

The tool runs for a long time before it fails, and its users do not read
Python tracebacks. These pin the behaviour that a failure is one readable
line naming the fix, and that required arguments are refused up front rather
than surfacing as an AttributeError deep inside the run.
"""

import json
from pathlib import Path

import polars as pl
import pytest
from click.testing import CliRunner

from conftest import flatten_cli_output

from mkprobes import cli
from mkprobes.codebook.finalconstruct import count_final_probes


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


class TestFriendlyErrors:
    def test_unexpected_failure_is_one_readable_message(self, runner: CliRunner, tmp_path: Path):
        bad = tmp_path / "codebook.json"
        bad.write_text("{not valid json")

        result = runner.invoke(cli.main, ["hash", str(bad)])

        assert result.exit_code != 0
        assert "JSONDecodeError" in flatten_cli_output(result.output)
        assert "--debug" in flatten_cli_output(result.output)
        assert "Traceback" not in flatten_cli_output(result.output)

    def test_debug_flag_reraises(self, runner: CliRunner, tmp_path: Path):
        bad = tmp_path / "codebook.json"
        bad.write_text("{not valid json")

        result = runner.invoke(cli.main, ["--debug", "hash", str(bad)])

        assert isinstance(result.exception, json.JSONDecodeError)

    def test_usage_errors_are_untouched(self, runner: CliRunner):
        result = runner.invoke(cli.main, ["hash"])

        assert result.exit_code != 0
        assert "Missing argument" in flatten_cli_output(result.output)


class TestRequiredOptions:
    """These used to fail with `AttributeError: 'NoneType' has no attribute 'read_text'`."""

    def test_filter_genes_requires_genes(self, runner: CliRunner, tmp_path: Path):
        result = runner.invoke(cli.main, ["filter-genes", str(tmp_path)])

        assert result.exit_code != 0
        assert "--genes" in flatten_cli_output(result.output)
        assert "AttributeError" not in flatten_cli_output(result.output)

    def test_construct_requires_target(self, runner: CliRunner, tmp_path: Path):
        result = runner.invoke(cli.main, ["construct", str(tmp_path), str(tmp_path)])

        assert result.exit_code != 0
        assert "--gene" in flatten_cli_output(result.output)
        assert "AttributeError" not in flatten_cli_output(result.output)

    def test_construct_requires_codebook(self, runner: CliRunner, tmp_path: Path):
        # Click reports only the first missing option, so supply --gene to reach it.
        result = runner.invoke(
            cli.main, ["construct", str(tmp_path), str(tmp_path), "--gene", "Sox2"]
        )

        assert result.exit_code != 0
        assert "--codebook" in flatten_cli_output(result.output)
        assert "AttributeError" not in flatten_cli_output(result.output)


class TestFilterGenesCountsFinalProbes:
    def _final(self, output: Path, gene: str, n: int) -> None:
        output.mkdir(parents=True, exist_ok=True)
        pl.DataFrame({"seq": ["ACGT"] * n}).write_parquet(
            output / f"{gene}_final_BamHIKpnI_1,2,3.parquet"
        )

    def _screened(self, output: Path, gene: str, n: int) -> None:
        output.mkdir(parents=True, exist_ok=True)
        pl.DataFrame({"seq": ["ACGT"] * n}).write_parquet(
            output / f"{gene}_screened_ol-2_BamHIKpnI.parquet"
        )

    def test_counts_final_not_screened(self, tmp_path: Path):
        # Screening yields the candidate pool; construct decides how many
        # survive readout attachment. Counting the former over-reports.
        self._screened(tmp_path, "GeneA", 90)
        self._final(tmp_path, "GeneA", 12)

        assert count_final_probes(tmp_path, "GeneA") == 12

    def test_missing_target_reads_as_none(self, tmp_path: Path):
        tmp_path.joinpath("output").mkdir()

        assert count_final_probes(tmp_path / "output", "Absent") is None

    def test_pass_list_uses_final_counts(self, runner: CliRunner, tmp_path: Path):
        output = tmp_path / "output"
        self._screened(output, "Rich", 90)
        self._final(output, "Rich", 60)
        self._screened(output, "Thin", 90)
        self._final(output, "Thin", 10)
        genes = tmp_path / "genes.txt"
        genes.write_text("Rich\nThin\n")
        out = tmp_path / "pass.txt"

        result = runner.invoke(
            cli.main,
            ["filter-genes", str(output), "--genes", str(genes), "--min-probes", "48",
             "--out", str(out)],
        )

        assert result.exit_code == 0, result.output
        assert out.read_text().split() == ["Rich"]

    def test_missing_targets_are_reported(self, runner: CliRunner, tmp_path: Path):
        output = tmp_path / "output"
        self._final(output, "Present", 60)
        genes = tmp_path / "genes.txt"
        genes.write_text("Present\nNeverRan\n")

        result = runner.invoke(cli.main, ["filter-genes", str(output), "--genes", str(genes)])

        assert result.exit_code == 0, result.output
        # The old code raised FileNotFoundError naming a path it never writes.
        assert "FileNotFoundError" not in flatten_cli_output(result.output)


class TestRestrictionIsFixedByChemistry:
    """
    `--restriction` looked configurable but assembly only ever looked for the
    default pair's filenames, so a different pair computed a whole panel and
    then failed to find a single gene.
    """

    @pytest.mark.parametrize(
        "command,args",
        [
            ("run-panel", ["--restriction", "EcoRI,XhoI"]),
            ("screen", ["--restriction", "EcoRI,XhoI"]),
        ],
    )
    def test_other_enzymes_are_refused_before_any_work(
        self, runner: CliRunner, tmp_path: Path, command: str, args: list[str]
    ):
        codebook = tmp_path / "codebook.json"
        codebook.write_text(json.dumps({"A": [1, 2, 3]}))
        base = (
            [command, str(tmp_path), str(codebook)]
            if command == "run-panel"
            else [command, str(tmp_path), "A"]
        )

        result = runner.invoke(cli.main, base + args)

        assert result.exit_code != 0
        assert "--restriction" in flatten_cli_output(result.output)
        assert "BamHI" in flatten_cli_output(result.output)

    def test_the_default_pair_is_accepted(self, runner: CliRunner, tmp_path: Path):
        codebook = tmp_path / "codebook.json"
        codebook.write_text(json.dumps({"A": [1, 2, 3]}))

        result = runner.invoke(
            cli.main,
            ["run-panel", str(tmp_path), str(codebook), "--restriction", "KpnI,BamHI",
             "--list-failed"],
        )

        # Order is not significant, and --list-failed exits before any design work.
        assert result.exit_code == 0, result.output

    def test_assembly_filenames_track_the_constant(self):
        # The token was hardcoded in three places; assembly must derive it, or
        # the two drift apart again.
        import importlib

        source = importlib.import_module("mkprobes.assembly").__loader__.get_source(
            "mkprobes.assembly"
        )
        assert "_final_BamHIKpnI_" not in source
        assert "RESTRICTION_TOKEN" in source


class TestHelpAlwaysSucceeds:
    """
    The friendly error handler wraps unexpected exceptions, and click signals
    `--help` by raising `click.exceptions.Exit` - which is a RuntimeError. An
    over-broad catch turned every subcommand's help into exit code 1.
    """

    @pytest.mark.parametrize("command", sorted(cli.main.commands))
    def test_subcommand_help_exits_zero(self, runner: CliRunner, command: str):
        result = runner.invoke(cli.main, [command, "--help"])

        assert result.exit_code == 0, f"`mkprobes {command} --help` failed:\n{result.output}"
        assert "Usage:" in flatten_cli_output(result.output)

    def test_group_help_exits_zero(self, runner: CliRunner):
        result = runner.invoke(cli.main, ["--help"])

        assert result.exit_code == 0


class TestExternalToolPreflight:
    """
    A missing aligner used to surface as a raw subprocess error inside sixteen
    worker processes, minutes into a run. These pin that it is caught up front,
    and that read-only paths are not gated on it.
    """

    def test_candidates_names_the_missing_tool_and_the_fix(
        self, runner: CliRunner, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        monkeypatch.setattr("shutil.which", lambda _: None)

        result = runner.invoke(cli.main, ["candidates", str(tmp_path), "--gene", "Sox2"])

        assert result.exit_code != 0
        assert "bowtie2" in flatten_cli_output(result.output)
        assert "conda install" in flatten_cli_output(result.output)

    def test_run_panel_checks_before_designing(
        self, runner: CliRunner, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        codebook = tmp_path / "codebook.json"
        codebook.write_text(json.dumps({"A": [1, 2, 3]}))
        monkeypatch.setattr("shutil.which", lambda _: None)

        result = runner.invoke(cli.main, ["run-panel", str(tmp_path), str(codebook)])

        assert result.exit_code != 0
        assert "bowtie2" in flatten_cli_output(result.output)

    def test_list_failed_does_not_need_the_aligner(
        self, runner: CliRunner, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        # Inspecting a finished panel is read-only; gating it on bowtie2 would
        # be gratuitous.
        codebook = tmp_path / "codebook.json"
        codebook.write_text(json.dumps({"A": [1, 2, 3]}))
        monkeypatch.setattr("shutil.which", lambda _: None)

        result = runner.invoke(
            cli.main, ["run-panel", str(tmp_path), str(codebook), "--list-failed"]
        )

        assert result.exit_code == 0, result.output
        assert "A" in flatten_cli_output(result.output)
