"""
Project scaffolding and manifest validation.

The manifest was the one file with no generator, no schema documentation and no
validation, sitting at the very end of a workflow that takes hours. These pin
that a scaffolded project is valid on the first try, and that the ways it can be
wrong are caught before assembly starts.
"""

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from conftest import flatten_cli_output

from mkprobes import cli
from mkprobes.init_project import check_manifest, manifest_stub, max_bcidx
from mkprobes.utils.targets import read_target_list


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


@pytest.fixture
def project(runner: CliRunner, tmp_path: Path) -> Path:
    target = tmp_path / "panel_a"
    result = runner.invoke(cli.main, ["init", str(target), "--species", "mouse"])
    assert result.exit_code == 0, result.output
    return target


class TestInit:
    def test_writes_the_three_files(self, project: Path):
        assert (project / "genes.txt").exists()
        assert (project / "manifest.json").exists()
        assert (project / "README.md").exists()

    def test_target_template_is_readable_by_the_pipeline(self, project: Path):
        # The template is mostly comments; the readers must cope with that.
        assert read_target_list(project / "genes.txt") == ["Sox2", "Pax6"]

    def test_manifest_validates_once_a_codebook_exists(self, project: Path):
        (project / "codebook.json").write_text(json.dumps({"Sox2": [1, 2, 3]}))

        probesets = check_manifest(project / "manifest.json")

        assert [p.name for p in probesets] == ["panel_a"]
        assert probesets[0].species == "mouse"

    def test_refuses_to_clobber(self, runner: CliRunner, project: Path):
        result = runner.invoke(cli.main, ["init", str(project)])

        assert result.exit_code != 0
        assert "--force" in flatten_cli_output(result.output)

    def test_force_overwrites(self, runner: CliRunner, project: Path):
        (project / "genes.txt").write_text("Edited\n")

        result = runner.invoke(cli.main, ["init", str(project), "--force"])

        assert result.exit_code == 0, result.output
        assert "Sox2" in (project / "genes.txt").read_text()

    def test_out_of_range_bcidx_is_refused(self, runner: CliRunner, tmp_path: Path):
        result = runner.invoke(
            cli.main, ["init", str(tmp_path / "p"), "--bcidx", str(max_bcidx() + 1)]
        )

        assert result.exit_code != 0
        assert "--bcidx" in flatten_cli_output(result.output)


class TestCheckManifest:
    def _write(self, tmp_path: Path, entries: object) -> Path:
        path = tmp_path / "manifest.json"
        path.write_text(json.dumps(entries))
        (tmp_path / "codebook.json").write_text(json.dumps({"A": [1, 2, 3]}))
        return path

    def test_rejects_out_of_range_bcidx(self, tmp_path: Path):
        # Out of range used to crash at the very last step, after the whole panel.
        entries = manifest_stub("p", "mouse")
        entries[0]["bcidx"] = max_bcidx() + 5
        path = self._write(tmp_path, entries)

        with pytest.raises(ValueError, match="bcidx"):
            check_manifest(path)

    def test_rejects_a_missing_codebook(self, tmp_path: Path):
        path = tmp_path / "manifest.json"
        path.write_text(json.dumps(manifest_stub("p", "mouse")))

        with pytest.raises(ValueError, match="make-codebook"):
            check_manifest(path)

    def test_rejects_duplicate_names(self, tmp_path: Path):
        entries = manifest_stub("p", "mouse") + manifest_stub("p", "mouse")
        path = self._write(tmp_path, entries)

        with pytest.raises(ValueError, match="unique"):
            check_manifest(path)

    def test_rejects_an_empty_manifest(self, tmp_path: Path):
        path = self._write(tmp_path, [])

        with pytest.raises(ValueError, match="no probe sets"):
            check_manifest(path)

    def test_rejects_malformed_json(self, tmp_path: Path):
        path = tmp_path / "manifest.json"
        path.write_text('[{"name": "p"}]')  # missing required fields

        with pytest.raises(ValueError, match="not a valid manifest"):
            check_manifest(path)

    def test_cli_reports_success(self, runner: CliRunner, tmp_path: Path):
        path = self._write(tmp_path, manifest_stub("p", "mouse"))

        result = runner.invoke(cli.main, ["check-manifest", str(path)])

        assert result.exit_code == 0, result.output
        assert "is valid" in flatten_cli_output(result.output)


class TestReadTargetList:
    def test_skips_comments_and_blanks(self, tmp_path: Path):
        path = tmp_path / "genes.txt"
        path.write_text("# header\n\nSox2\nPax6  # marker\n\n")

        assert read_target_list(path) == ["Sox2", "Pax6"]

    def test_preserves_order(self, tmp_path: Path):
        path = tmp_path / "genes.txt"
        path.write_text("Zic1\nAldoc\nSox2\n")

        assert read_target_list(path) == ["Zic1", "Aldoc", "Sox2"]

    def test_rejects_duplicates_by_name(self, tmp_path: Path):
        path = tmp_path / "genes.txt"
        path.write_text("Sox2\nPax6\nSox2\n")

        with pytest.raises(ValueError, match="Sox2"):
            read_target_list(path)

    def test_rejects_an_empty_list(self, tmp_path: Path):
        path = tmp_path / "genes.txt"
        path.write_text("# only comments\n\n")

        with pytest.raises(ValueError, match="no targets"):
            read_target_list(path)


class TestGeneratedReadme:
    """The README is the walkthrough a new user follows, so its commands have to
    be real and appropriate to the species."""

    @pytest.mark.parametrize(
        "species,expected",
        [("mouse", False), ("human", False), ("octopus", True)],
    )
    def test_transcript_mode_matches_the_dataset_kind(
        self, runner: CliRunner, tmp_path: Path, species: str, expected: bool
    ):
        # Reference datasets pick the canonical isoform from Ensembl; only
        # custom datasets need the longest-transcript fallback.
        target = tmp_path / species
        assert runner.invoke(cli.main, ["init", str(target), "--species", species]).exit_code == 0

        readme = target.joinpath("README.md").read_text()
        assert ("-m longest" in readme) is expected

    def test_every_command_named_exists(self, project: Path):
        import re

        named = set(re.findall(r"^mkprobes ([a-z][a-z-]+)", project.joinpath("README.md").read_text(), re.MULTILINE))
        assert named, "README names no commands"
        assert named <= set(cli.main.commands), f"unknown: {named - set(cli.main.commands)}"
