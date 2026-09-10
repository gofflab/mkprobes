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
from mkprobes.codebook.generate import ORDER, make_codebook
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

    def test_stub_carries_the_offset(self, project: Path):
        entries = json.loads((project / "manifest.json").read_text())
        assert entries[0]["offset"] == 0
        assert "offset" in entries[0]["_comment"]

    def test_offset_is_written_from_the_flag(self, runner: CliRunner, tmp_path: Path):
        result = runner.invoke(cli.main, ["init", str(tmp_path / "panel_b"), "--offset", "10"])

        assert result.exit_code == 0, result.output
        assert json.loads((tmp_path / "panel_b" / "manifest.json").read_text())[0]["offset"] == 10

    def test_out_of_range_offset_is_refused(self, runner: CliRunner, tmp_path: Path):
        result = runner.invoke(cli.main, ["init", str(tmp_path / "p"), "--offset", str(len(ORDER))])

        assert result.exit_code != 0
        assert "--offset" in flatten_cli_output(result.output)


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

    def test_rejects_out_of_range_offset(self, tmp_path: Path):
        entries = manifest_stub("p", "mouse")
        entries[0]["offset"] = len(ORDER)
        path = self._write(tmp_path, entries)

        with pytest.raises(ValueError, match="offset"):
            check_manifest(path)

    def test_rejects_an_unknown_field(self, tmp_path: Path):
        # A misspelt offset used to be ignored, which for pooled panels means
        # two codebooks silently sharing bits.
        entries = manifest_stub("p", "mouse")
        entries[0]["offest"] = 10
        path = self._write(tmp_path, entries)

        with pytest.raises(ValueError, match="offest"):
            check_manifest(path)

    def test_rejects_a_codebook_that_does_not_start_at_the_offset(self, tmp_path: Path):
        # The manifest was edited after the codebook was generated (or the
        # codebook was generated with --offset overriding it).
        entries = manifest_stub("p", "mouse")
        entries[0]["offset"] = 10
        path = self._write(tmp_path, entries)  # codebook starts at bit 1, position 0

        with pytest.raises(ValueError, match="starts at bit position 0"):
            check_manifest(path)

    def test_accepts_a_codebook_generated_at_the_offset(self, tmp_path: Path):
        entries = manifest_stub("p", "mouse")
        entries[0]["offset"] = 10
        path = self._write(tmp_path, entries)
        (tmp_path / "codebook.json").write_text(json.dumps(make_codebook(["A", "B"], n_bits=10, offset=10)))

        assert check_manifest(path)[0].offset == 10

    def test_rejects_bits_outside_the_readout_table(self, tmp_path: Path):
        path = self._write(tmp_path, manifest_stub("p", "mouse"))
        (tmp_path / "codebook.json").write_text(json.dumps({"A": [1, 2, len(ORDER) + 1]}))

        with pytest.raises(ValueError, match="only 1 to"):
            check_manifest(path)

    def test_rejects_pooled_panels_sharing_bits(self, tmp_path: Path):
        # One manifest, two panels, both generated at offset 0: the very
        # mistake `offset` exists to prevent.
        a, b = manifest_stub("a", "mouse")[0], manifest_stub("b", "mouse", bcidx=1)[0]
        b["codebook"] = "codebook_b.json"
        path = self._write(tmp_path, [a, b])
        (tmp_path / "codebook.json").write_text(json.dumps(make_codebook(["A"], n_bits=10)))
        (tmp_path / "codebook_b.json").write_text(json.dumps(make_codebook(["B"], n_bits=10)))

        with pytest.raises(ValueError, match="share readout bit"):
            check_manifest(path)

    def test_accepts_pooled_panels_with_disjoint_bits(self, tmp_path: Path):
        a, b = manifest_stub("a", "mouse")[0], manifest_stub("b", "mouse", bcidx=1, offset=10)[0]
        b["codebook"] = "codebook_b.json"
        path = self._write(tmp_path, [a, b])
        (tmp_path / "codebook.json").write_text(json.dumps(make_codebook(["A"], n_bits=10)))
        (tmp_path / "codebook_b.json").write_text(json.dumps(make_codebook(["B"], n_bits=10, offset=10)))

        assert [p.offset for p in check_manifest(path)] == [0, 10]

    def test_one_codebook_under_two_probe_sets_is_not_a_clash(self, tmp_path: Path):
        # run-panel allows the same codebook under several probe sets; that is
        # one panel described twice, not two panels colliding.
        a, b = manifest_stub("a", "mouse")[0], manifest_stub("b", "mouse", bcidx=1)[0]
        path = self._write(tmp_path, [a, b])

        assert len(check_manifest(path)) == 2

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
        assert "offset 0" in flatten_cli_output(result.output)


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
