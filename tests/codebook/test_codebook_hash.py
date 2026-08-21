"""
The codebook hash as a durable identifier.

It answers "which codebook produced this file?" months later, so it has to be
one value everywhere and it has to outlive the terminal it was printed in.
"""

import json
from pathlib import Path

import polars as pl
import pytest
from click.testing import CliRunner

from mkprobes import cli
from mkprobes.codebook.codebook import ProbeSet, hash_codebook, hash_codebook_file
from mkprobes.utils.provenance import read_provenance

CODEBOOK = {"Sox2": [1, 2, 3], "Pax6": [4, 5, 6], "Blank-1": [7, 8, 9], "Blank-2": [1, 4, 7]}


@pytest.fixture
def codebook_file(tmp_path: Path) -> Path:
    path = tmp_path / "codebook.json"
    path.write_text(json.dumps(CODEBOOK, indent=2))
    return path


class TestOneValueEverywhere:
    def test_file_hash_covers_blank_codes(self, codebook_file: Path):
        # make-codebook hashes the codebook as written; downstream loaders drop
        # the Blanks. Hashing whichever dict is in hand gave two answers.
        without_blanks = {k: v for k, v in CODEBOOK.items() if not k.startswith("Blank")}

        assert hash_codebook_file(codebook_file) == hash_codebook(CODEBOOK)
        assert hash_codebook_file(codebook_file) != hash_codebook(without_blanks)

    def test_probeset_resolves_the_same_file(self, codebook_file: Path):
        probeset = ProbeSet(
            name="p", species="mouse", codebook="codebook.json", bcidx=0, n_probes=4
        )

        assert probeset.codebook_path(codebook_file.parent) == codebook_file
        assert hash_codebook_file(probeset.codebook_path(codebook_file.parent)) == hash_codebook(
            CODEBOOK
        )

    def test_reordering_targets_does_not_change_it(self, tmp_path: Path):
        # Hashing is key-sorted, so a re-serialised codebook is the same codebook.
        a = tmp_path / "a.json"
        b = tmp_path / "b.json"
        a.write_text(json.dumps(CODEBOOK))
        b.write_text(json.dumps(dict(reversed(list(CODEBOOK.items()))), indent=4))

        assert hash_codebook_file(a) == hash_codebook_file(b)

    def test_changing_bits_changes_it(self, tmp_path: Path):
        a = tmp_path / "a.json"
        b = tmp_path / "b.json"
        a.write_text(json.dumps(CODEBOOK))
        b.write_text(json.dumps({**CODEBOOK, "Sox2": [1, 2, 9]}))

        assert hash_codebook_file(a) != hash_codebook_file(b)


class TestSurvivesTheTerminal:
    def test_make_codebook_writes_a_sidecar(self, tmp_path: Path):
        genes = tmp_path / "genes.txt"
        genes.write_text("Sox2\nPax6\nEomes\n")
        out = tmp_path / "cb.json"

        result = CliRunner().invoke(cli.main, ["make-codebook", str(tmp_path), str(genes), "-o", str(out)])

        assert result.exit_code == 0, result.output
        sidecar = out.with_suffix(".hash")
        assert sidecar.exists()
        assert sidecar.read_text().strip() == hash_codebook_file(out)

    def test_sidecar_agrees_with_the_hash_command(self, tmp_path: Path):
        genes = tmp_path / "genes.txt"
        genes.write_text("Sox2\nPax6\n")
        out = tmp_path / "cb.json"
        CliRunner().invoke(cli.main, ["make-codebook", str(tmp_path), str(genes), "-o", str(out)])

        printed = CliRunner().invoke(cli.main, ["hash", str(out)])

        assert printed.output.strip() == out.with_suffix(".hash").read_text().strip()

    def test_sidecar_is_not_mistaken_for_a_target(self, tmp_path: Path):
        # A `_hash` key inside the JSON would be read as an extra target by
        # every consumer that iterates it as {target: bits}.
        genes = tmp_path / "genes.txt"
        genes.write_text("Sox2\nPax6\n")
        out = tmp_path / "cb.json"
        CliRunner().invoke(cli.main, ["make-codebook", str(tmp_path), str(genes), "-o", str(out)])

        assert all(isinstance(v, list) for v in json.loads(out.read_text()).values())


class TestRecordedInOutputs:
    def test_construct_output_names_its_codebook(self, tmp_path: Path, monkeypatch):
        from mkprobes.utils.provenance import provenance_metadata

        digest = "abc123"
        path = tmp_path / "out.parquet"
        pl.DataFrame({"seq": ["ACGT"]}).write_parquet(
            path, metadata=provenance_metadata(stage="construct", codebook_hash=digest)
        )

        assert read_provenance(path)["codebook_hash"] == digest
