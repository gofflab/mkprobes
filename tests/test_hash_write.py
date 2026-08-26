"""
`mkprobes hash --write` backfills the sidecar.

make-codebook writes <codebook>.hash, but only for codebooks it creates. A
codebook made before that existed has no sidecar, and re-running make-codebook
to get one would reassign every bit - so the hash has to be writable from the
codebook already on disk.
"""

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from conftest import flatten_cli_output

from mkprobes import cli
from mkprobes.codebook.codebook import hash_codebook_file

CODEBOOK = {"Sox2": [1, 2, 3], "Pax6": [4, 5, 6], "Blank-1": [7, 8, 9]}


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


@pytest.fixture
def codebook(tmp_path: Path) -> Path:
    path = tmp_path / "codebook.json"
    path.write_text(json.dumps(CODEBOOK, indent=2))
    return path


@pytest.fixture
def warnings():
    from loguru import logger

    messages: list[str] = []
    sink = logger.add(messages.append, level="WARNING")
    yield messages
    logger.remove(sink)


def test_prints_without_writing_by_default(runner: CliRunner, codebook: Path):
    result = runner.invoke(cli.main, ["hash", str(codebook)])

    assert result.exit_code == 0, result.output
    assert result.output.strip() == hash_codebook_file(codebook)
    assert not codebook.with_suffix(".hash").exists()


def test_write_creates_the_sidecar(runner: CliRunner, codebook: Path):
    result = runner.invoke(cli.main, ["hash", str(codebook), "--write"])

    assert result.exit_code == 0, result.output
    assert codebook.with_suffix(".hash").read_text().strip() == hash_codebook_file(codebook)


def test_the_codebook_is_untouched(runner: CliRunner, codebook: Path):
    # Backfilling must not reassign anything; that is the whole point.
    before = codebook.read_text()

    runner.invoke(cli.main, ["hash", str(codebook), "--write"])

    assert codebook.read_text() == before


def test_matches_what_make_codebook_would_have_written(runner: CliRunner, tmp_path: Path):
    genes = tmp_path / "genes.txt"
    genes.write_text("Sox2\nPax6\n")
    generated = tmp_path / "cb.json"
    runner.invoke(cli.main, ["make-codebook", str(tmp_path), str(genes), "-o", str(generated)])
    from_make = generated.with_suffix(".hash").read_text()

    generated.with_suffix(".hash").unlink()
    runner.invoke(cli.main, ["hash", str(generated), "--write"])

    assert generated.with_suffix(".hash").read_text() == from_make


def test_rewriting_is_idempotent(runner: CliRunner, codebook: Path):
    runner.invoke(cli.main, ["hash", str(codebook), "--write"])
    first = codebook.with_suffix(".hash").read_text()

    runner.invoke(cli.main, ["hash", str(codebook), "--write"])

    assert codebook.with_suffix(".hash").read_text() == first


def test_warns_when_the_codebook_changed_under_an_existing_sidecar(
    runner: CliRunner, codebook: Path, warnings: list[str]
):
    runner.invoke(cli.main, ["hash", str(codebook), "--write"])
    stale = codebook.with_suffix(".hash").read_text().strip()
    codebook.write_text(json.dumps({**CODEBOOK, "Sox2": [1, 2, 9]}, indent=2))

    result = runner.invoke(cli.main, ["hash", str(codebook), "--write"])

    assert result.exit_code == 0, result.output
    assert stale in "".join(warnings)
    # Updated regardless: the sidecar identifies the codebook beside it.
    assert codebook.with_suffix(".hash").read_text().strip() == hash_codebook_file(codebook)


def test_missing_codebook_is_refused(runner: CliRunner, tmp_path: Path):
    result = runner.invoke(cli.main, ["hash", str(tmp_path / "absent.json"), "--write"])

    assert result.exit_code != 0
    assert "absent.json" in flatten_cli_output(result.output)
