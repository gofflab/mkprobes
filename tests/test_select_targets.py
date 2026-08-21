"""
Optional target suggestion.

The method picks genes whose expression is independent of the ones already
chosen. That has a sharp failure mode: a gene correlating with nothing looks
maximally independent, so unfiltered expression data yields suggestions that are
purely noise. These tests pin both the useful behaviour and the guard against
that failure.
"""

from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner

from conftest import flatten_cli_output

from mkprobes import cli

pytest.importorskip("anndata")


def synthetic(tmp_path: Path, n_noise: int = 0, name: str = "expr.h5ad", seed: int = 1) -> Path:
    """
    Expression with known structure: 8 latent programmes, 40 genes reading each,
    plus `n_noise` genes that correlate with nothing.
    """
    import anndata as ad

    rng = np.random.default_rng(seed)
    n_obs, n_prog, per_prog = 400, 8, 40
    factors = rng.normal(size=(n_prog, n_obs))
    columns, names = [], []
    for programme in range(n_prog):
        for i in range(per_prog):
            columns.append(
                factors[programme] * rng.uniform(0.6, 1.4) + rng.normal(scale=0.4, size=n_obs)
            )
            names.append(f"P{programme:02d}_{i:03d}")
    for j in range(n_noise):
        columns.append(rng.normal(size=n_obs))
        names.append(f"N_{j:03d}")

    adata = ad.AnnData(np.array(columns).T.astype("float64"))
    adata.var_names = names
    path = tmp_path / name
    adata.write_h5ad(path)
    return path


def programme_of(gene: str) -> int:
    """Which latent programme a synthetic gene reads, or -1 for a noise gene."""
    return -1 if gene.startswith("N_") else int(gene[1:3])


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


@pytest.fixture
def log_messages():
    """Captures loguru output, which goes to stderr rather than into result.output."""
    from loguru import logger

    messages: list[str] = []
    sink = logger.add(messages.append, level="WARNING")
    yield messages
    logger.remove(sink)


class TestSuggestions:
    def test_suggests_genes_from_programmes_not_already_covered(self, runner: CliRunner, tmp_path: Path):
        expr = synthetic(tmp_path)
        have = tmp_path / "have.txt"
        have.write_text("P00_000\nP00_001\nP00_002\n")  # three genes, all one programme

        result = runner.invoke(
            cli.main, ["suggest-targets", str(expr), "--add", "6", "--have", str(have)]
        )

        assert result.exit_code == 0, result.output
        suggested = [line for line in result.output.splitlines() if line.startswith("P")]
        assert len(suggested) == 6
        # The point of the method: spread across programmes rather than piling
        # onto the one already covered.
        assert len({programme_of(g) for g in suggested}) >= 3
        # Each gene carries its own noise, so the covered programme is not
        # driven exactly to zero and may still contribute one pick. It must not
        # dominate, which is what piling onto one axis would look like.
        assert sum(programme_of(g) == 0 for g in suggested) <= 1

    def test_works_with_no_prior_targets(self, runner: CliRunner, tmp_path: Path):
        expr = synthetic(tmp_path)

        result = runner.invoke(cli.main, ["suggest-targets", str(expr), "--add", "5"])

        assert result.exit_code == 0, result.output
        assert len([ln for ln in result.output.splitlines() if ln.startswith("P")]) == 5

    def test_is_reproducible_for_a_given_seed(self, runner: CliRunner, tmp_path: Path):
        expr = synthetic(tmp_path)
        args = ["suggest-targets", str(expr), "--add", "5", "--seed", "7"]

        first = runner.invoke(cli.main, args)
        second = runner.invoke(cli.main, args)

        assert first.exit_code == 0 and second.exit_code == 0
        assert first.output == second.output

    def test_writes_a_list_the_next_step_can_read(self, runner: CliRunner, tmp_path: Path):
        from mkprobes.utils.targets import read_target_list

        expr = synthetic(tmp_path)
        have = tmp_path / "have.txt"
        have.write_text("# my markers\nP00_000\n\nP01_000  # keep\n")
        out = tmp_path / "panel.txt"

        result = runner.invoke(
            cli.main,
            ["suggest-targets", str(expr), "--add", "4", "--have", str(have), "-o", str(out)],
        )

        assert result.exit_code == 0, result.output
        targets = read_target_list(out)
        assert targets[:2] == ["P00_000", "P01_000"]  # prior targets kept, in order
        assert len(targets) == 6


class TestUnfilteredInputGuard:
    def test_warns_when_suggestions_are_noise(
        self, runner: CliRunner, tmp_path: Path, log_messages: list[str]
    ):
        # Unstructured genes dominate the selection: they correlate with
        # nothing, which this method reads as maximal independence.
        expr = synthetic(tmp_path, n_noise=200)

        result = runner.invoke(cli.main, ["suggest-targets", str(expr), "--add", "6"])

        assert result.exit_code == 0, result.output
        # The failure mode itself: every pick is a noise gene.
        picked = [ln for ln in result.output.splitlines() if ln.strip()]
        assert all(programme_of(g) == -1 for g in picked)
        warning = "".join(log_messages)
        assert "uninformative genes" in warning
        assert "highly_variable_genes" in warning

    def test_stays_quiet_on_filtered_input(
        self, runner: CliRunner, tmp_path: Path, log_messages: list[str]
    ):
        expr = synthetic(tmp_path, n_noise=0)

        result = runner.invoke(cli.main, ["suggest-targets", str(expr), "--add", "6"])

        assert result.exit_code == 0, result.output
        assert "uninformative genes" not in "".join(log_messages)


class TestInputValidation:
    def test_rejects_unknown_prior_targets(self, runner: CliRunner, tmp_path: Path):
        expr = synthetic(tmp_path)
        have = tmp_path / "have.txt"
        have.write_text("P00_000\nNotAGene\n")

        result = runner.invoke(
            cli.main, ["suggest-targets", str(expr), "--add", "3", "--have", str(have)]
        )

        assert result.exit_code != 0
        assert "NotAGene" in flatten_cli_output(result.output)

    def test_rejects_asking_for_more_genes_than_exist(self, runner: CliRunner, tmp_path: Path):
        expr = synthetic(tmp_path)

        result = runner.invoke(cli.main, ["suggest-targets", str(expr), "--add", "99999"])

        assert result.exit_code != 0
        assert "Lower --add" in flatten_cli_output(result.output)

    def test_rejects_a_missing_layer(self, runner: CliRunner, tmp_path: Path):
        expr = synthetic(tmp_path)

        result = runner.invoke(
            cli.main, ["suggest-targets", str(expr), "--add", "3", "--layer", "lognorm"]
        )

        assert result.exit_code != 0
        assert "no layer" in flatten_cli_output(result.output)

    def test_rejects_a_file_that_is_not_anndata(self, runner: CliRunner, tmp_path: Path):
        not_h5ad = tmp_path / "expr.h5ad"
        not_h5ad.write_text("this is not an h5ad file")

        result = runner.invoke(cli.main, ["suggest-targets", str(not_h5ad), "--add", "3"])

        assert result.exit_code != 0
        assert "AnnData" in flatten_cli_output(result.output)
