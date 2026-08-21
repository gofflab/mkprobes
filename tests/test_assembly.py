"""
Golden tests for final oligo assembly.

`assembly.run()` builds the molecule that actually gets synthesized, so these
pin the exact splint and padlock sequences rather than just their shape. The
expected values in `tests/data/assembly/golden_oligos.json` were verified
byte-identical to `probegen/2_assemble_manifest.py` at `fishtools@906c15d`,
run single-threaded.

The fixture is six real probe pairs for each of two octopus transcripts, taken
from a production run. RepeatMasker is skipped so the tests need no external
tools.
"""

import json
import shutil
from pathlib import Path

import polars as pl
import pytest
from Bio import Seq
from Bio.Restriction import BamHI, KpnI

from mkprobes.assembly import run
from mkprobes.codebook.codebook import ProbeSet

# Aliased: pytest would otherwise collect this helper as a test case.
from mkprobes.starmap.starmap import test_splint_padlock as splint_padlock_pairs

FIXTURE = Path(__file__).parent / "data" / "assembly"
GOLDEN = json.loads((FIXTURE / "golden_oligos.json").read_text())


def assemble(tmp_path: Path) -> pl.DataFrame:
    """Runs the real assembly over a private copy of the fixture."""
    work = tmp_path / "panel"
    shutil.copytree(FIXTURE, work)
    probeset = ProbeSet(**json.loads((FIXTURE / "manifest.json").read_text())[0])
    run(work, probeset, n=6, skip_repeatmasker=True)
    return pl.read_parquet(work / "generated" / f"{probeset.name}.parquet")


@pytest.fixture(scope="module")
def assembled(tmp_path_factory: pytest.TempPathFactory) -> pl.DataFrame:
    return assemble(tmp_path_factory.mktemp("assembled"))


class TestGoldenOligos:
    def test_padlock_sequences_are_unchanged(self, assembled: pl.DataFrame):
        assert assembled["padlockcons"].to_list() == GOLDEN["padlockcons"]

    def test_splint_sequences_are_unchanged(self, assembled: pl.DataFrame):
        assert assembled["splintcons"].to_list() == GOLDEN["splintcons"]

    def test_probe_pair_count(self, assembled: pl.DataFrame):
        assert len(assembled) == GOLDEN["n_probe_pairs"]

    def test_deterministic_across_runs(self, tmp_path: Path, assembled: pl.DataFrame):
        # The head splint and the ATAAT filler are drawn from generators shared
        # across rows. Consuming them inside a polars UDF made the oligo pool
        # vary between runs, because map_elements is threaded.
        again = assemble(tmp_path)
        assert again["padlockcons"].to_list() == assembled["padlockcons"].to_list()
        assert again["splintcons"].to_list() == assembled["splintcons"].to_list()


class TestAssemblyInvariants:
    """Geometry and chemistry the assay depends on, asserted independently of the
    golden values so a regenerated golden file cannot quietly break them."""

    def test_padlock_length_within_synthesis_window(self, assembled: pl.DataFrame):
        assert assembled["padlockcons"].str.len_chars().is_between(139, 150).all()

    def test_splint_backfilled_to_target_length(self, assembled: pl.DataFrame):
        assert (assembled["splintcons"].str.len_chars() == 148).all()

    @pytest.mark.parametrize("column", ["spl_cut", "pad_cut"])
    def test_probe_core_is_free_of_restriction_sites(self, assembled: pl.DataFrame, column: str):
        # The header and footer deliberately carry BamHI/KpnI sites so the
        # construct can be excised. The core between them must not, or the
        # digest would cut the probe itself.
        for seq in assembled[column]:
            core = Seq.Seq(seq)
            assert not BamHI.search(core), f"BamHI site in {column}: {seq}"
            assert not KpnI.search(core), f"KpnI site in {column}: {seq}"

    def test_ligation_junction_pairs(self, assembled: pl.DataFrame):
        # The splint must template the padlock's two ends, 6 nt on each side.
        for splint, padlock in zip(assembled["spl_cut"], assembled["pad_cut"]):
            assert splint_padlock_pairs(splint, padlock, lengths=(6, 6))

    def test_ligation_junction_survives_double_digest(self, assembled: pl.DataFrame):
        # What matters at the bench is the fragment released by the KpnI/BamHI
        # double digest, not the ordered oligo.
        def excise(seq: str) -> str:
            return str(BamHI.catalyze(KpnI.catalyze(Seq.Seq(seq))[1])[0])

        for splint, padlock in zip(assembled["splintcons"], assembled["padlockcons"]):
            assert splint_padlock_pairs(excise(splint), excise(padlock), lengths=(6, 6))

    def test_sequences_are_unambiguous(self, assembled: pl.DataFrame):
        for column in ("splintcons", "padlockcons"):
            assert not assembled[column].str.to_uppercase().str.contains("N").any()
