"""Tests for annotation-driven transcript selection on custom datasets."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import polars as pl
import pytest

from mkprobes.ext.dataset import Dataset
from mkprobes.ext.external_data import ExternalData, MockGTF
from mkprobes.genes.chkgenes import get_transcripts, get_transcripts_generic


class FakeFasta(dict):
    """dict-like standing in for pyfastx.Fasta: keys() + len(fa[id])."""


@pytest.fixture
def generic_ds(tmp_path: Path) -> Dataset:
    gtf = pl.DataFrame({
        "feature": ["transcript"] * 4,
        "gene_name": ["Och.1", "Och.1", "Och.2", "Och.3"],
        "gene_id": ["Och.1", "Och.1", "Och.2", "Och.3"],
        "transcript_id": ["Och.1.1", "Och.1.2", "Och.2.1", "Och.3.1"],
        "transcript_name": ["Och.1.1", "Och.1.2", "Och.2.1", "Och.3.1"],
    })
    fa = FakeFasta({
        "Och.1.1": "A" * 2000,
        "Och.1.2": "A" * 3500,
        "Och.2.1": "A" * 1800,
        "Och.3.1": "A" * 2400,
    })
    mock_ed = MagicMock(spec=ExternalData)
    mock_ed.gtf = gtf
    mock_ed.fa = fa
    return Dataset(path=tmp_path, external_data=mock_ed, species="octopus")


class TestGetTranscriptsGeneric:
    def test_transcript_id_passthrough(self, generic_ds: Dataset):
        res = get_transcripts_generic(generic_ds, ["Och.1.1"], mode="longest")
        assert res["transcript_id"].to_list() == ["Och.1.1"]

    def test_longest_picks_longest_sequence(self, generic_ds: Dataset):
        res = get_transcripts_generic(generic_ds, ["Och.1"], mode="longest")
        assert res["transcript_id"].to_list() == ["Och.1.2"]

    def test_all_returns_every_isoform(self, generic_ds: Dataset):
        res = get_transcripts_generic(generic_ds, ["Och.1"], mode="all")
        assert sorted(res["transcript_id"].to_list()) == ["Och.1.1", "Och.1.2"]

    def test_multiple_tokens(self, generic_ds: Dataset):
        res = get_transcripts_generic(generic_ds, ["Och.1", "Och.3.1"], mode="longest")
        assert res["transcript_id"].to_list() == ["Och.1.2", "Och.3.1"]

    def test_unresolved_raises_with_suggestions(self, generic_ds: Dataset, caplog):
        with pytest.raises(ValueError, match="could not be resolved"):
            get_transcripts_generic(generic_ds, ["Och.99"], mode="longest")

    def test_annotation_table_resolution(self, generic_ds: Dataset, tmp_path: Path):
        ortho = tmp_path / "orthologs.tsv"
        ortho.write_text(
            "transcript_id\thuman_symbol\nOch.2.1\tSHANK3\nOch.3.1\tGRIN1\n"
        )
        generic_ds.annotation_paths = {"orthologs": ortho}
        res = get_transcripts_generic(generic_ds, ["SHANK3"], mode="longest")
        assert res["transcript_id"].to_list() == ["Och.2.1"]
        # Case-insensitive
        res = get_transcripts_generic(generic_ds, ["shank3"], mode="longest")
        assert res["transcript_id"].to_list() == ["Och.2.1"]

    def test_mock_gtf_passthrough_and_error(self, tmp_path: Path):
        mock_ed = MagicMock(spec=ExternalData)
        mock_ed.gtf = MockGTF()
        mock_ed.fa = FakeFasta({"seq1": "A" * 100})
        ds = Dataset(path=tmp_path, external_data=mock_ed, species="x")
        res = get_transcripts_generic(ds, ["seq1"], mode="longest")
        assert res["transcript_id"].to_list() == ["seq1"]
        with pytest.raises(ValueError, match="no GTF"):
            get_transcripts_generic(ds, ["nonexistent"], mode="longest")


class TestGetTranscriptsDispatch:
    def test_canonical_falls_back_to_longest_on_generic(self, generic_ds: Dataset):
        res = get_transcripts(generic_ds, ["Och.1"], mode="canonical")
        assert res["transcript_id"].to_list() == ["Och.1.2"]

    def test_reference_only_mode_rejected_on_generic(self, generic_ds: Dataset):
        with pytest.raises(ValueError, match="requires a human/mouse reference"):
            get_transcripts(generic_ds, ["Och.1"], mode="appris")


class TestCandidatesCliWiring:
    def test_candidates_cli_uses_load_dataset(self, tmp_path: Path):
        from click.testing import CliRunner

        from mkprobes import candidates as cand_mod

        sentinel = MagicMock(spec=Dataset)
        with (
            patch.object(cand_mod, "load_dataset", return_value=sentinel) as mock_load,
            patch.object(cand_mod, "get_candidates") as mock_get,
        ):
            res = CliRunner().invoke(
                cand_mod.candidates,
                [str(tmp_path), "--gene", "Och.1.1", "--output", str(tmp_path / "out")],
            )
        assert res.exit_code == 0, res.output
        mock_load.assert_called_once()
        assert mock_get.call_args[0][0] is sentinel
