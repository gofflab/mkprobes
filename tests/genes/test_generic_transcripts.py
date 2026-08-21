"""Tests for annotation-driven transcript selection on custom datasets."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import polars as pl
import pytest

from mkprobes.ext.dataset import Dataset
from mkprobes.ext.external_data import ExternalData, MockGTF
from mkprobes.genes.chkgenes import (
    _resolve_via_annotations,
    get_transcripts,
    get_transcripts_generic,
)


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
            # `candidates` pre-flights bowtie2 before doing any work. This test
            # is about wiring and mocks the work away, so it must not also
            # require the aligner to be installed. The pre-flight itself is
            # covered in tests/test_cli_errors.py.
            patch("mkprobes.ext.ingest.check_external_tools"),
        ):
            res = CliRunner().invoke(
                cand_mod.candidates,
                [str(tmp_path), "--gene", "Och.1.1", "--output", str(tmp_path / "out")],
            )
        assert res.exit_code == 0, res.output
        mock_load.assert_called_once()
        assert mock_get.call_args[0][0] is sentinel


class TestGeneNameColumn:
    """
    Target lists are written in whatever names the annotation carries. Naming
    the column that holds them keeps lookup off the rest of a wide table -
    which on a real ortholog table meant scanning protein sequences and ESM
    embeddings for something that looked like a gene name.
    """

    def _table(self) -> pl.DataFrame:
        return pl.DataFrame({
            "transcript_id": ["t1", "t2", "t3"],
            "gene_id": ["g1", "g2", "g3"],
            # Ortholog tables routinely map one transcript to several symbols.
            "human_name": ["SOX2", "UBE2A,UBE2B", "RNF185,,RNF5,,"],
            "local_name": ["och-sox2", "och-ube2", "och-rnf5"],
            "notes": ["SOX2", "", ""],
        })

    def _dataset(self, column: str | None):
        dataset = MagicMock(spec=Dataset)
        dataset.annotation_paths = {"annot": Path("x")}
        dataset.annotation = lambda name: self._table()
        dataset.gene_name_column = column
        return dataset

    def test_named_column_resolves(self):
        assert _resolve_via_annotations(self._dataset("human_name"), "SOX2") == ["t1", "g1"]

    def test_multi_valued_cells_match_per_entry(self):
        # An exact match against the whole cell would never fire on these.
        assert _resolve_via_annotations(self._dataset("human_name"), "UBE2B") == ["t2", "g2"]
        assert _resolve_via_annotations(self._dataset("human_name"), "RNF5") == ["t3", "g3"]

    def test_matching_is_case_insensitive(self):
        assert _resolve_via_annotations(self._dataset("human_name"), "sox2") == ["t1", "g1"]

    def test_other_columns_are_not_consulted(self):
        # `local_name` carries och-sox2, but the dataset says human_name.
        assert _resolve_via_annotations(self._dataset("human_name"), "och-sox2") == []

    def test_a_different_column_selects_different_names(self):
        assert _resolve_via_annotations(self._dataset("local_name"), "och-sox2") == ["t1", "g1"]
        assert _resolve_via_annotations(self._dataset("local_name"), "SOX2") == []

    def test_without_a_named_column_every_column_is_searched(self):
        # Historical behaviour, kept for datasets built before the option.
        every = self._dataset(None)
        assert _resolve_via_annotations(every, "SOX2") == ["t1", "g1"]
        assert _resolve_via_annotations(every, "och-sox2") == ["t1", "g1"]

    def test_unknown_column_matches_nothing_rather_than_raising(self):
        assert _resolve_via_annotations(self._dataset("no_such_column"), "SOX2") == []
