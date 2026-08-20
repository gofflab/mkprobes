"""Round-trip tests for defensive probe-name parsing (underscore-rich de novo IDs)."""

import polars as pl
import pytest

from mkprobes.utils.sequtils import probe_coord_exprs, probe_identity_exprs


def parse_coords(name: str) -> dict:
    return pl.DataFrame({"name": [name]}).with_columns(probe_coord_exprs()).to_dicts()[0]


def parse_identity(name: str) -> dict:
    return pl.DataFrame({"name": [name]}).with_columns(probe_identity_exprs()).to_dicts()[0]


class TestProbeCoordExprs:
    def test_reference_style(self):
        r = parse_coords("Sox2_ENSMUST00000123:100-150_splint")
        assert r["gene"] == "Sox2"
        assert r["transcript_ori"] == "ENSMUST00000123"
        assert (r["pos_start"], r["pos_end"]) == (100, 150)

    def test_reference_gene_with_dash(self):
        r = parse_coords("Nkx2-1_ENSMUST00000001:1-50_padlock")
        assert r["gene"] == "Nkx2-1"
        assert r["transcript_ori"] == "ENSMUST00000001"

    def test_stringtie_dotted_ids(self):
        r = parse_coords("Och.1.1_Och.1.1:10-55_splint")
        assert r["gene"] == "Och.1.1"
        assert r["transcript_ori"] == "Och.1.1"
        assert (r["pos_start"], r["pos_end"]) == (10, 55)

    def test_trinity_underscore_rich_ids(self):
        # The even-duplicate rule: generic names are always {X}_{X}, so any
        # underscore-rich X parses unambiguously.
        tid = "TRINITY_DN123_c0_g1_i1"
        r = parse_coords(f"{tid}_{tid}:5-60_padlock")
        assert r["gene"] == tid
        assert r["transcript_ori"] == tid
        assert (r["pos_start"], r["pos_end"]) == (5, 60)

    def test_maker_ids(self):
        tid = "maker-scaffold1-snap-gene-0.12-mRNA-1"
        r = parse_coords(f"{tid}_{tid}:1-44_splint")
        assert r["gene"] == tid
        assert r["transcript_ori"] == tid

    def test_ncbi_underscored_ids(self):
        tid = "rna-XM_012345.1"
        r = parse_coords(f"{tid}_{tid}:7-51_splint")
        assert r["gene"] == tid
        assert r["transcript_ori"] == tid

    def test_without_probe_type_suffix(self):
        r = parse_coords("Och.1.1_Och.1.1:10-55")
        assert r["gene"] == "Och.1.1"
        assert (r["pos_start"], r["pos_end"]) == (10, 55)


class TestProbeIdentityExprs:
    @pytest.mark.parametrize(
        "prefix",
        [
            "Sox2_ENSMUST00000123",
            "Och.1.1_Och.1.1",
            "TRINITY_DN123_c0_g1_i1_TRINITY_DN123_c0_g1_i1",
            "rna-XM_012345.1_rna-XM_012345.1",
        ],
    )
    @pytest.mark.parametrize("ptype", ["splint", "padlock"])
    def test_roundtrip(self, prefix: str, ptype: str):
        r = parse_identity(f"{prefix}:100-150_{ptype}")
        assert r["full_name"] == f"{prefix}:100-150"
        assert r["probe_type"] == ptype
