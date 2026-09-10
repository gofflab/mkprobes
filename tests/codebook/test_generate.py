"""Tests for package-integrated codebook generation (mkprobes make-codebook)."""

import json
from itertools import chain
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
from click.testing import CliRunner

from mkprobes.codebook.generate import (
    FORBIDDEN,
    ORDER,
    choose_bits,
    discover_matrices,
    make_codebook,
    make_codebook_cli,
    offset_for_codebook,
    resolve_expression,
    resolve_offset,
)

GENES = ["Och.687.1", "Och.958.1", "Och.576.10"]


class TestMakeCodebook:
    def test_auto_sizes_to_10_bits_for_small_panel(self):
        cb = make_codebook(GENES, seed=0)
        # 10-bit on3 dist2 code has capacity 120 (legacy parity).
        assert len(cb) == 120
        n_blanks = sum(k.startswith("Blank") for k in cb)
        assert n_blanks == 117
        assert set(GENES) <= set(cb)

    def test_deterministic_with_seed(self):
        assert make_codebook(GENES, seed=0) == make_codebook(GENES, seed=0)
        assert make_codebook(GENES, seed=0) != make_codebook(GENES, seed=1)

    def test_three_distinct_bits_per_target(self):
        cb = make_codebook(GENES, seed=0)
        for gene, bits in cb.items():
            assert len(bits) == 3 and len(set(bits)) == 3, gene
            assert bits == sorted(bits)

    def test_no_forbidden_codewords_on_genes(self):
        cb = make_codebook(GENES, seed=0)
        for gene, bits in cb.items():
            if not gene.startswith("Blank"):
                assert tuple(bits) not in FORBIDDEN, gene

    def test_explicit_n_bits(self):
        cb = make_codebook(GENES, n_bits=12, seed=0)
        used = set(chain.from_iterable(cb.values()))
        assert len(cb) > 120  # 12-bit capacity exceeds 10-bit's

    def test_existing_codebook_derives_offset_and_rejects_overlap(self):
        first = make_codebook(GENES, n_bits=10, seed=0)
        second = make_codebook(["Och.9.1", "Och.10.1"], n_bits=10, existing_codebook=first, seed=0)
        first_bits = set(chain.from_iterable(first.values()))
        second_bits = set(chain.from_iterable(second.values()))
        assert not first_bits & second_bits

        with pytest.raises(ValueError, match="overlap"):
            make_codebook(GENES, n_bits=10, existing_codebook=first, seed=0)

    def test_offset_is_a_position_in_the_readout_order(self):
        # A pooled second panel starts where the first stopped: offset 10
        # after a 10-bit panel. The offset counts bit positions, so the
        # readout IDs it lands on are the interleaved ones, not 11..20.
        cb = make_codebook(GENES, n_bits=10, offset=10, seed=0)
        used = set(chain.from_iterable(cb.values()))
        assert used == set(ORDER[10:20]) == {12, 20, 5, 13, 21, 6, 14, 22, 7, 15}
        assert not used & set(chain.from_iterable(make_codebook(GENES, n_bits=10, seed=0).values()))

    def test_offset_and_existing_mutually_exclusive(self):
        first = make_codebook(GENES, n_bits=10, seed=0)
        with pytest.raises(ValueError, match="either offset or existing"):
            make_codebook(["X1"], offset=5, existing_codebook=first)

    def test_expression_informed_beats_or_matches_seed0(self):
        # Skewed expression: the optimizer must find an assignment whose
        # per-bit load entropy is at least that of the plain seed-0 shuffle.
        rng = np.random.default_rng(7)
        genes = [f"G{i}" for i in range(40)]
        expression = {g: float(v) for g, v in zip(genes, rng.lognormal(3, 2, len(genes)))}

        from mkprobes.codebook.codebook import CodebookPicker, _entropy

        matrix = discover_matrices()[10]

        def load_entropy(codebook: dict) -> float:
            # entropy of per-bit expression load over gene (non-blank) entries
            loads: dict[int, float] = {}
            for g, bits in codebook.items():
                if g.startswith("Blank"):
                    continue
                for b in bits:
                    loads[b] = loads.get(b, 0.0) + expression[g]
            return _entropy(list(loads.values()))

        informed = make_codebook(genes, expression=expression, iterations=50)
        plain = make_codebook(genes, seed=0)
        assert load_entropy(informed) >= load_entropy(plain) - 1e-9

    def test_expression_missing_target_raises(self):
        with pytest.raises(ValueError, match="missing"):
            make_codebook(GENES, expression={"Och.687.1": 1.0})


class TestResolveExpression:
    def _table(self, tmp_path: Path, text: str, name: str = "expr.tsv") -> Path:
        p = tmp_path / name
        p.write_text(text)
        return p

    def test_from_file_single_numeric_column(self, tmp_path: Path):
        p = self._table(tmp_path, "transcript_id\tfpkm\nOch.687.1\t10.0\nOch.958.1\t5.0\nOch.576.10\t2.5\n")
        vals = resolve_expression(None, str(p), GENES)
        assert vals == {"Och.687.1": 10.0, "Och.958.1": 5.0, "Och.576.10": 2.5}

    def test_missing_targets_filled_with_median(self, tmp_path: Path):
        p = self._table(tmp_path, "transcript_id\tfpkm\nOch.687.1\t10.0\nOch.958.1\t2.0\n")
        vals = resolve_expression(None, str(p), GENES)
        assert vals["Och.576.10"] == pytest.approx(6.0)  # median of 10, 2

    def test_ambiguous_numeric_columns_require_choice(self, tmp_path: Path):
        p = self._table(tmp_path, "transcript_id\tfpkm\ttpm\nOch.687.1\t10.0\t8.0\n")
        with pytest.raises(ValueError, match="expression-column"):
            resolve_expression(None, str(p), ["Och.687.1"])
        vals = resolve_expression(None, str(p), ["Och.687.1"], column="tpm")
        assert vals["Och.687.1"] == 8.0

    def test_gene_id_column_also_matches(self, tmp_path: Path):
        p = self._table(tmp_path, "gene_id\tfpkm\nOch.687\t3.0\n")
        vals = resolve_expression(None, str(p), ["Och.687"])
        assert vals["Och.687"] == 3.0

    def test_registered_annotation_table_wins(self, tmp_path: Path):
        table = self._table(tmp_path, "transcript_id\tfpkm\nOch.687.1\t4.0\n")
        ds = MagicMock()
        ds.annotation_paths = {"fpkm": table}
        from mkprobes.ext.dataset import _read_annotation_table

        ds.annotation.side_effect = lambda name: _read_annotation_table(ds.annotation_paths[name])
        vals = resolve_expression(ds, "fpkm", ["Och.687.1"])
        assert vals["Och.687.1"] == 4.0

    def test_unknown_source_errors_with_available(self):
        ds = MagicMock()
        ds.annotation_paths = {"orthologs": Path("x")}
        with pytest.raises(ValueError, match="orthologs"):
            resolve_expression(ds, "nope", GENES)


def _manifest(path: Path, codebook: str = "codebook.json", offset: int = 10, name: str = "p") -> dict:
    return {"name": name, "species": "octopus", "codebook": codebook, "bcidx": 0, "offset": offset}


class TestOffsetFromManifest:
    """
    The manifest is where a panel says which bits it occupies. These pin that
    `make-codebook` takes the offset from there, matched by codebook path the
    way `run-panel` matches design settings.
    """

    def test_manifest_beside_output_supplies_offset(self, tmp_path: Path):
        (tmp_path / "manifest.json").write_text(json.dumps([_manifest(tmp_path)]))
        offset, source = offset_for_codebook(tmp_path / "codebook.json", None)
        assert (offset, source) == (10, tmp_path / "manifest.json")

    def test_no_manifest_means_no_offset(self, tmp_path: Path):
        assert offset_for_codebook(tmp_path / "codebook.json", None) == (None, None)

    def test_implicit_manifest_not_naming_the_codebook_is_ignored(self, tmp_path: Path):
        (tmp_path / "manifest.json").write_text(json.dumps([_manifest(tmp_path, codebook="other.json")]))
        assert offset_for_codebook(tmp_path / "codebook.json", None) == (None, None)

    def test_explicit_manifest_not_naming_the_codebook_is_an_error(self, tmp_path: Path):
        manifest = tmp_path / "m.json"
        manifest.write_text(json.dumps([_manifest(tmp_path, codebook="other.json")]))
        with pytest.raises(ValueError, match="no probe set whose codebook"):
            offset_for_codebook(tmp_path / "codebook.json", manifest)

    def test_disagreeing_duplicates_are_ambiguous(self, tmp_path: Path):
        entries = [_manifest(tmp_path, offset=0, name="a"), _manifest(tmp_path, offset=10, name="b")]
        (tmp_path / "manifest.json").write_text(json.dumps(entries))
        with pytest.raises(ValueError, match="different offsets"):
            offset_for_codebook(tmp_path / "codebook.json", None)

    def test_flag_overrides_manifest(self, tmp_path: Path):
        (tmp_path / "manifest.json").write_text(json.dumps([_manifest(tmp_path)]))
        assert resolve_offset(tmp_path / "codebook.json", None, 20, None) == 20
        assert resolve_offset(tmp_path / "codebook.json", None, None, None) == 10

    def test_existing_codebook_must_agree_with_manifest(self, tmp_path: Path):
        first = tmp_path / "first.json"
        first.write_text(json.dumps(make_codebook(GENES, n_bits=10, seed=0)))
        (tmp_path / "manifest.json").write_text(json.dumps([_manifest(tmp_path, offset=10)]))
        assert resolve_offset(tmp_path / "codebook.json", None, None, first) == 0  # derived later

        (tmp_path / "manifest.json").write_text(json.dumps([_manifest(tmp_path, offset=7)]))
        with pytest.raises(ValueError, match="implies 10"):
            resolve_offset(tmp_path / "codebook.json", None, None, first)

    def test_invalid_manifest_is_reported(self, tmp_path: Path):
        (tmp_path / "manifest.json").write_text('[{"name": "p"}]')
        with pytest.raises(ValueError, match="not a valid manifest"):
            offset_for_codebook(tmp_path / "codebook.json", None)


class TestMakeCodebookCli:
    def test_offset_comes_from_the_manifest(self, tmp_path: Path):
        genes_file = tmp_path / "genes.txt"
        genes_file.write_text("\n".join(GENES) + "\n")
        out = tmp_path / "codebook.json"
        (tmp_path / "manifest.json").write_text(json.dumps([_manifest(tmp_path, offset=10)]))

        res = CliRunner().invoke(make_codebook_cli, [str(tmp_path), str(genes_file), "-o", str(out)])

        assert res.exit_code == 0, res.output
        used = set(chain.from_iterable(json.loads(out.read_text()).values()))
        assert used == set(ORDER[10:20])

    def test_offset_flag_wins_over_the_manifest(self, tmp_path: Path):
        genes_file = tmp_path / "genes.txt"
        genes_file.write_text("\n".join(GENES) + "\n")
        out = tmp_path / "codebook.json"
        (tmp_path / "manifest.json").write_text(json.dumps([_manifest(tmp_path, offset=10)]))

        res = CliRunner().invoke(
            make_codebook_cli, [str(tmp_path), str(genes_file), "-o", str(out), "--offset", "20"]
        )

        assert res.exit_code == 0, res.output
        used = set(chain.from_iterable(json.loads(out.read_text()).values()))
        assert used == set(ORDER[20:30])

    def test_explicit_manifest_without_this_codebook_fails(self, tmp_path: Path):
        genes_file = tmp_path / "genes.txt"
        genes_file.write_text("\n".join(GENES) + "\n")
        manifest = tmp_path / "m.json"
        manifest.write_text(json.dumps([_manifest(tmp_path, codebook="other.json")]))

        res = CliRunner().invoke(
            make_codebook_cli, [str(tmp_path), str(genes_file), "--manifest", str(manifest)]
        )

        assert res.exit_code != 0
        assert "no probe set" in res.output

    def test_uninformed_run_writes_json(self, tmp_path: Path):
        genes_file = tmp_path / "genes.txt"
        genes_file.write_text("\n".join(GENES) + "\n")
        res = CliRunner().invoke(make_codebook_cli, [str(tmp_path), str(genes_file)])
        assert res.exit_code == 0, res.output
        cb = json.loads((tmp_path / "genes.codebook.json").read_text())
        assert set(GENES) <= set(cb)
        assert len(cb) == 120

    def test_expression_via_file(self, tmp_path: Path):
        genes_file = tmp_path / "genes.txt"
        genes_file.write_text("\n".join(GENES) + "\n")
        expr = tmp_path / "expr.tsv"
        expr.write_text("transcript_id\tfpkm\nOch.687.1\t10\nOch.958.1\t5\nOch.576.10\t1\n")
        with patch("mkprobes.codebook.generate.load_dataset") as mock_load:
            mock_load.return_value.annotation_paths = {}
            res = CliRunner().invoke(
                make_codebook_cli,
                [str(tmp_path), str(genes_file), "-e", str(expr), "--iterations", "10"],
            )
        assert res.exit_code == 0, res.output
        assert (tmp_path / "genes.codebook.json").exists()

    def test_duplicate_targets_rejected(self, tmp_path: Path):
        genes_file = tmp_path / "genes.txt"
        genes_file.write_text("A\nA\n")
        res = CliRunner().invoke(make_codebook_cli, [str(tmp_path), str(genes_file)])
        assert res.exit_code != 0
        # Naming the offending target is the point; a bare "duplicate" is not
        # actionable on a list of several hundred.
        assert "A" in res.output

    def test_comments_and_blank_lines_are_ignored(self, tmp_path: Path):
        genes_file = tmp_path / "genes.txt"
        genes_file.write_text("# my panel\nA\n\nB  # keep\n\n")
        out = tmp_path / "cb.json"
        res = CliRunner().invoke(make_codebook_cli, [str(tmp_path), str(genes_file), "-o", str(out)])
        assert res.exit_code == 0, res.output
        assert set(json.loads(out.read_text())) >= {"A", "B"}
