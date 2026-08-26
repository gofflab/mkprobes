"""
Tests for the vendored seqspec files and the tooling that reads them.

The specs describe the molecule that gets synthesized, so the important test is
not that they parse but that the pool `assembly.run()` actually produces
satisfies them - and that a pool with something wrong with it does not. The
golden oligos in `tests/data/assembly/` are the same production probe pairs the
assembly tests pin, so the two cannot drift apart silently.
"""

import importlib.util
import json
from pathlib import Path

import pytest

from mkprobes.oligospec import (
    POOL_SPEC,
    READOUT_SPEC,
    SPEC_DIR,
    available_bcidx,
    detect_bcidx,
    exact_length,
    load_spec,
    match_regions,
    pool_spec_path,
    representative_widths,
    validate_pool,
    working_probe,
)

FIXTURE = Path(__file__).parent / "data" / "assembly"
GOLDEN = json.loads((FIXTURE / "golden_oligos.json").read_text())
CODEBOOK = json.loads((FIXTURE / "codebook.json").read_text())

#: seqspec reports this on any spec with no sequencing reads, because it
#: compares read file counts with `len(set([])) != 1`. SOLAR is an imaging
#: assay and has no reads at all, so it is expected rather than a defect.
NO_READS_ERROR = "check_read_file_count"


@pytest.fixture
def pool(tmp_path: Path) -> Path:
    """The golden probe pairs written out the way `assemble gen` writes a pool."""
    lines = [
        oligo
        for splint, padlock in zip(GOLDEN["splintcons"], GOLDEN["padlockcons"])
        for oligo in (splint, padlock)
    ]
    path = tmp_path / "golden_final.txt"
    path.write_text("\n".join(lines) + "\n")
    return path


def _generator():
    """Loads scripts/generate_seqspec.py, which is not importable as a module."""
    path = Path(__file__).parent.parent / "scripts" / "generate_seqspec.py"
    spec = importlib.util.spec_from_file_location("generate_seqspec", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# --------------------------------------------------------------------------- #
# The specs themselves
# --------------------------------------------------------------------------- #


def test_specs_are_checked_in_as_generated():
    """
    The vendored specs must match what the generator produces from the
    chemistry tables. If someone edits headerfooter.csv or the readout table
    without re-running the generator, this is what catches it - for every
    index, not just the one the tests happen to exercise.
    """
    rendered = _generator().render()
    on_disk = {path.name for path in SPEC_DIR.iterdir()}
    assert on_disk == set(rendered), (
        "The set of files in data/seqspec differs from what "
        "scripts/generate_seqspec.py produces. Re-run it and commit the result."
    )
    differing = [name for name, text in rendered.items() if (SPEC_DIR / name).read_text() != text]
    assert not differing, (
        f"{len(differing)} vendored file(s) differ from the generator, first: {differing[0]}. "
        "Re-run `python scripts/generate_seqspec.py` and commit the result."
    )


def test_one_pool_spec_per_bcidx():
    """`bcidx` is a manifest field, so every value it accepts needs a spec."""
    from mkprobes.init_project import max_bcidx

    assert available_bcidx() == list(range(max_bcidx() + 1))


#: Every region `bcidx` changes. Wider than the primer sites alone, because the
#: scars those sites leave are part of the working probe, and the splint's clamp
#: has to template the padlock's ends - which are those scars.
BCIDX_DEPENDENT = {
    "splint_primer_5p",
    "splint_spacer_5p",
    "splint_clamp_fix",
    "splint_clamp_3p",
    "splint_primer_3p",
    "padlock_primer_5p",
    "padlock_spacer_5p",
    "padlock_ligation_3p",
    "padlock_ligation_end",
    "padlock_primer_3p",
}

#: The recognition sites and the bases they leave on the probe. These are enzyme
#: sequence, not design sequence, so unlike everything around them they must be
#: identical at every index.
RESTRICTION_REGIONS = {
    "splint_kpni_site": "GGTAC",
    "splint_kpni_retained": "C",
    "splint_bamhi_site": "GATCC",
    "padlock_kpni_site": "GGTAC",
    "padlock_kpni_retained": "C",
    "padlock_bamhi_site": "GATCC",
}


def test_restriction_sites_are_labelled_and_constant():
    """
    Each site is its own region, spelled out, and the same at every index -
    it is the enzyme's sequence, not the panel's.
    """
    for bcidx in (0, 1, 7, max(available_bcidx())):
        leaves = {leaf.region_id: leaf.sequence for leaf in load_spec(pool_spec_path(bcidx)).root.leaves()}
        for region_id, sequence in RESTRICTION_REGIONS.items():
            assert leaves[region_id] == sequence, f"{region_id} at bcidx {bcidx}"


def test_each_site_straddles_the_cut_with_one_base_left_on_the_probe():
    """
    The point of the layout: five of each site's six bases leave with the handle
    that is thrown away, and exactly one stays on the working probe. If a change
    ever moved a cut, the probe would carry more scar than intended.
    """
    spec = load_spec(pool_spec_path(0))
    for kind, site, enzyme in (
        ("splint", "GGTACC", "kpni"),
        ("padlock", "GGTACC", "kpni"),
        ("splint", "GGATCC", "bamhi"),
        ("padlock", "GGATCC", "bamhi"),
    ):
        oligo = spec.root.find(f"{kind}_oligo")
        ordered = [leaf.region_id for leaf in oligo.leaves()]
        sequences = {leaf.region_id: leaf.sequence for leaf in oligo.leaves()}
        site_id = f"{kind}_{enzyme}_site"
        neighbour = ordered[ordered.index(site_id) + (1 if enzyme == "kpni" else -1)]

        # Reassembled across the cut, the six bases spell the whole site.
        pieces = (
            sequences[site_id] + sequences[neighbour][:1]
            if enzyme == "kpni"
            else sequences[neighbour][-1:] + sequences[site_id]
        )
        assert pieces == site, f"{kind} {enzyme}: {pieces}"
        assert len(sequences[site_id]) == 5, "five bases leave with the handle"


def test_the_retained_base_is_inside_the_working_probe():
    """The one base that survives the digest has to be part of the probe."""
    spec = load_spec(pool_spec_path(0))
    for kind in ("splint", "padlock"):
        oligo = spec.root.find(f"{kind}_oligo")
        probe_ids = {leaf.region_id for leaf in oligo.find(f"{kind}_probe").leaves()}
        assert f"{kind}_kpni_retained" in probe_ids
        # ...and the site itself is not, since it goes with the handle.
        assert f"{kind}_kpni_site" not in probe_ids
        assert f"{kind}_bamhi_site" not in probe_ids


def test_bcidx_changes_exactly_the_documented_regions():
    """
    Pins what an index actually controls. The specs say ten regions depend on
    `bcidx`; if that set ever grew, a per-index spec would be silently wrong
    about the rest, and if it shrank the extra files would be pure weight.
    """
    baseline = {leaf.region_id: leaf.sequence for leaf in load_spec(pool_spec_path(0)).root.leaves()}
    varying: set[str] = set()
    for bcidx in available_bcidx()[1:]:
        other = {leaf.region_id: leaf.sequence for leaf in load_spec(pool_spec_path(bcidx)).root.leaves()}
        assert set(other) == set(baseline), f"bcidx {bcidx} has a different set of regions"
        varying |= {rid for rid, seq in other.items() if seq != baseline[rid]}
    assert varying == BCIDX_DEPENDENT


def test_detects_bcidx(pool: Path):
    assert detect_bcidx(pool) == 0


def test_detects_no_bcidx_for_a_foreign_pool(tmp_path: Path):
    path = tmp_path / "not_a_pool.txt"
    path.write_text("ACGT\nTGCA\n")
    assert detect_bcidx(path) is None


def test_wrong_bcidx_is_rejected(pool: Path):
    """
    A pool checked against the wrong index must fail, not quietly pass. This is
    the whole reason there is a spec per index.
    """
    report = validate_pool(pool, load_spec(pool_spec_path(7)))
    assert len(report.problems) >= report.n_pairs


def test_onlist_is_the_vendored_readout_table():
    from mkprobes.codebook.finalconstruct import READOUTS

    onlist = (SPEC_DIR / "solar_readouts.txt").read_text().split()
    assert onlist == [READOUTS[i] for i in sorted(READOUTS)]
    assert len(onlist) == 49
    assert {len(seq) for seq in onlist} == {20}


@pytest.mark.parametrize("path", [POOL_SPEC, READOUT_SPEC, "SPEC_BCIDX_7"])
def test_conforms_to_seqspec(path):
    """Validates against the real seqspec implementation, when it is installed."""
    seqspec_check = pytest.importorskip("seqspec.seqspec_check", reason="seqspec not installed")
    from seqspec.utils import load_spec as seqspec_load

    path = pool_spec_path(7) if path == "SPEC_BCIDX_7" else path
    errors = seqspec_check.check(seqspec_load(path, strict=False))
    unexpected = [e for e in errors if e["error_type"] != NO_READS_ERROR]
    assert not unexpected, unexpected


@pytest.mark.parametrize("oligo_id", ["splint_oligo", "padlock_oligo"])
def test_oligos_are_148_nt(oligo_id: str):
    """
    Folding the coupled regions back in has to reproduce the length assembly
    asserts on, or the spec describes a molecule nobody ordered.
    """
    oligo = load_spec(POOL_SPEC).root.find(oligo_id)
    assert exact_length(oligo) == 148
    assert sum(representative_widths(oligo).values()) == 148
    assert oligo.min_len <= 148 <= oligo.max_len


# --------------------------------------------------------------------------- #
# Validating a real pool
# --------------------------------------------------------------------------- #


def test_golden_pool_validates(pool: Path):
    report = validate_pool(pool, load_spec(POOL_SPEC), CODEBOOK)
    assert report.problems == []
    assert report.n_pairs == GOLDEN["n_probe_pairs"]


def test_readouts_are_a_permutation_of_the_codebook_entry(pool: Path):
    """
    `construct_encoding` cycles through permutations of a gene's three bits, so
    position on the oligo carries no meaning and only the set is fixed. If that
    ever became a fixed order, the spec's region descriptions would be wrong.
    """
    report = validate_pool(pool, load_spec(POOL_SPEC), CODEBOOK)
    assigned = {tuple(sorted(bits)) for name, bits in CODEBOOK.items() if not name.startswith("Blank-")}
    assert {tuple(sorted(codes)) for codes in report.readouts} == assigned
    # Same gene, more than one arrangement: the permutation cycling is real.
    assert len({codes for codes in report.readouts}) > len(assigned)


def test_working_probe_excludes_the_amplification_handles():
    spec = load_spec(POOL_SPEC)
    padlock_region = spec.root.find("padlock_oligo")
    padlock = GOLDEN["padlockcons"][0].upper()
    parts = match_regions(padlock, padlock_region)
    assert parts is not None

    body = working_probe(parts, padlock_region)
    assert body in padlock
    assert len(body) < len(padlock)
    # The handles carry the restriction sites; the released probe must not.
    assert "GGATCC" not in body and "GGTACC" not in body


# --------------------------------------------------------------------------- #
# Rejecting a broken pool
# --------------------------------------------------------------------------- #


def _write(tmp_path: Path, splints, padlocks) -> Path:
    lines = [oligo for pair in zip(splints, padlocks) for oligo in pair]
    path = tmp_path / "pool.txt"
    path.write_text("\n".join(lines) + "\n")
    return path


def test_rejects_wrong_length(tmp_path: Path):
    splints = [GOLDEN["splintcons"][0][:-1], *GOLDEN["splintcons"][1:]]
    report = validate_pool(_write(tmp_path, splints, GOLDEN["padlockcons"]))
    assert any("exactly 148 nt" in problem for problem in report.problems)


def test_rejects_mutated_constant_region(tmp_path: Path):
    padlock = GOLDEN["padlockcons"][0]
    # One base of the 5' amplification primer, kept the same length.
    broken = padlock[:5] + ("T" if padlock[5].upper() != "T" else "A") + padlock[6:]
    report = validate_pool(_write(tmp_path, GOLDEN["splintcons"], [broken, *GOLDEN["padlockcons"][1:]]))
    assert any("does not match the padlock_oligo structure" in problem for problem in report.problems)


def test_rejects_readout_not_in_the_table(tmp_path: Path):
    padlock = GOLDEN["padlockcons"][0]
    broken = padlock[:60] + "G" * 20 + padlock[80:]
    report = validate_pool(_write(tmp_path, GOLDEN["splintcons"], [broken, *GOLDEN["padlockcons"][1:]]))
    assert report.problems


def test_rejects_restriction_site_in_the_probe_body(tmp_path: Path):
    padlock = GOLDEN["padlockcons"][0]
    # A BamHI site inside the homology arm, which the digest would cut.
    broken = padlock[:25] + "GGATCC" + padlock[31:]
    report = validate_pool(_write(tmp_path, GOLDEN["splintcons"], [broken, *GOLDEN["padlockcons"][1:]]))
    assert any("BamHI" in problem for problem in report.problems)


def test_rejects_corrupted_filler(tmp_path: Path):
    """
    The padlock's filler varies in length but not in content, so arbitrary
    bases there are a real defect - and one that keeps the oligo 148 nt, so
    only matching the filler's actual sequence catches it.
    """
    padlock_region = load_spec(POOL_SPEC).root.find("padlock_oligo")
    # A padlock with a full-length arm has no filler at all; find one that does.
    index, fill = next(
        (i, parts["padlock_fill"])
        for i, padlock in enumerate(GOLDEN["padlockcons"])
        if (parts := match_regions(padlock.upper(), padlock_region)) and parts["padlock_fill"]
    )
    padlock = GOLDEN["padlockcons"][index].upper()

    start = padlock.rindex(fill)
    broken = padlock[:start] + "G" * len(fill) + padlock[start + len(fill) :]
    assert len(broken) == len(padlock)

    padlocks = list(GOLDEN["padlockcons"])
    padlocks[index] = broken
    report = validate_pool(_write(tmp_path, GOLDEN["splintcons"], padlocks))
    assert report.problems


def test_rejects_odd_oligo_count(tmp_path: Path):
    path = tmp_path / "pool.txt"
    path.write_text("\n".join([GOLDEN["splintcons"][0], GOLDEN["padlockcons"][0], GOLDEN["splintcons"][1]]))
    assert any("odd" in problem for problem in validate_pool(path).problems)


def test_rejects_a_blank_codeword(pool: Path):
    """Blanks measure the false-positive rate, so synthesising one is a serious error."""
    codebook = {"Real": CODEBOOK["Och.687.1"], "Blank-1": CODEBOOK["Och.958.1"]}
    report = validate_pool(pool, load_spec(POOL_SPEC), codebook)
    assert any("Blank" in problem for problem in report.problems)


def test_rejects_codeword_absent_from_the_codebook(pool: Path):
    report = validate_pool(pool, load_spec(POOL_SPEC), {"Only": [10, 11, 12]})
    assert any("not assigned to any target" in problem for problem in report.problems)


# --------------------------------------------------------------------------- #
# Orientation
# --------------------------------------------------------------------------- #


def test_recognition_sites_are_palindromes():
    """
    Why the sites are drawn without a direction: GGATCC and GGTACC each read the
    same on both strands, so there is no orientation to show, and the enzyme
    cuts the duplex rather than one strand.
    """
    complement = str.maketrans("ACGT", "TGCA")
    for site in ("GGATCC", "GGTACC"):
        assert site.translate(complement)[::-1] == site


def test_sites_read_5_to_3_as_written():
    """
    The site has to spell itself out in the direction the oligo is written, or
    the enzyme would not cut where the spec claims. Reassembled across the cut,
    each one appears exactly once per oligo.
    """
    spec = load_spec(pool_spec_path(0))
    for kind in ("splint", "padlock"):
        oligo = spec.root.find(f"{kind}_oligo")
        written = "".join(leaf.sequence for leaf in oligo.leaves())
        assert written.count("GGTACC") == 1, f"{kind}: KpnI site not readable 5'->3'"
        assert written.count("GGATCC") == 1, f"{kind}: BamHI site not readable 5'->3'"


def test_strand_assignment():
    """Sites undirected, reverse-primer sites on the opposite strand, rest forward."""
    from mkprobes.oligospec import _strand

    oligo = load_spec(pool_spec_path(0)).root.find("splint_oligo")
    strands = {leaf.region_id: _strand(leaf) for leaf in oligo.leaves()}
    assert strands["splint_kpni_site"] == 0
    assert strands["splint_bamhi_site"] == 0
    assert strands["splint_primer_3p"] == -1
    assert strands["splint_primer_5p"] == +1
    assert strands["splint_arm"] == +1
    # The retained base is a probe base, not part of the site: it keeps a direction.
    assert strands["splint_kpni_retained"] == +1


# --------------------------------------------------------------------------- #
# Drawing
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("path", [POOL_SPEC, READOUT_SPEC])
def test_draws(path: Path, tmp_path: Path):
    pytest.importorskip("dna_features_viewer", reason="viz extra not installed")
    from mkprobes.oligospec import draw_spec

    out = draw_spec(load_spec(path), tmp_path / "spec.png")
    assert out.exists() and out.stat().st_size > 0
