"""
Read the vendored seqspec files, validate assembled oligos against them, and draw them.

The specs in `data/seqspec/` describe the molecules `mkprobes assemble` emits.
This module is the code that uses them: `mkprobes validate-pool` checks an
ordered pool against the spec before it goes to a vendor, and
`mkprobes draw-spec` renders the spec as a labelled diagram.

Reading is done with plain YAML, so validating a pool does not require the
`seqspec` package to be installed. `seqspec` itself is only needed to re-check
that our files still conform to the specification, which is a developer task
(`tests/test_oligospec.py`).

CLI: ``mkprobes validate-pool <pool.txt>`` and ``mkprobes draw-spec``.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from importlib.resources import files
from pathlib import Path

import rich_click as click
import yaml
from loguru import logger

from .constants import SOLAR_RESTRICTION
from .starmap.starmap import test_splint_padlock

#: Vendored specs, written by `scripts/generate_seqspec.py`.
SPEC_DIR = Path(str(files("mkprobes") / "data" / "seqspec"))
#: The readout probes do not depend on `bcidx`, so there is one of these.
READOUT_SPEC = SPEC_DIR / "solar_readout_probes.seqspec.yaml"


def pool_spec_path(bcidx: int = 0) -> Path:
    """
    The pool spec for one `bcidx`.

    There is one per index because `bcidx` reaches further than it looks: it
    sets both primer binding sites on each oligo, the restriction scars those
    sites leave on the working probe, and the splint's clamp, which templates
    the padlock's ends - and those ends are the scars. Ten regions in all. A
    pool validated against the wrong index therefore fails on every oligo, for
    reasons that have nothing to do with the probes.
    """
    path = SPEC_DIR / f"solar_bcidx{bcidx}.seqspec.yaml"
    if not path.exists():
        available = available_bcidx()
        raise FileNotFoundError(
            f"No spec for bcidx {bcidx}. The vendored specs cover "
            f"{min(available)}-{max(available)}."
        )
    return path


def available_bcidx() -> list[int]:
    """Every `bcidx` a vendored spec exists for, in order."""
    return sorted(
        int(path.name.removeprefix("solar_bcidx").removesuffix(".seqspec.yaml"))
        for path in SPEC_DIR.glob("solar_bcidx*.seqspec.yaml")
    )


#: The default pool spec. `mkprobes init` and the test panel both use bcidx 0.
POOL_SPEC = SPEC_DIR / "solar_bcidx0.seqspec.yaml"

#: What each variable region is allowed to contain.
#:
#: seqspec spells every variable region `X`, so the alphabet lives here instead.
#: The defaults are deliberately tight - they are what `assembly.py` can
#: actually emit - because a pool that drifts outside them is the thing this
#: check exists to catch. Anything unlisted falls back to plain DNA.
ALPHABETS: dict[str, str] = {
    # itertools.cycle("ATAAT") in assembly.splint_pad
    "splint_pad": "AT",
    # rc of the padlock's head splint, which is drawn from A/C/T
    "splint_clamp_var": "AGT",
    # generate_head_splint draws from A/T/C, weighted to C
    "padlock_head_splint": "ACT",
    # the "AA"/"TA"/"AT"/"TT" separators in finalconstruct.stitch
    "padlock_spacer_1": "AT",
    "padlock_spacer_2": "AT",
}

#: Variable regions whose lengths are coupled, and to what total.
#:
#: `padpad` and `splint_pad` fill exactly as much as the arm leaves free, so the
#: two always sum to a constant. seqspec sums child ranges independently and so
#: cannot say this; the pool is nevertheless wrong if it does not hold.
COUPLED_LENGTHS: dict[str, tuple[tuple[str, str], int]] = {
    "splint_oligo": (("splint_pad", "splint_arm"), 33),
    "padlock_oligo": (("padlock_arm", "padlock_fill"), 27),
}

#: Variable regions whose content is a prefix of a fixed string.
#:
#: `assembly.padpad` right-truncates a constant to whatever length the arm
#: leaves free, so the filler varies in length but never in content. seqspec
#: has one notion of "variable" and spells it `X`, which would let any bases at
#: all through here; matching the prefixes instead keeps the check tight.
PREFIX_CONSTANTS: dict[str, str] = {"padlock_fill": "AATCACATAAAT"}


@dataclass
class Region:
    """One region of a spec, with its children."""

    region_id: str
    region_type: str
    name: str
    sequence_type: str
    sequence: str
    min_len: int
    max_len: int
    onlist: list[str] = field(default_factory=list)
    regions: list[Region] = field(default_factory=list)

    @property
    def is_leaf(self) -> bool:
        return not self.regions

    def leaves(self) -> list[Region]:
        if self.is_leaf:
            return [self]
        return [leaf for child in self.regions for leaf in child.leaves()]

    def find(self, region_id: str) -> Region:
        if self.region_id == region_id:
            return self
        for child in self.regions:
            try:
                return child.find(region_id)
            except KeyError:
                continue
        raise KeyError(f"No region {region_id!r} under {self.region_id!r}.")


@dataclass
class Spec:
    """A parsed seqspec file: its metadata and its top-level region."""

    name: str
    assay_id: str
    description: str
    root: Region
    path: Path


def _load_onlist(region: dict, spec_dir: Path) -> list[str]:
    onlist = region.get("onlist")
    if not onlist:
        return []
    target = spec_dir / onlist["url"]
    if not target.exists():
        raise FileNotFoundError(
            f"Region {region['region_id']!r} lists its allowed sequences in {onlist['url']}, "
            f"which is missing from {spec_dir}. Re-run scripts/generate_seqspec.py."
        )
    return [line.strip().upper() for line in target.read_text().splitlines() if line.strip()]


def _parse_region(raw: dict, spec_dir: Path) -> Region:
    return Region(
        region_id=raw["region_id"],
        region_type=raw["region_type"],
        name=raw.get("name", raw["region_id"]),
        sequence_type=raw["sequence_type"],
        sequence=raw.get("sequence", ""),
        min_len=raw["min_len"],
        max_len=raw["max_len"],
        onlist=_load_onlist(raw, spec_dir),
        regions=[_parse_region(child, spec_dir) for child in raw.get("regions") or []],
    )


def load_spec(path: Path | str = POOL_SPEC) -> Spec:
    """Reads a seqspec YAML file into a `Spec`."""
    path = Path(path)
    raw = yaml.safe_load(path.read_text())
    library = raw.get("library_spec") or []
    if len(library) != 1:
        raise ValueError(
            f"{path.name} has {len(library)} top-level regions; these specs carry exactly one "
            "(the `dna` modality)."
        )
    return Spec(
        name=raw.get("name", path.stem),
        assay_id=raw.get("assay_id", path.stem),
        description=raw.get("description", ""),
        root=_parse_region(library[0], path.parent),
        path=path,
    )


# --------------------------------------------------------------------------- #
# Matching an oligo against a spec
# --------------------------------------------------------------------------- #


def _leaf_pattern(region: Region) -> str:
    if region.sequence_type == "fixed":
        return re.escape(region.sequence)
    if region.sequence_type == "onlist":
        if not region.onlist:
            raise ValueError(f"Region {region.region_id!r} is `onlist` but lists no sequences.")
        return "(?:" + "|".join(re.escape(seq) for seq in region.onlist) + ")"
    if (constant := PREFIX_CONSTANTS.get(region.region_id)) is not None:
        # Longest first, so the regex prefers the longer filler when both a long
        # arm and a long filler could explain the same stretch.
        lengths = range(min(region.max_len, len(constant)), region.min_len - 1, -1)
        return "(?:" + "|".join(re.escape(constant[:n]) for n in lengths) + ")"
    alphabet = ALPHABETS.get(region.region_id, "ACGT")
    return f"[{alphabet}]{{{region.min_len},{region.max_len}}}"


def compile_pattern(region: Region) -> re.Pattern[str]:
    """
    Builds one anchored regex matching any oligo the region admits.

    Each leaf becomes a named group, so a match reports where every region
    landed. Adjacent variable regions are genuinely ambiguous - nothing in the
    sequence says where the splint's A/T pad stops and an A/T-initial arm
    begins - so the boundary a match reports between two such regions is one
    valid reading, not the only one. Their combined length is not ambiguous,
    which is what `COUPLED_LENGTHS` checks.
    """
    parts = [f"(?P<{leaf.region_id}>{_leaf_pattern(leaf)})" for leaf in region.leaves()]
    return re.compile("^" + "".join(parts) + "$", re.IGNORECASE)


def match_regions(sequence: str, region: Region) -> dict[str, str] | None:
    """Returns each leaf's slice of `sequence`, or None if it does not conform."""
    matched = compile_pattern(region).match(sequence.strip())
    return None if matched is None else matched.groupdict()


# --------------------------------------------------------------------------- #
# Validation
# --------------------------------------------------------------------------- #


def exact_length(oligo: Region) -> int | None:
    """
    The one length every oligo of this kind has, or None if it genuinely varies.

    seqspec sums each child's range independently, so `splint_oligo` advertises
    137-159 nt. Folding in `COUPLED_LENGTHS` collapses that to the 148 nt the
    synthesiser is actually asked for, which makes a length complaint useful
    rather than merely true.
    """
    coupling = COUPLED_LENGTHS.get(oligo.region_id)
    coupled = set(coupling[0]) if coupling else set()
    total = coupling[1] if coupling else 0
    for leaf in oligo.leaves():
        if leaf.region_id in coupled:
            continue
        if leaf.min_len != leaf.max_len:
            return None
        total += leaf.min_len
    return total


def _check_coupled(oligo_id: str, parts: dict[str, str]) -> list[str]:
    coupling = COUPLED_LENGTHS.get(oligo_id)
    if coupling is None:
        return []
    (first, second), total = coupling
    got = len(parts[first]) + len(parts[second])
    if got == total:
        return []
    return [
        (
            f"{first} + {second} is {got} nt, but the two always fill {total} nt between them "
            f"({first}={len(parts[first])}, {second}={len(parts[second])})"
        )
    ]


def working_probe(parts: dict[str, str], oligo: Region) -> str:
    """
    The part of an ordered oligo that survives the KpnI/BamHI double digest.

    Everything else - the 5' amplification primer, the 3' constant region, the
    splint's backfill - exists only to get the pool synthesized and amplified,
    and is cut away before the probe is used.
    """
    probe = oligo.find(f"{oligo.region_id.split('_')[0]}_probe")
    return "".join(parts[leaf.region_id] for leaf in probe.leaves())


def _check_restriction(oligo: Region, parts: dict[str, str]) -> list[str]:
    """
    The working probe must contain no site of the enzymes used to release it.

    Checked on the working probe alone: the handles that get cut away contain
    these sites on purpose, which is how the probe is released at all.
    """
    from Bio import Restriction, Seq

    body = Seq.Seq(working_probe(parts, oligo))
    kind = oligo.region_id.split("_")[0]
    return [
        f"the working {kind} contains a {name} site, so the digest that releases it would cut it in half"
        for name in SOLAR_RESTRICTION
        if getattr(Restriction, name).search(body)
    ]


def validate_oligo(sequence: str, oligo: Region) -> tuple[dict[str, str] | None, list[str]]:
    """Checks one oligo against one spec region. Returns its parts and any problems."""
    sequence = sequence.strip().upper()
    expected = exact_length(oligo)
    if expected is not None and len(sequence) != expected:
        return None, [f"is {len(sequence)} nt; every {oligo.region_id} is exactly {expected} nt"]

    parts = match_regions(sequence, oligo)
    if parts is None:
        return None, [
            (
                f"is the right length but does not match the {oligo.region_id} structure. "
                f"Expected, 5' to 3': {', '.join(leaf.region_id for leaf in oligo.leaves())}"
            )
        ]
    problems = _check_coupled(oligo.region_id, parts)
    problems += _check_restriction(oligo, parts)
    return parts, problems


@dataclass
class PoolReport:
    """What `validate_pool` found."""

    n_pairs: int
    problems: list[str]
    readouts: list[tuple[int, int, int]]

    @property
    def ok(self) -> bool:
        return not self.problems


def detect_bcidx(pool: Path | str) -> int | None:
    """
    Works out which `bcidx` a pool was built with, or None if none of them fit.

    `bcidx` only changes the primer binding sites, so validating against the
    wrong one fails on every oligo with a structure complaint that says nothing
    about the probes. Matching the first pair against each candidate spec turns
    that into either the right answer or one clear message. Indices are tried in
    order and the common case is 0, so this usually reads a single file.
    """
    lines = _read_pool(pool)
    if len(lines) < 2:
        return None
    for bcidx in available_bcidx():
        spec = load_spec(pool_spec_path(bcidx))
        if match_regions(lines[0], spec.root.find("splint_oligo")) and match_regions(
            lines[1], spec.root.find("padlock_oligo")
        ):
            return bcidx
    return None


def _read_pool(pool: Path | str) -> list[str]:
    return [line.strip() for line in Path(pool).read_text().splitlines() if line.strip()]


def validate_pool(
    pool: Path | str,
    spec: Spec | None = None,
    codebook: dict[str, list[int]] | None = None,
) -> PoolReport:
    """
    Validates an assembled oligo pool against the SOLAR spec.

    `pool` is a `<name>_final.txt` from `mkprobes assemble` - one oligo per
    line, splint and padlock alternating. Every structural rule comes from the
    spec; the pairing rules (splint and padlock alternate, and the splint's
    clamp templates its partner's ends) come from the chemistry, which seqspec
    has no way to express.
    """
    spec = spec or load_spec(POOL_SPEC)
    splint_region = spec.root.find("splint_oligo")
    padlock_region = spec.root.find("padlock_oligo")
    readout_ids = _readout_ids(spec)

    lines = _read_pool(pool)
    problems: list[str] = []
    if len(lines) % 2:
        problems.append(
            f"{len(lines)} oligos in the pool, which is odd. Every probe is a splint/padlock "
            "pair, so the count is always even."
        )

    readouts: list[tuple[int, int, int]] = []
    for index in range(0, len(lines) - 1, 2):
        splint, padlock = lines[index], lines[index + 1]
        pair = index // 2 + 1

        splint_parts, splint_problems = validate_oligo(splint, splint_region)
        padlock_parts, padlock_problems = validate_oligo(padlock, padlock_region)
        problems += [f"pair {pair} splint: {problem}" for problem in splint_problems]
        problems += [f"pair {pair} padlock: {problem}" for problem in padlock_problems]

        if splint_parts is None or padlock_parts is None:
            continue

        # The reason the pair exists: the splint's last 12 nt have to hold the
        # padlock's two ends together, or nothing circularises and the probe is
        # silent. Checked on the working probes, since that is what the clamp
        # acts on once the handles are digested away.
        if not test_splint_padlock(
            working_probe(splint_parts, splint_region),
            working_probe(padlock_parts, padlock_region),
            lengths=(6, 6),
        ):
            problems.append(
                f"pair {pair}: the splint's clamp does not template the padlock's ends, so this "
                "padlock cannot circularise"
            )

        codes = tuple(readout_ids[padlock_parts[f"padlock_readout_{n}"].upper()] for n in (1, 2, 3))
        readouts.append(codes)  # type: ignore[arg-type]
        if len(set(codes)) != 3:
            problems.append(f"pair {pair}: readouts {codes} repeat a bit; a codeword uses three distinct bits")

    if codebook is not None:
        problems += _check_against_codebook(readouts, codebook)

    return PoolReport(n_pairs=len(lines) // 2, problems=problems, readouts=readouts)


def _readout_ids(spec: Spec) -> dict[str, int]:
    """Maps each readout sequence back to its 1-based ID in the vendored table."""
    onlist = spec.root.find("padlock_readout_1").onlist
    return {seq: idx for idx, seq in enumerate(onlist, start=1)}


def _check_against_codebook(
    readouts: list[tuple[int, int, int]], codebook: dict[str, list[int]]
) -> list[str]:
    """Every codeword in the pool has to be one the codebook actually assigned."""
    valid = {tuple(sorted(bits)) for name, bits in codebook.items() if not name.startswith("Blank-")}
    blanks = {tuple(sorted(bits)) for name, bits in codebook.items() if name.startswith("Blank-")}
    problems = []
    for codes in {tuple(sorted(c)) for c in readouts}:
        if codes in blanks:
            problems.append(
                f"codeword {codes} is a Blank in the codebook. Blanks measure the false-positive "
                "rate and must never be synthesised, so this pool would destroy that estimate."
            )
        elif codes not in valid:
            problems.append(f"codeword {codes} is in the pool but is not assigned to any target in the codebook")
    return problems


# --------------------------------------------------------------------------- #
# Drawing
# --------------------------------------------------------------------------- #

#: One colour per region_type, so the same kind of part reads the same across
#: both diagrams. Chosen to stay distinguishable in greyscale print.
REGION_COLORS: dict[str, str] = {
    "custom_primer": "#b0bec5",
    "linker": "#eceff1",
    "cdna": "#4db6ac",
    "barcode": "#7e57c2",
    "named": "#cfd8dc",
    "dna": "#cfd8dc",
}


#: Regions read on the opposite strand: the reverse primer is their reverse
#: complement and extends back toward the 5' end. seqspec has no strand field,
#: so the direction lives here and in the region descriptions.
REVERSE_STRAND_REGIONS = frozenset({"splint_primer_3p", "padlock_primer_3p"})

#: Regions with no orientation to draw. Both recognition sites are palindromes -
#: GGATCC and GGTACC are their own reverse complements - so they read the same
#: on either strand, and the enzyme cuts the duplex rather than a strand. An
#: arrow would assert a direction the feature does not have, so these are drawn
#: as plain boxes.
UNDIRECTED_REGIONS = frozenset({
    "splint_kpni_site",
    "padlock_kpni_site",
    "splint_bamhi_site",
    "padlock_bamhi_site",
})


def _strand(region: Region) -> int:
    if region.region_id in UNDIRECTED_REGIONS:
        return 0
    return -1 if region.region_id in REVERSE_STRAND_REGIONS else +1


#: Overrides `REGION_COLORS` for particular regions. seqspec has no region_type
#: for a restriction site, so the sites are typed `linker` like every other
#: piece of constant machinery - but they are the one thing on the oligo that
#: says where it gets cut, so they are worth picking out by eye.
REGION_ID_COLORS: dict[str, str] = {
    "splint_kpni_site": "#ef9a9a",
    "padlock_kpni_site": "#ef9a9a",
    "splint_bamhi_site": "#ef9a9a",
    "padlock_bamhi_site": "#ef9a9a",
    "splint_kpni_retained": "#ffcc80",
    "padlock_kpni_retained": "#ffcc80",
}


def _color(region: Region) -> str:
    return REGION_ID_COLORS.get(region.region_id) or REGION_COLORS.get(region.region_type, "#cfd8dc")


def _feature_label(region: Region) -> str:
    if region.min_len == region.max_len:
        return f"{region.region_id} ({region.min_len})"
    return f"{region.region_id} ({region.min_len}-{region.max_len})"


def representative_widths(molecule: Region) -> dict[str, int]:
    """
    A width per region for drawing: one plausible oligo rather than a worst case.

    Drawing every variable region at its maximum would total 159 nt for a
    splint that is always 148, and would imply the arm and its filler can both
    be long at once when in fact one shortens as the other grows. Coupled
    regions therefore split their shared budget in proportion to how much each
    can vary, which keeps both visible and makes the figure add up to the
    length actually synthesized.
    """
    widths = {leaf.region_id: leaf.max_len for leaf in molecule.leaves()}
    coupling = COUPLED_LENGTHS.get(molecule.region_id)
    if coupling is None:
        return widths

    (first_id, second_id), total = coupling
    first, second = molecule.find(first_id), molecule.find(second_id)
    spare = total - first.min_len - second.min_len
    give = round(spare * (first.max_len - first.min_len) / ((first.max_len - first.min_len) + (second.max_len - second.min_len)))
    widths[first_id] = first.min_len + give
    widths[second_id] = total - widths[first_id]
    return widths


def _records(spec: Spec):
    """One DnaFeaturesViewer record per separately synthesized molecule."""
    try:
        from dna_features_viewer import GraphicFeature, GraphicRecord
    except ImportError as exc:  # pragma: no cover - depends on the install
        raise click.ClickException(
            "Drawing needs DnaFeaturesViewer, which is an optional extra.\n"
            "Install it with:  pip install 'mkprobes[viz]'   (or: pip install dna_features_viewer)"
        ) from exc

    # The top-level `dna` region nests one child per molecule when there are
    # several (splint and padlock); a single-molecule spec is drawn on its own.
    molecules = [r for r in spec.root.regions if not r.is_leaf] or [spec.root]
    records = []
    for molecule in molecules:
        widths = representative_widths(molecule)
        features, cursor = [], 0
        for leaf in molecule.leaves():
            width = widths[leaf.region_id]
            features.append(
                GraphicFeature(
                    start=cursor,
                    end=cursor + width,
                    strand=_strand(leaf),
                    color=_color(leaf),
                    label=_feature_label(leaf),
                    linewidth=0.8,
                )
            )
            cursor += width
        records.append((molecule, GraphicRecord(sequence_length=cursor, features=features)))
    return records


def draw_spec(spec: Spec, output: Path, width: float = 14.0, dpi: int = 200) -> Path:
    """Renders every molecule in `spec` as a stacked figure of labelled regions."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    records = _records(spec)
    # Labels stack above the molecule, so a 15-region padlock needs several
    # times the height of a single-region readout probe. A fixed height per
    # panel either crops the busy one or strands the simple one in whitespace.
    heights = [1.0 + 0.24 * len(record.features) for _, record in records]
    # The orientation note gets its own (empty) row rather than a `fig.text`,
    # which constrained_layout does not reserve space for - it lands on top of
    # the ruler instead.
    note = any(leaf.region_id in REVERSE_STRAND_REGIONS for leaf in spec.root.leaves())
    if note:
        heights.append(0.6)
    fig, axes = plt.subplots(
        len(heights),
        1,
        figsize=(width, sum(heights) + 0.6),
        height_ratios=heights,
        constrained_layout=True,
    )
    axes = [axes] if len(heights) == 1 else list(axes)
    if note:
        caption = axes.pop()
        caption.axis("off")
        caption.text(
            0,
            0.5,
            "Arrows run 5'\u2192 3'. The 3' primer sites point left: the reverse primer is "
            "their reverse complement and primes back toward the 5' end.\n"
            "Red marks a restriction site, drawn square-ended because GGATCC and GGTACC are "
            "palindromes - they read the same either way, and the enzyme cuts the duplex.\n"
            "Each site straddles its cut: the five bases shown go with the discarded handle, "
            "and only the single amber base stays on the probe.",
            fontsize=8,
            va="center",
            linespacing=1.6,
            color="#546e7a",
        )

    for ax, (molecule, record) in zip(axes, records):
        record.plot(ax=ax, with_ruler=True, annotate_inline=False)
        # Each label carries its region's permitted range, so the caption only
        # needs to say what the whole molecule comes to.
        exact = exact_length(molecule)
        span = f"always {exact} nt" if exact else f"{molecule.min_len}-{molecule.max_len} nt"
        ax.set_title(f"{molecule.name}  [{span}]", fontsize=10, loc="left")

    fig.suptitle(spec.name, fontsize=12, x=0.01, ha="left")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #


def _resolve_spec(spec: Path | None, which: str, bcidx: int | None) -> Spec:
    if spec is not None:
        return load_spec(spec)
    if which == "readout":
        return load_spec(READOUT_SPEC)
    return load_spec(pool_spec_path(bcidx or 0))


# fmt: off
@click.command("validate-pool")
@click.argument("pool", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option("--codebook", "-c", type=click.Path(exists=True, dir_okay=False, path_type=Path), help="Also check that every codeword in the pool is one this codebook assigns.")
@click.option("--bcidx", type=int, default=None, help="Header/footer pair the panel was built with (the manifest field). Detected from the pool when omitted.")
@click.option("--spec", type=click.Path(exists=True, dir_okay=False, path_type=Path), default=None, help="Validate against this seqspec file instead of a vendored one.")
@click.option("--limit", type=int, default=20, show_default=True, help="Stop listing problems after this many.")
# fmt: on
def validate_pool_cli(pool: Path, codebook: Path | None, bcidx: int | None, spec: Path | None, limit: int):
    """Check an assembled oligo pool against the SOLAR seqspec.

    POOL is a `<name>_final.txt` from `mkprobes assemble`: one oligo per line,
    splint and padlock alternating. Every oligo is checked for the right
    regions in the right order, the right total length, no KpnI or BamHI site
    inside the working probe, and a clamp that can actually circularise its
    partner. With `--codebook`, the readouts on each padlock are checked
    against the codewords that codebook assigns.

    `bcidx` is worked out from the pool itself, so you only need to pass it to
    assert that a pool is the index you meant it to be.
    """
    if spec is None and bcidx is None:
        # Without this, a bcidx-7 pool checked against the bcidx-0 spec fails on
        # every oligo, and the message blames the probes rather than the index.
        bcidx = detect_bcidx(pool)
        if bcidx is None:
            raise click.ClickException(
                f"{pool.name}: could not tell which bcidx this pool was built with - its first "
                "probe pair matches none of them. Either the pool is not a SOLAR oligo pool, or "
                "it was built against a header/footer table this version does not ship. Pass "
                "--bcidx to check against one anyway and see what differs."
            )
        logger.info(f"Detected bcidx {bcidx}.")

    try:
        loaded = _resolve_spec(spec, "pool", bcidx)
    except FileNotFoundError as exc:
        raise click.BadParameter(str(exc), param_hint="--bcidx") from exc

    book = json.loads(codebook.read_text()) if codebook else None
    report = validate_pool(pool, loaded, book)

    if report.ok:
        logger.info(f"{pool.name}: {report.n_pairs} probe pairs, all conform to {loaded.name}.")
        return

    for problem in report.problems[:limit]:
        logger.error(problem)
    if len(report.problems) > limit:
        logger.error(f"... and {len(report.problems) - limit} more (raise --limit to see them).")
    raise click.ClickException(
        f"{pool.name}: {len(report.problems)} problem(s) across {report.n_pairs} probe pairs, "
        f"checked against {loaded.name}. This pool does not match the SOLAR spec - do not order it."
    )


# fmt: off
@click.command("draw-spec")
@click.option("--which", type=click.Choice(("pool", "readout")), default="pool", show_default=True, help="Which vendored spec to draw.")
@click.option("--bcidx", type=int, default=0, show_default=True, help="Which pool spec to draw. Only the primer binding sites differ between indices.")
@click.option("--spec", type=click.Path(exists=True, dir_okay=False, path_type=Path), default=None, help="Draw this seqspec file instead of a vendored one.")
@click.option("--out", "-o", type=click.Path(dir_okay=False, path_type=Path), default=None, help="Where to write the figure. Defaults to <assay_id>.png in the working directory.")
@click.option("--dpi", type=int, default=200, show_default=True)
# fmt: on
def draw_spec_cli(which: str, bcidx: int, spec: Path | None, out: Path | None, dpi: int):
    """Draw a seqspec as a labelled diagram of its regions.

    Produces a generalized view: every region is drawn at the length it would
    plausibly have and labelled with its permitted range, so the picture
    describes the whole library rather than any one oligo. Primer binding sites
    read on the opposite strand are drawn pointing the other way. Needs the
    `viz` extra (`pip install 'mkprobes[viz]'`).
    """
    try:
        loaded = _resolve_spec(spec, which, bcidx)
    except FileNotFoundError as exc:
        raise click.BadParameter(str(exc), param_hint="--bcidx") from exc

    target = out or Path(f"{loaded.assay_id}.png")
    draw_spec(loaded, target, dpi=dpi)
    logger.info(f"Wrote {target}")
