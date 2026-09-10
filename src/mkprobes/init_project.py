"""
Project scaffolding: ``mkprobes init``.

The documented workflow used to end at a hand-authored `manifest.json` that no
walkthrough showed, whose `bcidx` field indexes rows of an internal table, and
whose `n_probes` accepts two magic strings. A user could complete every design
step and only then discover their manifest was wrong.

`init` writes a project that is already valid: a commented target list, a
manifest stub with real values and an explanation of each field, and a README
listing the commands to run in order. `check_manifest` validates one before
assembly spends hours proving it wrong.
"""

import json
from itertools import chain
from pathlib import Path
from typing import Any

import rich_click as click
from loguru import logger
from pydantic import TypeAdapter, ValidationError

from .assembly import hfs
from .codebook.codebook import ProbeSet
from .codebook.generate import ORDER
from .constants import GOOD_SPECIES, SOLAR_RESTRICTION
from .design import DesignParameters

GENES_TEMPLATE = """\
# One target per line: a gene name, or a transcript ID for a custom dataset.
#
# Blank lines are ignored, and everything after a '#' is a comment, so you can
# note why a target is here. Order does not matter.
#
# Replace the examples below with your own targets, then run the commands in
# README.md in order.

Sox2
Pax6   # dorsal telencephalon marker
"""

README_TEMPLATE = """\
# {name}

A SOLAR probe design project, created by `mkprobes init`.

## Before you start

You need a dataset for your species. If you have not built one yet:

```bash
mkprobes prepare {dataset_parent} --species mouse     # mouse or human
mkprobes ingest {dataset} --genome genome.fa --gtf annotation.gtf --species <name>
```

See the {docs_before} page for what to download and what to expect.

## Steps

Run these from this directory, in order. Each one is a single command, and each
checks the previous step's output before it starts.

```bash
# 1. Resolve your target names, then pick one transcript per gene
mkprobes chkgenes {dataset} genes.txt
mkprobes convert-to-transcripts {dataset} genes.converted.txt{transcript_mode}

# 2. Assign readout bits to each target, starting at the manifest's offset
mkprobes make-codebook {dataset} genes.converted.tss.txt -o codebook.json

# 3. Design probes for every target (this is the long one)
mkprobes run-panel {dataset} codebook.json

# 4. Check how many probes each target ended up with
mkprobes filter-genes output --genes genes.converted.tss.txt --min-probes 48

# 5. Triage targets that came up short, then re-run step 3 to apply the result
mkprobes assemble manifest.json short 12

# 6. Build the orderable oligo pool
mkprobes assemble manifest.json gen
```

The pool lands in `generated/`, alongside a `.provenance.json` recording the
version, dataset and parameters that produced it. `mkprobes provenance <file>`
prints the same record from any output.

## Files here

| File | What it is |
| --- | --- |
| `genes.txt` | Your targets, one per line. Edit this first. |
| `manifest.json` | Describes this panel: design settings for step 3, assembly fields for step 6. |
| `codebook.json` | Written by step 2. Do not edit by hand. |
| `output/` | Per-target design output, written by step 3. |
| `generated/` | The orderable pool, written by step 6. |
"""

MANIFEST_COMMENT = {
    "name": "Names the output files in generated/.",
    "species": "Used to pick the RepeatMasker taxon; any name is accepted.",
    "codebook": "Written by `mkprobes make-codebook`, relative to this file.",
    "bcidx": (
        "Which header/footer pair to build against. Each index uses two rows of the "
        "internal table, so valid values are 0 to {max_bcidx}. Use a different index "
        "for each panel you will pool together."
    ),
    "offset": (
        "Where this panel's codebook starts in the readout order, as a 0-based bit "
        "position (not a readout ID). `mkprobes make-codebook` reads it from here. "
        "Panels hybridised together must not share bits: leave the first panel at 0 "
        "and give each further panel the number of bits already taken (a 10-bit "
        "panel takes 10, so the next starts at 10). {n_readouts} readout IDs exist "
        "in all. Pair it with a distinct bcidx."
    ),
    "n_probes": (
        'Maximum probes per target in the pool. A number, or "high" (34) or '
        '"low" (16). Omit to let the species decide.'
    ),
    "design": (
        "How probes are designed. `mkprobes run-panel` reads these from here, and "
        "assembly warns if the outputs were designed under anything else. The values "
        "written are the defaults; delete a field to keep its default. They change how "
        "probes hybridise, so read the design_probes guide before editing."
    ),
    "design.tm_range": (
        "Tm window in °C (at the design formamide) a probe must fall in. Lower the first "
        "number on AT-rich transcripts to admit probes that bind less tightly."
    ),
    "design.length_range": (
        "Probe length window in nt. Raise the second number (60 at most) so AT-rich "
        "windows can grow long enough to reach the Tm floor."
    ),
    "design.split_tm": (
        "Tm in °C each arm of the split probe must reach. The largest lever on AT-rich "
        "transcripts, and the one to lower with the most care: both arms must bind."
    ),
    "design.min_probes": "Probes per target the screen aims for.",
    "design.max_overlap": (
        "How far neighbouring probes may overlap (nt, multiples of 5) to reach min_probes. "
        "0 keeps probes disjoint."
    ),
}


def max_bcidx() -> int:
    """Highest usable `bcidx`: each one consumes two header/footer rows."""
    return len(hfs) // 2 - 1


def manifest_stub(
    name: str, species: str, bcidx: int = 0, n_probes: int = 24, offset: int = 0
) -> list[dict[str, Any]]:
    """
    A manifest that is valid on the first try.

    The design block is written out in full, with the defaults for the species,
    so the numbers a user might need to change are in front of them rather
    than buried in code. Writing them explicitly designs exactly what leaving
    them out would.
    """
    design = DesignParameters().resolve(reference=species in GOOD_SPECIES)
    return [
        {
            "_comment": {
                key: value.format(max_bcidx=max_bcidx(), n_readouts=len(ORDER))
                for key, value in MANIFEST_COMMENT.items()
            },
            "name": name,
            "species": species,
            "codebook": "codebook.json",
            "bcidx": bcidx,
            "offset": offset,
            "n_probes": n_probes,
            "design": design.model_dump(mode="json"),
        }
    ]


def check_manifest(path: Path) -> list[ProbeSet]:
    """
    Validates a manifest, checking the things Pydantic cannot.

    Raises `ValueError` naming the problem and the fix.
    """
    try:
        probesets = TypeAdapter(list[ProbeSet]).validate_json(path.read_text())
    except ValidationError as e:
        raise ValueError(f"{path} is not a valid manifest:\n{e}") from e

    if not probesets:
        raise ValueError(f"{path} describes no probe sets. It needs at least one entry.")

    limit = max_bcidx()
    for probeset in probesets:
        if not 0 <= probeset.bcidx <= limit:
            raise ValueError(
                f"{path}: probe set {probeset.name!r} has bcidx {probeset.bcidx}, but only "
                f"0 to {limit} exist. Each index uses two rows of the header/footer table."
            )
        if not 0 <= probeset.offset < len(ORDER):
            raise ValueError(
                f"{path}: probe set {probeset.name!r} has offset {probeset.offset}, but only "
                f"{len(ORDER)} readout IDs exist, so offsets run from 0 to {len(ORDER) - 1}."
            )
        codebook = path.parent / probeset.codebook
        if not codebook.exists():
            raise ValueError(
                f"{path}: probe set {probeset.name!r} refers to {probeset.codebook}, which does "
                f"not exist. Run `mkprobes make-codebook` first, or correct the path."
            )
        check_codebook_offset(path, probeset)

    names = [p.name for p in probesets]
    if len(set(names)) != len(names):
        raise ValueError(f"{path}: probe set names must be unique, got {names}.")
    check_disjoint_bits(path, probesets)
    return probesets


def codebook_bits(path: Path, probeset: ProbeSet) -> set[int]:
    """Every readout bit a probe set's codebook occupies, Blanks included."""
    codebook = probeset.load_codebook(path.parent, include_blank=True)
    return set(chain.from_iterable(codebook.values()))


def check_codebook_offset(path: Path, probeset: ProbeSet) -> None:
    """
    The codebook has to start where the manifest says it does.

    A codebook generated before the offset was edited, or with `--offset`
    overriding the manifest, would otherwise sail through to a pooled hyb
    with the wrong bits. The first position in the readout order that the
    codebook uses must be the manifest's offset.
    """
    bits = codebook_bits(path, probeset)
    if not bits:
        return
    if unknown := sorted(bits - set(ORDER)):
        raise ValueError(
            f"{path}: probe set {probeset.name!r}'s codebook {probeset.codebook} uses readout "
            f"bit(s) {unknown}, but only 1 to {len(ORDER)} exist."
        )
    first = min(ORDER.index(bit) for bit in bits)
    if first != probeset.offset:
        raise ValueError(
            f"{path}: probe set {probeset.name!r} has offset {probeset.offset}, but its codebook "
            f"{probeset.codebook} starts at bit position {first}. Re-run `mkprobes make-codebook` "
            f"so the codebook follows the manifest, or set offset to {first}."
        )


def check_disjoint_bits(path: Path, probesets: list[ProbeSet]) -> None:
    """
    Panels in one manifest are panels meant to be pooled, so their codebooks
    must occupy disjoint readout bits. Two probe sets naming the same codebook
    file are one panel described twice (as `run-panel` allows) and are not
    compared.
    """
    seen: dict[str, tuple[str, set[int]]] = {}
    for probeset in probesets:
        key = str((path.parent / probeset.codebook).resolve())
        if key in seen:
            continue
        bits = codebook_bits(path, probeset)
        for other_name, other_bits in seen.values():
            if shared := sorted(bits & other_bits):
                raise ValueError(
                    f"{path}: probe sets {other_name!r} and {probeset.name!r} share readout "
                    f"bit(s) {shared}, so they cannot be hybridised together. Give one a "
                    "different offset and re-run `mkprobes make-codebook` for it."
                )
        seen[key] = (probeset.name, bits)


@click.command("init")
@click.argument("project", type=click.Path(file_okay=False, path_type=Path))
@click.option("--species", default="mouse", show_default=True, help="Species name recorded in the manifest.")
@click.option(
    "--dataset",
    type=click.Path(path_type=Path),
    default=None,
    help="Path to the dataset to design against. Used in the generated README.",
)
@click.option("--bcidx", type=int, default=0, show_default=True, help="Header/footer pair for this panel.")
@click.option("--offset", type=int, default=0, show_default=True,
              help="Bit position this panel's codebook starts at. 0 for a first panel; the number of "
              "bits already taken for a panel pooled with earlier ones.")
@click.option("--force", is_flag=True, help="Overwrite files that already exist.")
def init(project: Path, species: str, dataset: Path | None, bcidx: int, offset: int, force: bool):
    """Create a probe design project, ready to run.

    Writes a commented target list, a valid manifest, and a README listing the
    commands to run in order. Edit `genes.txt`, then follow the README.
    """
    limit = max_bcidx()
    if not 0 <= bcidx <= limit:
        raise click.BadParameter(
            f"only 0 to {limit} exist; each index uses two rows of the header/footer table.",
            param_hint="--bcidx",
        )
    if not 0 <= offset < len(ORDER):
        raise click.BadParameter(
            f"only {len(ORDER)} readout IDs exist, so offsets run from 0 to {len(ORDER) - 1}.",
            param_hint="--offset",
        )

    dataset_path = dataset or Path("../data") / species
    files = {
        "genes.txt": GENES_TEMPLATE,
        "manifest.json": json.dumps(manifest_stub(project.name, species, bcidx, offset=offset), indent=2) + "\n",
        "README.md": README_TEMPLATE.format(
            name=project.name,
            dataset=dataset_path,
            dataset_parent=dataset_path.parent,
            docs_before="`before_you_start`",
            # Reference datasets pick the canonical isoform from Ensembl; custom
            # ones have no such annotation and fall back to the longest.
            transcript_mode="" if species in GOOD_SPECIES else " -m longest",
        ),
    }

    existing = [name for name in files if (project / name).exists()]
    if existing and not force:
        raise click.ClickException(
            f"{project} already contains {', '.join(existing)}. "
            "Pass --force to overwrite, or choose another directory."
        )

    project.mkdir(parents=True, exist_ok=True)
    for name, content in files.items():
        (project / name).write_text(content)

    logger.info(f"Created {project}/ with {', '.join(files)}.")
    click.echo(
        f"\nProject ready at {project}/\n"
        f"  1. Edit {project / 'genes.txt'} - one target per line\n"
        f"  2. Follow the steps in {project / 'README.md'}\n"
    )


@click.command("check-manifest")
@click.argument("manifest", type=click.Path(exists=True, dir_okay=False, path_type=Path))
def check_manifest_cli(manifest: Path):
    """Check a manifest before assembly spends hours proving it wrong."""
    try:
        probesets = check_manifest(manifest)
    except ValueError as e:
        raise click.ClickException(str(e)) from e
    for probeset in probesets:
        click.echo(
            f"{probeset.name}: {probeset.species}, codebook {probeset.codebook}, "
            f"bcidx {probeset.bcidx}, offset {probeset.offset}, enzymes {'+'.join(SOLAR_RESTRICTION)}, "
            f"design {probeset.design.describe()}"
        )
    click.echo(f"{manifest} is valid.")
