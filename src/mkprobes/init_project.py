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
from pathlib import Path
from typing import Any

import rich_click as click
from loguru import logger
from pydantic import TypeAdapter, ValidationError

from .assembly import hfs
from .codebook.codebook import ProbeSet
from .constants import GOOD_SPECIES, SOLAR_RESTRICTION

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

# 2. Assign readout bits to each target
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
| `manifest.json` | Describes this panel for the assembly step. |
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
    "n_probes": (
        'Maximum probes per target in the pool. A number, or "high" (34) or '
        '"low" (16). Omit to let the species decide.'
    ),
}


def max_bcidx() -> int:
    """Highest usable `bcidx`: each one consumes two header/footer rows."""
    return len(hfs) // 2 - 1


def manifest_stub(name: str, species: str, bcidx: int = 0, n_probes: int = 24) -> list[dict[str, Any]]:
    """A manifest that is valid on the first try."""
    return [
        {
            "_comment": {key: value.format(max_bcidx=max_bcidx()) for key, value in MANIFEST_COMMENT.items()},
            "name": name,
            "species": species,
            "codebook": "codebook.json",
            "bcidx": bcidx,
            "n_probes": n_probes,
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
        codebook = path.parent / probeset.codebook
        if not codebook.exists():
            raise ValueError(
                f"{path}: probe set {probeset.name!r} refers to {probeset.codebook}, which does "
                f"not exist. Run `mkprobes make-codebook` first, or correct the path."
            )

    names = [p.name for p in probesets]
    if len(set(names)) != len(names):
        raise ValueError(f"{path}: probe set names must be unique, got {names}.")
    return probesets


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
@click.option("--force", is_flag=True, help="Overwrite files that already exist.")
def init(project: Path, species: str, dataset: Path | None, bcidx: int, force: bool):
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

    dataset_path = dataset or Path("../data") / species
    files = {
        "genes.txt": GENES_TEMPLATE,
        "manifest.json": json.dumps(manifest_stub(project.name, species, bcidx), indent=2) + "\n",
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
            f"bcidx {probeset.bcidx}, enzymes {'+'.join(SOLAR_RESTRICTION)}"
        )
    click.echo(f"{manifest} is valid.")
