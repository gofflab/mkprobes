"""
Panel-level probe design driver: candidates -> screen -> construct for every
target in a codebook, in parallel.

This is the package integration of the lab's batch drivers
(scripts/probegen/1_run_codebook.py and 1_run_codebook_generic.py) and keeps
their flow: the codebook is the work list (Blank-* entries excluded), each
gene runs the full three-stage pipeline in its own process with a per-gene
log file, finished genes are skipped on re-runs, an `.acceptable.json`
allow-list feeds `candidates --allow` (and forces a re-screen/re-construct
for those genes), and failures are collected into `<codebook>.failed.txt`
without stopping the rest of the panel.

Exposed on the CLI as ``mkprobes run-panel``.
"""

import json
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool
from multiprocessing import get_context
from pathlib import Path

import polars as pl
import rich_click as click
from loguru import logger
from rich.progress import Progress

from .design import (
    DEFAULT_MAX_OVERLAP,
    DEFAULT_MIN_PROBES,
    DesignParameters,
    design_from_flags,
    design_options,
    matches_recorded,
)
from .utils.provenance import read_provenance

# Production defaults, inherited from the original batch drivers.
DEFAULT_RESTRICTION = ("BamHI", "KpnI")
DEFAULT_MINIMUM = DEFAULT_MIN_PROBES
DEFAULT_MAXOVERLAP = DEFAULT_MAX_OVERLAP
DEFAULT_TARGET_PROBES = 48
DEFAULT_WORKERS = 16


def final_parquet(output: Path, gene: str, bits: list[int], restriction: tuple[str, ...]) -> Path:
    """Path of the construct output for a gene, as written by `construct`."""
    return output / f"{gene}_final_{''.join(restriction)}_{','.join(map(str, sorted(bits)))}.parquet"


def load_worklist(codebook_path: Path) -> dict[str, list[int]]:
    """
    Loads a codebook and returns the gene work list (Blank-* excluded).

    Duplicate keys are detected while the JSON is still a list of pairs: by the
    time it is a dict, a repeated target has already silently taken the last
    value, and every later stage would design against the wrong bits.
    """

    def reject_duplicates(pairs: list[tuple[str, object]]) -> dict[str, object]:
        seen: set[str] = set()
        duplicated = sorted({key for key, _ in pairs if key in seen or seen.add(key)})
        if duplicated:
            raise ValueError(
                f"{codebook_path} lists {len(duplicated)} target(s) more than once: "
                f"{', '.join(duplicated)}. Each target needs exactly one set of bits."
            )
        return dict(pairs)

    codebook = json.loads(codebook_path.read_text(), object_pairs_hook=reject_duplicates)
    return {k: v for k, v in codebook.items() if not k.startswith("Blank")}


def load_acceptable(codebook_path: Path, allow_file: Path | None) -> dict[str, list[str]]:
    """
    Loads the per-gene allow-list: `--allow-file` if given, else the
    `<codebook>.acceptable.json` convention (written by the manifest-assembly
    off-target triage).
    """
    path = allow_file or codebook_path.with_suffix(".acceptable.json")
    return json.loads(path.read_text()) if path.exists() else {}


def find_missing_final(codebook: dict[str, list[int]], output: Path, restriction: tuple[str, ...]) -> list[str]:
    """Genes whose final construct output is missing."""
    return [g for g, bits in sorted(codebook.items()) if not final_parquet(output, g, bits, restriction).exists()]


def designed_under(final: Path) -> dict | None:
    """The design settings a construct output records, or `None` if it records none."""
    return (read_provenance(final) or {}).get("design")


def design_for_codebook(codebook: Path, manifest: Path | None) -> tuple[DesignParameters, Path | None]:
    """
    The design block of the probe set that names this codebook.

    The manifest is where a panel's design settings live, so `run-panel` reads
    them from there: the `--manifest` given, or else a `manifest.json` beside
    the codebook, which is where `mkprobes init` puts it. An explicit manifest
    that does not list the codebook is an error; an implicit one is ignored
    with a warning, since it may describe a different panel altogether.

    Returns the settings and the manifest they came from (`None` if none).
    """
    from pydantic import ValidationError

    from .codebook.codebook import ProbeSet

    explicit = manifest is not None
    manifest = manifest or codebook.parent / "manifest.json"
    if not manifest.exists():
        return DesignParameters(), None

    try:
        probesets = ProbeSet.from_manifest(manifest)
    except ValidationError as e:
        raise ValueError(f"{manifest} is not a valid manifest:\n{e}") from e

    matches = [ps for ps in probesets if (manifest.parent / ps.codebook).resolve() == codebook.resolve()]
    if not matches:
        if explicit:
            raise ValueError(
                f"{manifest} has no probe set whose codebook is {codebook}. Add one, or point "
                "--manifest at the manifest that describes this panel."
            )
        logger.warning(
            f"{manifest} sits beside {codebook.name} but none of its probe sets name it, so its "
            "design settings are not used. Pass --manifest to use a manifest from elsewhere."
        )
        return DesignParameters(), None

    designs = {json.dumps(ps.design.explicit(), sort_keys=True) for ps in matches}
    if len(designs) > 1:
        raise ValueError(
            f"{manifest} lists {codebook.name} under {len(matches)} probe sets with different "
            "design settings, so it is ambiguous which to design under. Make them agree."
        )
    return matches[0].design, manifest


def run_gene(
    dataset_path: Path,
    output: Path,
    codebook: dict[str, list[int]],
    gene: str,
    acceptable: list[str] | None,
    overwrite: bool = False,
    log_level: str = "DEBUG",
    restriction: tuple[str, ...] = DEFAULT_RESTRICTION,
    target_probes: int = DEFAULT_TARGET_PROBES,
    codebook_hash: str | None = None,
    design: DesignParameters | None = None,
    **kwargs,
):
    """
    Runs candidates -> screen -> construct for one gene.

    Runs in a worker process: logs to `output/<gene>.log`, skips finished
    genes, and reuses an existing `<gene>_crawled.parquet` unless
    overwriting. An allow-list forces a re-screen/re-construct so accepted
    off-targets take effect, and so do candidates that were tiled under other
    design settings. Exceptions are re-raised tagged with the gene.
    """
    # Deferred imports keep worker startup (forkserver) lean.
    from .candidates import candidates_match_design, get_candidates
    from .codebook.finalconstruct import construct
    from .ext.dataset import ReferenceDataset, load_dataset
    from .screen import run_screen

    design = design or DesignParameters()

    logger.remove()
    logger.add(sys.stderr, level=log_level)
    logger.add(output / f"{gene}.log", level=log_level, colorize=False, backtrace=True, diagnose=True)

    if final_parquet(output, gene, codebook[gene], restriction).exists() and not overwrite:
        return

    ds = load_dataset(dataset_path)
    resolved = design.resolve(reference=isinstance(ds, ReferenceDataset))
    try:
        crawled = output / f"{gene}_crawled.parquet"
        if crawled.exists() and not overwrite and not candidates_match_design(output, gene, resolved):
            # Reusing them would make a changed setting do nothing, silently.
            logger.info(f"{gene}: existing candidates were tiled under other design settings; redoing.")
            overwrite = True
        if overwrite or not crawled.exists():
            get_candidates(
                ds,
                transcript=gene,
                output=output,
                ignore_revcomp=False,
                allow=acceptable,
                overwrite=overwrite,
                design=design,
                **kwargs,
            )
            time.sleep(1)  # let parquet writes settle before the next stage reads them
        overwrite = overwrite or acceptable is not None
        run_screen(
            output,
            gene,
            minimum=resolved.min_probes,
            restriction=list(restriction),
            maxoverlap=resolved.max_overlap,
            overwrite=overwrite,
        )
        # No overlap given: construct takes the screened file at the overlap
        # the search settled on, so `max_overlap` reaches the pool.
        construct(
            ds,
            output,
            transcript=gene,
            codebook=codebook,
            restriction=list(restriction),
            target_probes=target_probes,
            codebook_hash=codebook_hash,
            overwrite=overwrite,
            design=resolved,
        )
    except (Exception, SystemExit) as e:
        # Keep the cause in the message: the panel driver reports this line to a
        # user who cannot see the worker's traceback, and "No probes left after
        # filtering" is the actionable part, not the gene name.
        #
        # `SystemExit` is included because a process exit inside a stage is not
        # an `Exception`: it used to pass through the driver's error handling
        # and end the whole panel, silently, with nothing recorded.
        raise RuntimeError(f"{gene}: {type(e).__name__}: {e}") from e


def run_panel(
    dataset_path: Path,
    codebook_path: Path,
    output: Path,
    *,
    gene: str | None = None,
    allow_file: Path | None = None,
    workers: int = DEFAULT_WORKERS,
    overwrite: bool = False,
    restriction: tuple[str, ...] = DEFAULT_RESTRICTION,
    target_probes: int = DEFAULT_TARGET_PROBES,
    design: DesignParameters | None = None,
) -> dict[str, list[str]]:
    """
    Designs probes for every gene in the codebook, in parallel.

    `design` holds the thermodynamic and tiling settings; unset fields are
    the defaults. Returns {"done": [...], "skipped": [...], "failed": [...]}.
    Failed genes are also appended to `<codebook>.failed.txt` (recreated per
    run).
    """
    from .codebook.codebook import hash_codebook_file

    design = design or DesignParameters()
    codebook = load_worklist(codebook_path)
    # From the file, not the worklist: load_worklist drops Blank codes, and
    # hashing the filtered dict yields a different value from the one
    # make-codebook reported.
    codebook_hash = hash_codebook_file(codebook_path)
    acceptable = load_acceptable(codebook_path, allow_file)
    if acceptable:
        logger.info(f"Acceptable off-targets loaded for {len(acceptable)} gene(s).")
    output.mkdir(parents=True, exist_ok=True)

    if gene is not None:
        if gene not in codebook:
            raise ValueError(f"{gene!r} is not in the codebook.")
        genes = [gene]
    else:
        genes = sorted(codebook)

    # Single-gene mode and allow-listed genes force a re-run (parity with the
    # original drivers: accepted off-targets must take effect).
    todo = {
        g: (overwrite or gene is not None or g in acceptable)
        for g in genes
        if overwrite
        or gene is not None
        or g in acceptable
        or not final_parquet(output, g, codebook[g], restriction).exists()
    }
    skipped = [g for g in genes if g not in todo]
    if skipped:
        logger.info(f"Skipping {len(skipped)} finished gene(s); pass --overwrite to redo.")
        # Finished genes are skipped by design, but a finished gene designed
        # under other settings than the ones asked for now is worth naming:
        # otherwise editing the manifest's design block would appear to do
        # nothing for most of the panel.
        stale = [
            g
            for g in skipped
            if not matches_recorded(design, designed_under(final_parquet(output, g, codebook[g], restriction)))
        ]
        if stale:
            logger.warning(
                f"{len(stale)} finished gene(s) were designed under other design settings than "
                f"requested now ({design.describe()}): {', '.join(stale[:10])}"
                f"{' ...' if len(stale) > 10 else ''}. They are kept as they are; pass --overwrite "
                "to redesign them."
            )
    if not todo:
        return {"done": [], "skipped": skipped, "failed": []}

    failed_path = codebook_path.parent / (codebook_path.stem + ".failed.txt")
    failed_path.unlink(missing_ok=True)

    failed: list[str] = []
    # Genes the pool never ran because a worker process died outright (a
    # crash in C code, a memory kill). They did not fail; they were not
    # attempted, and a re-run picks them up. Kept apart from `failed` so the
    # failure file names real failures only.
    unattempted: list[str] = []
    with (
        # "spawn" rather than the original scripts' "forkserver": forkserver
        # deadlocks on macOS (and is Linux-only in practice); spawn is portable
        # and worker startup cost is negligible next to per-gene compute.
        ProcessPoolExecutor(min(workers, len(todo)), mp_context=get_context("spawn")) as exc,
        Progress() as progress,
    ):
        task = progress.add_task("Designing probes", total=len(todo))
        futs = {
            exc.submit(
                run_gene,
                dataset_path,
                output=output,
                codebook=codebook,
                gene=g,
                acceptable=acceptable.get(g),
                overwrite=force,
                restriction=restriction,
                target_probes=target_probes,
                design=design,
                codebook_hash=codebook_hash,
            ): g
            for g, force in todo.items()
        }
        for fut in as_completed(futs):
            g = futs[fut]
            progress.advance(task)
            try:
                fut.result()
            except BrokenProcessPool:
                unattempted.append(g)
            except (Exception, SystemExit) as e:
                # One target must never take the panel down with it.
                failed.append(g)
                logger.critical(f"{g} failed: {e}")
                traceback.print_exception(type(e), e, e.__traceback__, file=sys.stderr)
                with failed_path.open("a") as fh:
                    fh.write(g + "\n")

    done = [g for g in todo if g not in failed and g not in unattempted]
    logger.info(
        f"Panel run complete: {len(done)} done, {len(skipped)} skipped, {len(failed)} failed"
        + (f", {len(unattempted)} not attempted." if unattempted else ".")
    )
    if failed:
        logger.critical(f"Failed genes written to {failed_path}: {failed}")
    if unattempted:
        logger.critical(
            f"A worker process died, so {len(unattempted)} gene(s) were never attempted "
            f"(not written to {failed_path.name}; they have no log and no error of their own). "
            "Look for a crash report in ~/Library/Logs/DiagnosticReports (macOS) or `dmesg` "
            "(Linux), then re-run this command: finished genes are skipped and these are picked up."
        )
    return {"done": done, "skipped": skipped, "failed": failed, "unattempted": unattempted}


@click.command("run-panel")
@click.argument("path", type=click.Path(exists=True, file_okay=False, path_type=Path))
@click.argument("codebook", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.argument("gene", type=str, default=None, required=False)
@click.option("--output", "-o", type=click.Path(file_okay=False, path_type=Path), default=None,
              help="Output directory (default: <codebook dir>/output).")
@click.option("--workers", "-j", type=int, default=DEFAULT_WORKERS, show_default=True,
              help="Parallel worker processes.")
@click.option("--allow-file", type=click.Path(exists=True, dir_okay=False, path_type=Path), default=None,
              help="Per-gene acceptable off-targets JSON "
              "(default: <codebook>.acceptable.json when present).")
@click.option("--manifest", type=click.Path(exists=True, dir_okay=False, path_type=Path), default=None,
              help="Manifest holding the panel's design settings (default: manifest.json beside "
              "the codebook, when it lists it).")
@click.option("--minimum", type=int, default=None,
              help=f"Probes per gene the screen aims for [default: {DEFAULT_MIN_PROBES}].")
@click.option("--maxoverlap", type=int, default=None,
              help="Maximum probe overlap (nt, multiples of 5) tried to reach --minimum "
              f"[default: {DEFAULT_MAX_OVERLAP}]. The overlap that reaches it is the one constructed.")
@design_options
@click.option("--restriction", type=str, default=",".join(DEFAULT_RESTRICTION), show_default=True,
              help="Restriction enzymes, comma-separated.")
@click.option("--target-probes", type=int, default=DEFAULT_TARGET_PROBES, show_default=True,
              help="Maximum probes per gene at the construct stage.")
@click.option("--overwrite", is_flag=True, help="Redo genes whose outputs already exist.")
@click.option("--list-failed", is_flag=True, help="List genes without a final output, then exit.")
@click.option("--list-failed-all", is_flag=True,
              help="Like --list-failed, plus each gene's top off-target counts.")
def run_panel_cli(
    path: Path,
    codebook: Path,
    gene: str | None,
    output: Path | None,
    workers: int,
    allow_file: Path | None,
    manifest: Path | None,
    minimum: int | None,
    maxoverlap: int | None,
    tm_range: tuple[float, float] | None,
    length_range: tuple[int, int] | None,
    split_tm: float | None,
    restriction: str,
    target_probes: int,
    overwrite: bool,
    list_failed: bool,
    list_failed_all: bool,
):
    """Design probes for every target in CODEBOOK (candidates -> screen -> construct, in parallel).

    Give an optional GENE to re-run just that target (forces overwrite for it).

    Design settings (Tm window, length window, split-arm Tm, probes per gene,
    overlap) come from the manifest's `design` block, and a flag given here
    overrides the manifest for this run. Anything set nowhere is the default
    every panel has been designed under.
    """
    from .constants import validate_restriction
    from .ext.ingest import DESIGN_TOOLS, check_external_tools

    output = output or codebook.parent / "output"
    try:
        from_manifest, manifest_used = design_for_codebook(codebook, manifest)
    except ValueError as e:
        raise click.ClickException(str(e)) from e
    from_flags = design_from_flags(tm_range, length_range, split_tm, minimum, maxoverlap)
    design = from_manifest.merged(from_flags)
    if manifest_used is not None:
        logger.info(f"Design settings from {manifest_used}: {from_manifest.describe()}.")
    if not from_flags.is_default():
        logger.info(f"Design settings from the command line: {from_flags.describe()}.")
    enzymes = tuple(e.strip() for e in restriction.split(",") if e.strip())
    try:
        validate_restriction(enzymes)
    except ValueError as e:
        # Fails here rather than after every gene has been designed: assembly
        # only looks for the default pair's filenames.
        raise click.BadParameter(str(e), param_hint="--restriction") from e

    if list_failed or list_failed_all:
        cb = load_worklist(codebook)
        for g in find_missing_final(cb, output, enzymes):
            click.echo(g)
            if list_failed_all:
                counts_path = output / f"{g}_offtarget_counts.csv"
                if counts_path.exists():
                    click.echo(pl.read_csv(counts_path)[:5])
        return

    # A missing aligner would otherwise surface as a raw subprocess error inside
    # every worker process, minutes into the run.
    check_external_tools(DESIGN_TOOLS)

    summary = run_panel(
        path,
        codebook,
        output,
        gene=gene,
        allow_file=allow_file,
        workers=workers,
        overwrite=overwrite,
        restriction=enzymes,
        target_probes=target_probes,
        design=design,
    )
    if summary["failed"] or summary.get("unattempted"):
        raise SystemExit(1)
