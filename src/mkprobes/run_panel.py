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
from multiprocessing import get_context
from pathlib import Path

import polars as pl
import rich_click as click
from loguru import logger
from rich.progress import Progress

# Production defaults, inherited from the original batch drivers.
DEFAULT_RESTRICTION = ("BamHI", "KpnI")
DEFAULT_MINIMUM = 60
DEFAULT_MAXOVERLAP = 0
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


def run_gene(
    dataset_path: Path,
    output: Path,
    codebook: dict[str, list[int]],
    gene: str,
    acceptable: list[str] | None,
    overwrite: bool = False,
    log_level: str = "DEBUG",
    minimum: int = DEFAULT_MINIMUM,
    maxoverlap: int = DEFAULT_MAXOVERLAP,
    restriction: tuple[str, ...] = DEFAULT_RESTRICTION,
    target_probes: int = DEFAULT_TARGET_PROBES,
    **kwargs,
):
    """
    Runs candidates -> screen -> construct for one gene.

    Runs in a worker process: logs to `output/<gene>.log`, skips finished
    genes, and reuses an existing `<gene>_crawled.parquet` unless
    overwriting. An allow-list forces a re-screen/re-construct so accepted
    off-targets take effect. Exceptions are re-raised tagged with the gene.
    """
    # Deferred imports keep worker startup (forkserver) lean.
    from .candidates import get_candidates
    from .codebook.finalconstruct import construct
    from .ext.dataset import load_dataset
    from .screen import run_screen

    logger.remove()
    logger.add(sys.stderr, level=log_level)
    logger.add(output / f"{gene}.log", level=log_level, colorize=False, backtrace=True, diagnose=True)

    if final_parquet(output, gene, codebook[gene], restriction).exists() and not overwrite:
        return

    ds = load_dataset(dataset_path)
    try:
        if overwrite or not (output / f"{gene}_crawled.parquet").exists():
            get_candidates(
                ds,
                transcript=gene,
                output=output,
                ignore_revcomp=False,
                allow=acceptable,
                overwrite=overwrite,
                **kwargs,
            )
            time.sleep(1)  # let parquet writes settle before the next stage reads them
        overwrite = overwrite or acceptable is not None
        run_screen(
            output,
            gene,
            minimum=minimum,
            restriction=list(restriction),
            maxoverlap=maxoverlap,
            overwrite=overwrite,
        )
        construct(
            ds,
            output,
            transcript=gene,
            codebook=codebook,
            restriction=list(restriction),
            target_probes=target_probes,
            overwrite=overwrite,
        )
    except Exception as e:
        # Keep the cause in the message: the panel driver reports this line to a
        # user who cannot see the worker's traceback, and "No probes left after
        # filtering" is the actionable part, not the gene name.
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
    minimum: int = DEFAULT_MINIMUM,
    maxoverlap: int = DEFAULT_MAXOVERLAP,
    restriction: tuple[str, ...] = DEFAULT_RESTRICTION,
    target_probes: int = DEFAULT_TARGET_PROBES,
) -> dict[str, list[str]]:
    """
    Designs probes for every gene in the codebook, in parallel.

    Returns {"done": [...], "skipped": [...], "failed": [...]}. Failed genes
    are also appended to `<codebook>.failed.txt` (recreated per run).
    """
    codebook = load_worklist(codebook_path)
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
    if not todo:
        return {"done": [], "skipped": skipped, "failed": []}

    failed_path = codebook_path.parent / (codebook_path.stem + ".failed.txt")
    failed_path.unlink(missing_ok=True)

    failed: list[str] = []
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
                minimum=minimum,
                maxoverlap=maxoverlap,
                restriction=restriction,
                target_probes=target_probes,
            ): g
            for g, force in todo.items()
        }
        for fut in as_completed(futs):
            g = futs[fut]
            progress.advance(task)
            try:
                fut.result()
            except Exception as e:
                failed.append(g)
                logger.critical(f"{g} failed: {e}")
                traceback.print_exception(type(e), e, e.__traceback__, file=sys.stderr)
                with failed_path.open("a") as fh:
                    fh.write(g + "\n")

    done = [g for g in todo if g not in failed]
    logger.info(f"Panel run complete: {len(done)} done, {len(skipped)} skipped, {len(failed)} failed.")
    if failed:
        logger.critical(f"Failed genes written to {failed_path}: {failed}")
    return {"done": done, "skipped": skipped, "failed": failed}


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
@click.option("--minimum", type=int, default=DEFAULT_MINIMUM, show_default=True,
              help="Minimum probes per gene at the screen stage.")
@click.option("--maxoverlap", type=int, default=DEFAULT_MAXOVERLAP, show_default=True,
              help="Maximum probe overlap tried to reach --minimum.")
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
    minimum: int,
    maxoverlap: int,
    restriction: str,
    target_probes: int,
    overwrite: bool,
    list_failed: bool,
    list_failed_all: bool,
):
    """Design probes for every target in CODEBOOK (candidates -> screen -> construct, in parallel).

    Give an optional GENE to re-run just that target (forces overwrite for it).
    """
    from .constants import validate_restriction
    from .ext.ingest import DESIGN_TOOLS, check_external_tools

    output = output or codebook.parent / "output"
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
        minimum=minimum,
        maxoverlap=maxoverlap,
        restriction=enzymes,
        target_probes=target_probes,
    )
    if summary["failed"]:
        raise SystemExit(1)
