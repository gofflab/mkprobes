import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Literal

import rich_click as click
from loguru import logger

from .assembly import cli as assemble
from .candidates import candidates
from .codebook.finalconstruct import click_construct, filter_genes
from .codebook.generate import make_codebook_cli
from .ext.dataset import Dataset, create_dataset
from .ext.ingest import ingest
from .genes.chkgenes import chkgenes, convert_to_transcripts, transcripts
from .init_project import check_manifest_cli, init
from .run_panel import run_panel_cli
from .screen import screen
from .select_targets import suggest_targets
from .utils._alignment import bowtie_build
from .utils.logging import setup_logging

setup_logging()
click.rich_click.SHOW_ARGUMENTS = True
click.rich_click.GROUP_ARGUMENTS_OPTIONS = True
click.rich_click.USE_MARKDOWN = True
click.rich_click.STYLE_HELPTEXT = ""


class FriendlyGroup(click.RichGroup):
    """
    Reports failures as one actionable line instead of a traceback.

    Probe design runs for a long time and its users are biologists, not Python
    developers; a wall of polars frames tells them nothing about what to fix.
    The traceback is still one flag away, and always reaches the log file.
    """

    def invoke(self, ctx: click.Context):
        try:
            return super().invoke(ctx)
        except (click.ClickException, click.Abort, click.exceptions.Exit, SystemExit):
            # click.exceptions.Exit is a RuntimeError, and carries the exit code
            # for `--help` and for any command that exits deliberately.
            raise
        except Exception as exc:
            if ctx.params.get("debug"):
                raise
            logger.opt(exception=exc).debug("Command failed")
            # Name the whole invocation: --debug belongs to the group, so
            # appending it to the subcommand - the obvious reading of a bare
            # "re-run with --debug" - fails with "No such option".
            retry = f"mkprobes --debug {ctx.invoked_subcommand} ..." if ctx.invoked_subcommand else "mkprobes --debug ..."
            raise click.ClickException(
                f"{type(exc).__name__}: {exc}\n\n"
                f"For the full traceback, re-run as `{retry}` "
                "(--debug goes before the command name, not after)."
            ) from exc


@click.group(cls=FriendlyGroup)
@click.option("--debug", is_flag=True, help="Show the full traceback when a command fails.")
def main(debug: bool):
    """Design SOLAR probe sets, from a reference to an orderable oligo pool.

    The workflow is six steps, each one command:

    1. **dataset** - `mkprobes prepare` (mouse/human) or `mkprobes ingest` (any species)
    2. **targets** - `mkprobes chkgenes`, then `mkprobes convert-to-transcripts`
    3. **codebook** - `mkprobes make-codebook`
    4. **probes** - `mkprobes run-panel` (candidates, screen and construct for every target)
    5. **panel QC** - `mkprobes filter-genes`
    6. **assembly** - `mkprobes assemble`

    Run any command with `--help` for its arguments. Full guide:
    https://www.gofflab.org/mkprobes/
    """
    if debug:
        setup_logging("DEBUG")


# fmt: off
@main.command()
@click.argument("path", type=click.Path(dir_okay=True, file_okay=False, path_type=Path))
@click.option("--species", "-s", type=click.Choice(("human", "mouse")), default="mouse", help="Species to use for probe design")
@click.option("--threads", "-t", type=int, default=16, help="Number of threads to use")
# fmt: on
def prepare(path: Path, species: Literal["human", "mouse"], threads: int = 16):
    """Prepare genomic database"""
    from .ext.prepare import download_gtf_fasta, run_jellyfish

    path = path.resolve()
    download_gtf_fasta(path / species, species)
    with ThreadPoolExecutor() as exc:
        futs = [
            exc.submit(run_jellyfish, path / species),
            exc.submit(bowtie_build, path / species / "cdna_ncrna_trna.fasta", "txome"),
        ]
        for fut in as_completed(futs):
            fut.result()
    Dataset(path / species)  # test all components


@main.command()
@click.argument("path", type=click.Path(exists=True, dir_okay=False, file_okay=True, path_type=Path))
@click.option(
    "--write",
    "write_sidecar",
    is_flag=True,
    help="Also write the hash to <codebook>.hash, as make-codebook does. Use this on a "
    "codebook made before mkprobes wrote the sidecar; it does not alter the codebook.",
)
def hash(path: Path, write_sidecar: bool):
    """Print a codebook's hash, the stable identifier for that set of bit assignments."""
    from .codebook.codebook import hash_codebook_file

    digest = hash_codebook_file(path)
    click.echo(digest)

    if not write_sidecar:
        return
    sidecar = path.with_suffix(".hash")
    if sidecar.exists() and (previous := sidecar.read_text().strip()) != digest:
        # The sidecar is meant to identify the codebook beside it, so it is
        # updated - but a mismatch means the codebook was edited after the
        # sidecar was written, and anything designed in between belongs to
        # neither hash.
        logger.warning(
            f"{sidecar.name} held {previous}, but {path.name} now hashes to {digest}. "
            "The codebook changed after the sidecar was written; outputs designed against "
            "the old one carry the old hash."
        )
    sidecar.write_text(digest + "\n")
    logger.info(f"Hash written to {sidecar}.")


@main.command()
@click.argument("path", type=click.Path(exists=True, dir_okay=False, path_type=Path))
def provenance(path: Path):
    """Show how an output file was made: version, parameters, and dataset."""
    from .utils.provenance import read_provenance

    record = read_provenance(path)
    if record is None:
        raise click.ClickException(
            f"{path} carries no provenance. Files written before mkprobes recorded it, "
            "and files not written by mkprobes, have none. Re-run the step that produced "
            "it to stamp a fresh copy."
        )
    click.echo(json.dumps(record, indent=2, sort_keys=True))


main.add_command(candidates)
main.add_command(screen)
main.add_command(chkgenes)
main.add_command(filter_genes)
main.add_command(transcripts)
main.add_command(convert_to_transcripts)
main.add_command(click_construct)
main.add_command(create_dataset)
main.add_command(ingest)
main.add_command(make_codebook_cli)
main.add_command(run_panel_cli)
main.add_command(assemble)
main.add_command(init)
main.add_command(check_manifest_cli)
main.add_command(suggest_targets)


if __name__ == "__main__":
    main()
