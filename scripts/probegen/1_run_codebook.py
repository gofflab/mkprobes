# %%
"""
DEPRECATED shim: the panel batch driver now lives in the package as
`mkprobes.run_panel` (CLI: `mkprobes run-panel`). This module re-exports the
old surface so existing invocations keep working.
"""

from pathlib import Path

import click
from loguru import logger

from mkprobes.run_panel import run_gene, run_panel  # noqa: F401  (re-exports)


@click.group()
def cli(): ...


@cli.command()
@click.argument("path", type=click.Path(exists=True, dir_okay=True, file_okay=False, path_type=Path))
@click.argument("codebook_path", metavar="CODEBOOK", type=click.Path(exists=True, file_okay=True, path_type=Path))
@click.argument("gene", type=str, default=None, required=False)
@click.option("--overwrite", is_flag=True)
@click.option("--listfailed", is_flag=True)
@click.option("--listfailedall", is_flag=True)
def single(
    path: Path,
    codebook_path: Path,
    gene: str | None,
    overwrite: bool = False,
    listfailed: bool = False,
    listfailedall: bool = False,
):
    logger.warning("This script is deprecated; use `mkprobes run-panel`.")
    if listfailed or listfailedall:
        from mkprobes.run_panel import DEFAULT_RESTRICTION, find_missing_final, load_worklist

        for g in find_missing_final(
            load_worklist(codebook_path), codebook_path.parent / "output", DEFAULT_RESTRICTION
        ):
            print(g)
        return
    run_panel(
        path,
        codebook_path,
        codebook_path.parent / "output",
        gene=gene,
        overwrite=overwrite,
    )


if __name__ == "__main__":
    cli()
