# %%
"""
DEPRECATED shim: codebook generation now lives in the package as
`mkprobes.codebook.generate` (CLI: `mkprobes make-codebook`, which also
supports optional expression-informed assignment). This module re-exports
the old surface so existing scripts keep working.
"""

from pathlib import Path

import rich_click as click
from loguru import logger

from mkprobes.codebook.generate import (  # noqa: F401  (re-exports)
    FORBIDDEN,
    MHD_CACHE,
    ORDER as order,
    VENDORED_MHD,
    _gen_mhd,
    _generate_matrix as _generate,
    discover_matrices,
    make_codebook,
)
from mkprobes.codebook.generate import _capacities as __capacities

# Legacy module-level names.
static = VENDORED_MHD
matrix_paths = discover_matrices()
ns = __capacities(matrix_paths)


def gen_codebook(tss: list[str], offset: int = 0, n_bits: int | None = None, seed: int = 0):
    logger.warning(
        "scripts/probegen/o_codebook.py is deprecated; use `mkprobes make-codebook` "
        "(mkprobes.codebook.generate.make_codebook)."
    )
    return make_codebook(tss, n_bits=n_bits, offset=offset, seed=seed)


@click.command()
@click.argument("path", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option("--offset", type=int, default=0)
@click.option("--n-bits", type=int, default=None)
@click.option(
    "--existing-codebook", type=click.Path(exists=True, dir_okay=False, path_type=Path), default=None
)
def run(path: Path, offset: int = 0, existing_codebook: Path | None = None, n_bits: int | None = None):
    import json

    genes = path.read_text().splitlines()
    existing = json.loads(existing_codebook.read_text()) if existing_codebook else None
    if offset and existing is not None:
        raise ValueError("Must specify either offset or existing codebook")
    generated = make_codebook(genes, offset=offset, n_bits=n_bits, existing_codebook=existing)
    path.with_suffix(".json").write_text(json.dumps(generated, indent=2))


# %%
if __name__ == "__main__":
    run()
