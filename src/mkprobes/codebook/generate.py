"""
Codebook generation for SOLAR panels.

Builds an MHD (minimum-Hamming-distance) codebook mapping each target to a
set of readout bits, from the code matrices vendored inside the package.
Codeword assignment is seeded-random by default and **optionally
expression-informed**: given per-target expression values, assignment is
optimized (via :class:`~mkprobes.codebook.codebook.CodebookPicker.find_optimalish`)
to balance total expression load across readout bits, so no bit's
fluorescence is dominated by a few highly expressed genes. Expression data
is never required — without it, generation falls back to a plain seeded
shuffle.

Exposed on the CLI as ``mkprobes make-codebook``.
"""

import json
import os
import re
from collections.abc import Mapping, Sequence
from importlib.resources import files
from itertools import chain
from pathlib import Path

import numpy as np
import numpy.typing as npt
import polars as pl
import rich_click as click
from loguru import logger

from ..ext.dataset import Dataset, _read_annotation_table, load_dataset
from .codebook import CodebookPicker, bit_count, hash_codebook, n_to_bit

# MHD code matrices are vendored inside the package; any matrix not vendored
# is generated into a user cache dir (never into the installed package).
VENDORED_MHD = Path(str(files("mkprobes") / "data" / "mhd"))
MHD_CACHE = Path(os.environ.get("MKPROBES_MHD_CACHE", str(Path.home() / ".cache" / "mkprobes" / "mhd")))

# Maps 0-based codeword bit positions to readout IDs, interleaved across
# imaging rounds (bits 1-24 cycle rounds in strides of 8; 25-49 sequential).
ORDER: tuple[int, ...] = tuple(
    chain.from_iterable([[i, i + 8, i + 16] for i in range(1, 9)])
) + tuple(range(25, 50))

# Readout-ID triplets that perfectly confound imaging rounds; genes assigned
# one of these are swapped with a blank codeword.
FORBIDDEN: frozenset[tuple[int, ...]] = frozenset(
    {
        *[(x, x + 8, x + 16) for x in range(1, 9)],
        *[(3 * i, 3 * i + 1, 3 * i + 2) for i in range(9, 13)],
    }
)

# Auto-selection considers codes of at least this many bits (8/9-bit matrices
# are special-purpose and only used when requested explicitly via n_bits).
MIN_AUTO_BITS = 10


def _gen_mhd(n: int, on: int, min_dist: int = 4, seed: int = 0) -> npt.NDArray[np.uint32]:
    """Generates an n-bit code with `on` ones per word and pairwise Hamming distance >= min_dist."""
    assert n < 32
    rand = np.random.default_rng(seed)

    while (s := rand.integers(0, 2**n - 1)).bit_count() != on:
        ...

    out = np.zeros((2 ** (n - 1)), dtype=np.uint32)
    out[0] = s
    cnt = 1
    for i in range(2**n):
        if i.bit_count() != on:
            continue
        if np.any(bit_count(out[:cnt] ^ i) < min_dist):
            continue
        out[cnt] = i
        cnt += 1
    return out[:cnt]


def _generate_matrix(path: Path, n: int) -> np.ndarray:
    x = _gen_mhd(n, 3, seed=0, min_dist=2)
    np.savetxt(path, out := n_to_bit(x, n, 3), fmt="%d", delimiter=",")
    return out


def discover_matrices() -> dict[int, Path]:
    """
    Discovers MHD code matrices: vendored ones win, the user cache holds any
    generated extras. Missing sizes in 10..30 are generated into the cache.
    """
    matrix_paths: dict[int, Path] = {}
    for source in (VENDORED_MHD, MHD_CACHE):
        if source.exists():
            for path in sorted(source.glob("*bit_on3_dist2.csv")):
                matrix_paths.setdefault(int(re.search(r"(\d+)", path.stem).group(1)), path)
    for n in range(MIN_AUTO_BITS, 31):
        if n not in matrix_paths:
            MHD_CACHE.mkdir(parents=True, exist_ok=True)
            matrix_paths[n] = MHD_CACHE / f"{n}bit_on3_dist2.csv"
            logger.info(f"Generating {n}-bit MHD matrix into {MHD_CACHE}.")
            _generate_matrix(matrix_paths[n], n)
    return matrix_paths


def _capacities(matrix_paths: Mapping[int, Path]) -> dict[int, int]:
    return {n: len(path.read_text().splitlines()) for n, path in matrix_paths.items()}


def choose_bits(n_genes: int, n_bits: int | None, capacities: Mapping[int, int]) -> int:
    """
    Chooses the code size: explicit `n_bits`, or the smallest code (>= 10
    bits) whose capacity exceeds the gene count by at least 5%.
    """
    if n_bits is not None:
        n = int(n_bits)
        if n not in capacities:
            raise ValueError(f"No {n}-bit code matrix available (have: {sorted(capacities)}).")
        if (capacities[n] - n_genes) / capacities[n] < 0.05:
            logger.warning("Less than 5% of coding capacity is blank")
        return n
    for n in sorted(capacities):
        if n < MIN_AUTO_BITS:
            continue
        if n_genes * 1.05 < capacities[n]:
            if (capacities[n] - n_genes) / capacities[n] < 0.05:
                logger.warning("Less than 5% of coding capacity is blank")
            return n
    raise ValueError(f"No suitable codebook found. {n_genes} genes found.")


def make_codebook(
    genes: Sequence[str],
    *,
    expression: Mapping[str, float] | None = None,
    n_bits: int | None = None,
    offset: int = 0,
    existing_codebook: Mapping[str, Sequence[int]] | None = None,
    seed: int = 0,
    iterations: int = 200,
) -> dict[str, list[int]]:
    """
    Generates a SOLAR codebook: {target: [bit, bit, bit]}, plus Blank-N
    decoy codewords filling the remaining capacity.

    Codeword assignment is a seeded shuffle by default. When `expression`
    (target -> expression value, e.g. FPKM/TPM) is given, the assignment is
    chosen from `iterations` seeded shuffles to maximize the entropy of
    per-bit expression load — balancing fluorescence across readout bits.
    Expression is strictly optional; omitting it changes nothing else.

    `existing_codebook` extends a previous panel: its bit range sets the
    offset, and gene/bit overlap is rejected.

    Readout-round confounder codewords (FORBIDDEN) are swapped onto blanks,
    matching the historical generator.
    """
    genes = list(genes)
    if existing_codebook is not None:
        if offset:
            raise ValueError("Specify either offset or existing_codebook, not both.")
        existing_bits = set(chain.from_iterable(existing_codebook.values()))
        offset = len(existing_bits)
        logger.info(f"Using offset {offset} from existing codebook")
        if overlap := set(existing_codebook) & set(genes):
            raise ValueError(f"Genes in existing codebook and input overlap: {sorted(overlap)}")

    matrix_paths = discover_matrices()
    capacities = _capacities(matrix_paths)
    n = choose_bits(len(genes), n_bits, capacities)
    logger.info(f"Using {n}-bit codebook with capacity {capacities[n]}.")
    if n + offset > len(ORDER):
        raise ValueError(
            f"offset {offset} too large for a {n}-bit code: only {len(ORDER)} readout IDs available."
        )

    cb = CodebookPicker(matrix_paths[n], genes=genes)

    if expression is not None:
        missing = [g for g in genes if g not in expression]
        if missing:
            raise ValueError(f"Expression values missing for {len(missing)} targets: {missing[:10]}")
        vec = np.asarray([float(expression[g]) for g in genes])
        chosen_seed, loads = cb.find_optimalish(vec, iterations=iterations)
        logger.info(
            f"Expression-informed assignment: per-bit load range "
            f"{loads.min():.1f}-{loads.max():.1f} (mean {loads.mean():.1f})."
        )
    else:
        chosen_seed = seed

    c = cb.export_codebook(chosen_seed, offset=0)
    out = {k: sorted(ORDER[x + offset] for x in v) for k, v in c.items()}

    to_swap = [k for k, v in out.items() if tuple(v) in FORBIDDEN and not k.startswith("Blank")]
    for i, k in enumerate(to_swap):
        out[f"Blank-{i + 1}"], out[k] = out[k], out[f"Blank-{i + 1}"]

    if existing_codebook is not None:
        new_bits = set(chain.from_iterable(out.values()))
        if new_bits & existing_bits:
            raise ValueError(f"Bits overlap with existing codebook: {sorted(new_bits & existing_bits)}")

    logger.info("Bits used: " + str(sorted(set(chain.from_iterable(out.values())))))
    return {k: out[k] for k in sorted(out, key=lambda x: (x.startswith("Blank"), x))}


def resolve_expression(
    dataset: Dataset | None,
    spec: str,
    genes: Sequence[str],
    column: str | None = None,
) -> dict[str, float]:
    """
    Resolves an expression source into {target: value} aligned to `genes`.

    `spec` is either the name of an annotation table registered in the
    dataset (preferred) or a path to a parquet/csv/tsv file. The table must
    carry a `transcript_id` and/or `gene_id` column; targets are matched
    against both. The value column is `column` when given, otherwise the
    table's single numeric column (ambiguity is an error).

    Targets absent from the table are filled with the table's median value
    (neutral for load balancing) with a warning.
    """
    if dataset is not None and spec in dataset.annotation_paths:
        table = dataset.annotation(spec)
        source = f"annotation table {spec!r}"
    elif Path(spec).exists():
        table = _read_annotation_table(Path(spec))
        source = str(spec)
    else:
        registered = sorted(dataset.annotation_paths) if dataset is not None else []
        raise ValueError(
            f"Expression source {spec!r} is neither a registered annotation table "
            f"(available: {registered or 'none'}) nor an existing file."
        )

    id_cols = [c for c in ("transcript_id", "gene_id") if c in table.columns]
    if not id_cols:
        raise ValueError(f"Expression table {source} has no transcript_id/gene_id column.")

    if column is None:
        numeric = [c for c in table.columns if c not in id_cols and table[c].dtype.is_numeric()]
        if len(numeric) != 1:
            raise ValueError(
                f"Expression table {source} has {len(numeric)} numeric columns ({numeric}); "
                "pass --expression-column to choose one."
            )
        column = numeric[0]
    elif column not in table.columns:
        raise ValueError(f"Column {column!r} not in expression table {source} ({table.columns}).")

    values: dict[str, float] = {}
    for id_col in id_cols:
        for key, value in table.select(id_col, column).drop_nulls().iter_rows():
            values.setdefault(key, float(value))

    matched = {g: values[g] for g in genes if g in values}
    missing = [g for g in genes if g not in values]
    if missing:
        fill = float(np.median(list(values.values()))) if values else 0.0
        logger.warning(
            f"{len(missing)}/{len(genes)} targets absent from {source} "
            f"(e.g. {missing[:5]}); filling with the table median ({fill:.2f})."
        )
        matched |= {g: fill for g in missing}
    logger.info(f"Expression from {source}, column {column!r}: {len(genes) - len(missing)}/{len(genes)} targets matched.")
    return matched


@click.command("make-codebook")
@click.argument("path", type=click.Path(exists=True, file_okay=False, path_type=Path))
@click.argument("genes", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option("--out", "-o", type=click.Path(dir_okay=False, path_type=Path), default=None,
              help="Output JSON path (default: <genes>.codebook.json).")
@click.option("--expression", "-e", "expression_spec", type=str, default=None,
              help="OPTIONAL: expression source for load-balanced assignment - the name of an "
              "annotation table registered in the dataset, or a parquet/csv/tsv file path. "
              "Omit for plain seeded assignment.")
@click.option("--expression-column", type=str, default=None,
              help="Value column in the expression table (needed only when ambiguous).")
@click.option("--n-bits", type=int, default=None, help="Code size; auto-sized from gene count if omitted.")
@click.option("--offset", type=int, default=0, help="Readout-ID offset (mutually exclusive with --existing-codebook).")
@click.option("--existing-codebook", type=click.Path(exists=True, dir_okay=False, path_type=Path),
              default=None, help="Extend this codebook: derives the offset and rejects gene/bit overlap.")
@click.option("--iterations", type=int, default=200, show_default=True,
              help="Assignments tried when --expression is given.")
@click.option("--seed", type=int, default=0, show_default=True, help="Assignment seed (uninformed mode).")
def make_codebook_cli(
    path: Path,
    genes: Path,
    out: Path | None,
    expression_spec: str | None,
    expression_column: str | None,
    n_bits: int | None,
    offset: int,
    existing_codebook: Path | None,
    iterations: int,
    seed: int,
):
    """Generate a codebook for a target list, optionally expression-informed."""
    targets = [g for g in genes.read_text().split() if g]
    if len(targets) != len(set(targets)):
        raise click.ClickException("Duplicate targets in the gene list.")

    expression = None
    if expression_spec is not None:
        dataset = load_dataset(path)
        expression = resolve_expression(dataset, expression_spec, targets, column=expression_column)

    codebook = make_codebook(
        targets,
        expression=expression,
        n_bits=n_bits,
        offset=offset,
        existing_codebook=json.loads(existing_codebook.read_text()) if existing_codebook else None,
        seed=seed,
        iterations=iterations,
    )

    out = out or genes.with_suffix(".codebook.json")
    out.write_text(json.dumps(codebook, indent=2))
    n_blanks = sum(k.startswith("Blank") for k in codebook)
    logger.info(
        f"Codebook written to {out}: {len(codebook) - n_blanks} targets + {n_blanks} blanks, "
        f"hash {hash_codebook(codebook)}"
        + (" (expression-informed)" if expression is not None else " (seeded assignment)")
    )
