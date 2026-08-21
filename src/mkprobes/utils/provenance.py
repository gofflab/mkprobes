"""
Provenance stamped into every output file.

A probe set is a physical order: someone synthesizes these sequences and spends
real money and bench time on them. Months later the question "which version of
which dataset, with which thresholds, produced this file?" has to be answerable
from the file itself, not from whoever happened to run it.

Every parquet the pipeline writes therefore carries a `mkprobes` key in its
parquet key/value metadata holding a JSON record: the package version, a UTC
timestamp, the command line, the dataset it was designed against, and the
parameters that shaped the output. Read it back with :func:`read_provenance`,
or from the shell with ``mkprobes provenance <file>``.
"""

import json
import sys
from datetime import UTC, datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import polars as pl

#: Parquet key/value metadata key holding the JSON provenance record.
PROVENANCE_KEY = "mkprobes"


def mkprobes_version() -> str:
    """Installed package version, or `"unknown"` when running from a source tree."""
    try:
        return version("mkprobes")
    except PackageNotFoundError:  # pragma: no cover - only when not installed
        return "unknown"


def provenance_record(dataset: Path | str | None = None, **parameters: Any) -> dict[str, Any]:
    """
    Builds the provenance record for one output.

    `parameters` should carry whatever actually shaped this file - thresholds,
    enzymes, overlap, probe counts. Values of `None` are dropped so callers can
    pass optional arguments through without special-casing them.
    """
    record: dict[str, Any] = {
        "mkprobes_version": mkprobes_version(),
        "generated_at": datetime.now(UTC).isoformat(timespec="seconds"),
        "command": " ".join(sys.argv),
    }
    if dataset is not None:
        record["dataset"] = str(Path(dataset).resolve())
    record.update({key: value for key, value in parameters.items() if value is not None})
    return record


def encode(record: dict[str, Any]) -> dict[str, str]:
    """A provenance record encoded for ``DataFrame.write_parquet(metadata=...)``."""
    return {PROVENANCE_KEY: json.dumps(record, default=str, sort_keys=True)}


def provenance_metadata(dataset: Path | str | None = None, **parameters: Any) -> dict[str, str]:
    """Builds and encodes a record in one step, for writes that need nothing else."""
    return encode(provenance_record(dataset, **parameters))


def read_provenance(path: Path | str) -> dict[str, Any] | None:
    """
    Provenance recorded in a parquet file, or `None` for files written before
    provenance was added (or by another tool).
    """
    try:
        metadata = pl.read_parquet_metadata(path)
    except (OSError, pl.exceptions.PolarsError):
        # Missing, unreadable, or not a parquet file at all.
        return None
    raw = metadata.get(PROVENANCE_KEY)
    if raw is None:
        return None
    try:
        return json.loads(raw)
    except json.JSONDecodeError:
        return None
