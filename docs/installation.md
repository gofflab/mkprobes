# Installation

`mkprobes` lives in its own repository ([github.com/gofflab/mkprobes](https://github.com/gofflab/mkprobes)); the in-tree copy inside the fishtools repository is frozen legacy. Start with a reproducible environment, then install with `uv`.

## Requirements

- OS: Linux/macOS recommended.
- Python: `>=3.12` (from `pyproject.toml`).
- External tools (needed on `PATH` for CLI probe design/prep; all available on bioconda):
  - `bowtie2`
  - `jellyfish`
  - `gffread`
  - `RepeatMasker` (optional; only used at final manifest assembly)

## Install with uv

From the repository root:

```bash
uv sync
```

This creates/updates the project virtual environment and installs `mkprobes` with all dependencies.

## Verify installation

Run these checks from repo root:

```bash
uv run mkprobes --help
uv run python -c "import mkprobes; print('mkprobes import OK')"
```

If `mkprobes` is not found, ensure you are running through `uv run` (or have activated the project's `.venv`).

To run the test suite:

```bash
uv run --extra dev pytest
```

## HPC/cluster add-on guidance

Recommended pattern for cluster runs:

1. Build and test one environment interactively.
2. Freeze/package it for compute nodes.
3. Keep large intermediate outputs on node-local scratch, then copy final outputs to shared storage.

Practical tips:

- Put reference dataset prep (`mkprobes prepare` or `mkprobes ingest`) in a shared, versioned location to avoid repeated downloads/indexing.
- Stage per-job output on local scratch (for high I/O steps) then sync back.
- Cap thread/process settings to match scheduler allocation:
  - `mkprobes prepare --threads <N>`

## Common install failures

- `ModuleNotFoundError: click` or `loguru`:
  - You are likely outside the intended environment.
  - Re-run commands through `uv run`, or re-run `uv sync`.
- `mkprobes prepare` or `mkprobes ingest` fails on external binaries:
  - Confirm `bowtie2`, `jellyfish`, and `gffread` are available in `PATH` on both login and compute nodes.
