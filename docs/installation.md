# Installation

`mkprobes` lives in its own repository ([github.com/gofflab/mkprobes](https://github.com/gofflab/mkprobes)); the in-tree copy inside the fishtools repository is frozen legacy. Two supported setups: `uv` (Option A) or a conda/mamba environment (Option B). Both need the same external tools; Option B installs them for you from bioconda.

## Requirements

- OS: Linux/macOS recommended.
- Python: `>=3.12` (from `pyproject.toml`).
- External tools (needed on `PATH` for CLI probe design/prep; all available on bioconda):
  - `bowtie2`
  - `jellyfish` (bioconda package `kmer-jellyfish`)
  - `gffread`
  - `RepeatMasker` (optional; only used at final manifest assembly)

All commands in this documentation are written bare (`mkprobes ...`), assuming your environment is **activated**. Activate once per shell session as shown below.

## Option A: install with uv

From the repository root:

```bash
uv sync
```

This creates/updates the project virtual environment and installs `mkprobes` with all Python dependencies. External tools must already be on `PATH` (e.g. via Homebrew, module system, or a separate conda env).

Activate the environment, then run commands directly:

```bash
source .venv/bin/activate
mkprobes --help
```

## Option B: install with conda/mamba

The repository ships an [`environment.yml`](https://github.com/gofflab/mkprobes/blob/main/environment.yml) that builds one environment containing the Python package **and** the external tools (bowtie2, jellyfish, gffread) from bioconda. From the repository root:

```bash
mamba env create -f environment.yml
```

(`conda env create -f environment.yml` works identically, just slower to solve.) Then:

```bash
conda activate mkprobes
mkprobes --help
```

Notes:

- Create the environment **from the repository root** — the package is installed editable (`pip: -e .`), so the path is relative.
- RepeatMasker is commented out in `environment.yml` (it is optional and pulls large repeat libraries); enable it there or add it later:

  ```bash
  mamba install -n mkprobes -c bioconda repeatmasker
  ```

- After a `git pull` the editable install picks up code changes automatically; re-run `mamba env update -f environment.yml` only when dependencies change.

## Verify installation

With your environment activated:

```bash
mkprobes --help
python -c "import mkprobes; print('mkprobes import OK')"
gffread --version && bowtie2 --version | head -1 && jellyfish --version
```

To run the test suite:

```bash
pip install -e ".[dev]"
pytest
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
  - Activate it (`source .venv/bin/activate` or `conda activate mkprobes`) and retry; re-run `uv sync` / `mamba env update -f environment.yml` if it is stale.
- `mkprobes prepare` or `mkprobes ingest` fails on external binaries:
  - Confirm `bowtie2`, `jellyfish`, and `gffread` are available in `PATH` on both login and compute nodes.
