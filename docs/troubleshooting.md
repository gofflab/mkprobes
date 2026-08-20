# Troubleshooting

## Installation and environment

- If the `mkprobes` command is missing, activate your environment first (`source .venv/bin/activate` for uv installs, `conda activate mkprobes` for conda installs); re-run `uv sync` or `mamba env update -f environment.yml` if the environment is stale.
- If external tools fail, confirm `bowtie2`, `jellyfish`, and `gffread` are on `PATH` (all available on bioconda).

## Probe design

- `Transcript not found in ensembl`: confirm whether input is gene name vs transcript and run `mkprobes transcripts`. On custom datasets, use the offline selection modes (`--longest`/`--all`).
- Very low final probe counts: rerun with adjusted acceptable off-targets and inspect generated `*.stats.json` files.
- Existing outputs skipped unexpectedly: rerun with `--overwrite`.
