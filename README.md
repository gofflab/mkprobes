# mkprobes

[![test](https://github.com/gofflab/mkprobes/actions/workflows/test.yml/badge.svg)](https://github.com/gofflab/mkprobes/actions/workflows/test.yml)

`mkprobes` is the combinatorial FISH probe-design toolkit used to build **SOLAR**
probesets (splint/padlock, STARmap-style). It provides a CLI covering the full
workflow, from a reference to an orderable oligo pool:

```text
0. project    mkprobes init                  scaffold a project
1. dataset    mkprobes prepare | ingest      reference (mouse/human) | any species
2. targets    mkprobes chkgenes / convert-to-transcripts
3. codebook   mkprobes make-codebook
4. probes     mkprobes run-panel             candidates -> screen -> construct, all targets
5. panel QC   mkprobes filter-genes
6. assembly   mkprobes assemble              -> orderable oligos
```

Start with [Getting started](https://www.gofflab.org/mkprobes/getting_started.html),
which walks through all of it. `mkprobes --help` prints the same order, and
every command takes `--help`.

## Provenance

This repository was extracted from
[chaichontat/fishtools](https://github.com/chaichontat/fishtools) (commit
`cd91ef7`, "Add standalone mkprobes package"; extracted 2026-08-20). The full
pre-extraction development history remains browsable there under
`fishtools/mkprobes/`, which is now frozen legacy — all probe-design
development happens here. The MHD codebook matrices (`src/mkprobes/data/mhd/`)
and the splint/padlock header/footer table
(`src/mkprobes/data/headerfooter.csv`) were vendored from the same repository
(`static/` and `data/`).

## Requirements

- Python ≥ 3.12
- External tools on `PATH` for dataset preparation and screening:
  [Bowtie 2](https://bowtie-bio.sourceforge.net/bowtie2/),
  [Jellyfish](https://github.com/gmarcais/Jellyfish), and (for final panel
  assembly) [RepeatMasker](https://www.repeatmasker.org/). All are available
  via bioconda.

## Installation

With `uv` (external tools must already be on `PATH`):

```bash
uv sync --extra dev
source .venv/bin/activate
mkprobes --help
```

Or with conda/mamba — one environment containing the package *and* the
external tools (bowtie2, jellyfish, gffread) from bioconda:

```bash
mamba env create -f environment.yml
conda activate mkprobes
mkprobes --help
```

See the [installation docs](https://www.gofflab.org/mkprobes/installation.html)
for details.

## Tests

With the environment activated:

```bash
pip install -e ".[dev]"
pytest
```

## Notes

- The experimental `picker` and `readouts` notebooks from the legacy package
  are excluded (import-time side effects, stale dependencies).
- Project-specific panel generation and assembly drivers are retained under
  `scripts/`. They are distributed in the source archive, not installed as
  console commands. `scripts/generate_mhd_matrices.py` is the (scratch)
  generator for the vendored MHD matrices.
