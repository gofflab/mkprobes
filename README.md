# mkprobes

`mkprobes` is the combinatorial FISH probe-design toolkit used to build **SOLAR**
probesets (splint/padlock, STARmap-style). It provides a CLI covering the full
workflow: reference/dataset preparation, per-target candidate generation,
off-target screening, probe construction against a codebook, and panel QC.

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

```bash
uv sync --extra dev
uv run mkprobes --help
```

## Tests

```bash
uv run --extra dev pytest
```

## Notes

- The experimental `picker` and `readouts` notebooks from the legacy package
  are excluded (import-time side effects, stale dependencies).
- Project-specific panel generation and assembly drivers are retained under
  `scripts/`. They are distributed in the source archive, not installed as
  console commands. `scripts/generate_mhd_matrices.py` is the (scratch)
  generator for the vendored MHD matrices.
