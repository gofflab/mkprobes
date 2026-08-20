# Phase 4: Panel QC and export

Phase 4 performs panel-level quality control after per-target construction of SOLAR (splint/padlock, STARmap-style) probes.

## Why this phase matters

Per-target success does not guarantee panel success. This phase enforces panel-level acceptance criteria so low-yield targets are identified and fixed before downstream ordering/analysis.

## Inputs

- Completed `output/*_final_*.parquet` files from phase 3.
- Original target list file.

## 4A. Probe-count QC by target

Filter targets by minimum probe count:

```bash
uv run mkprobes filter-genes panel_a/output --genes panel_a/genes.converted.txt --min-probes 48 --out panel_a/genes.pass.txt
```

Use this to identify targets requiring reruns/parameter tuning.

Key arguments explained:

- `OUTPUT_PATH`:
  - directory containing `*_screened_ol*.parquet` outputs used for counting. For each gene, the highest-overlap file `{gene}_screened_ol{N}[_{enzymes}].parquet` is read (`N` may be negative; the default is `-2`).
- `--genes`:
  - authoritative target list for QC comparison.
- `--min-probes`:
  - threshold defining pass/fail at panel level. Genes with exactly `--min-probes` probes pass (the comparison is `>=`).
- `--out`:
  - file receiving passing targets for reproducible downstream handoff.

## 4B. Iterate low-yield targets

Recommended loop:

1. Review low-count genes from QC output.
2. Re-run selected genes in phase 3 with adjusted options (`--allow`, `--minimum`, overlap strategy, etc.).
3. Re-run `filter-genes` until panel thresholds are met.

**Intention:** convert QC into a controlled rerun cycle rather than manual ad hoc fixes.

## 4C. Prepare downstream exports

The primary output files for downstream ordering/processing are the per-target final parquet files:

- `output/<target>_final_<restriction>_<bits>.parquet`

Final assembly into an orderable oligo pool is phase 5: {doc}`phase_5_manifest_assembly`.

**Intention:** keep the CLI workflow responsible for validated design outputs, while allowing export formatting to remain modular.

## HPC add-on guidance

- Run QC as a lightweight post-job step.
- Archive: codebook JSON, target list, QC outputs, and final parquet files as a single panel release bundle.
- Keep rerun jobs target-scoped to minimize wasted compute.

## Failure modes

- Missing `*_final_*.parquet` for some targets:
  - rerun phase 3 for missing targets.
- Probe counts remain below threshold after retries:
  - revisit target transcript choice and off-target acceptance policy.
