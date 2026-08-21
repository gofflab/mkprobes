# Probe-design dataflow map

This page gives a compact, execution-oriented map of inputs and outputs through the `mkprobes` CLI workflow for SOLAR (splint/padlock, STARmap-style) probesets.

## Phase-by-phase input/output map

| Phase | Primary command(s) | Required inputs | Primary outputs | Purpose |
|---|---|---|---|---|
| 1. Dataset prep | `mkprobes prepare`, `mkprobes ingest`, or `mkprobes create-dataset` | species/ref selection, genome + GTF/GFF3, or FASTA | indexed dataset folder (`txome`, `cdna18.jf`, etc.) | define sequence universe + search indices |
| 2. Targets + codebook | `mkprobes chkgenes`, `mkprobes convert-to-transcripts`, `mkprobes make-codebook` | target list (+ optional expression table) | transcript-resolved target list + codebook JSON + hash | lock panel identity, naming, and encoding |
| 3A. Candidate gen | `mkprobes run-panel` (batch) or `mkprobes candidates` | dataset path, target, output path | `*_all.parquet`, `*_bowtie.parquet`, `*_crawled.parquet` | broad candidate search + alignment annotation |
| 3B. Screening | `mkprobes screen` | candidate files + target | `*_screened_ol*.parquet` (+ stats) | filter and overlap-select high-quality probes |
| 3C. Construction | `mkprobes construct` | screened files + codebook | `*_final_*.parquet` | emit final encoded probe constructs |
| 4. Panel QC | `mkprobes filter-genes` | output path + target list | pass-list text output | enforce panel-level minimum probe count |
| 5. Manifest assembly | `mkprobes assemble` | manifest.json + final parquet files | assembled parquet/FASTA, `<panel>_final.txt`, provenance JSON | export orderable oligo pool |

## Key file names to watch

Per target `T`:

1. `output/T_crawled.parquet`
2. `output/T_screened_ol*.parquet`
3. `output/T_final_<restriction>_<bits>.parquet`

Panel-level:

1. `panel_a/codebook.json`
2. `panel_a/genes.converted.txt`
3. `panel_a/genes.pass.txt` (from `filter-genes`)

## Recommended run checkpoints

1. After phase 1: dataset index files exist and are readable.
2. After phase 2: codebook validation passes and hash recorded.
3. During phase 3: per-target final parquet appears for every target.
4. After phase 4: pass-list reaches project threshold.
5. After phase 5: assembled outputs and `<panel>.provenance.json` written under `generated/`.

## Common breakpoints and immediate checks

1. Missing `*_crawled.parquet`:
   - check `candidates` command args and dataset path.
2. Missing `*_screened_ol*.parquet`:
   - inspect `screen` constraints (`--minimum`, overlap, restriction filters).
3. Missing `*_final_*.parquet`:
   - verify target exists in codebook and screened inputs exist.
4. Low pass rate in `filter-genes`:
   - rerun selected targets with tuned screening parameters.
