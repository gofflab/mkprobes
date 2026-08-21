# CLI reference

This page summarizes the `mkprobes` command surface and where each command fits. The package lives in its own repository ([github.com/gofflab/mkprobes](https://github.com/gofflab/mkprobes)); run commands with `mkprobes ...` from the repo root.

For an end-to-end walkthrough on a non-model species, see {doc}`../workflows/solar_new_species`.

## `mkprobes` command group

### Dataset preparation

- `mkprobes prepare PATH --species {human,mouse} --threads N`
  - downloads and indexes the curated human/mouse reference.
- `mkprobes ingest PATH --genome GENOME.fa --gtf ANNOTATION.gtf --species NAME [options]`
  - builds a dataset for any species from a genome FASTA + GTF/GFF3.
  - options: `--extract {transcripts,cds}`, `--rrna-fasta`, `--trna-fasta`, `--blocklist-biotypes`, `--annotation-table NAME=PATH`, `--keep-genome`, `--fasta-key-regex`, `--strip-version/--no-strip-version` (default: no strip), `--validate-only`, `--overwrite`.
- `mkprobes create-dataset PATH --fasta FASTA --species NAME [options]`
  - builds a dataset directly from a transcriptome FASTA.
  - options: `--gtf`, `--blocklist-fasta`, `--fasta-key-regex`, `--strip-version/--no-strip-version` (default: strip), `--annotation NAME=PATH`, `--overwrite`.

### Target selection

- `mkprobes chkgenes PATH genes.txt`
  - validate/normalize gene names against the dataset.
- `mkprobes transcripts PATH --gene GENE [--canonical|--gencode|--ensembl|--appris|--longest|--all]`
  - `--longest`/`--all` are the offline modes for custom datasets; `--canonical` falls back to `--longest` on custom datasets. Ensembl/mygene/APPRIS modes are reference-only (human/mouse).
- `mkprobes convert-to-transcripts PATH genes.txt --mode MODE`
  - `-m/--mode` accepts `canonical`, `gencode`, `ensembl`, `appris`, `apprisalt`, `longest`, `all`.

### Codebook generation

- `mkprobes make-codebook PATH genes.tss.txt [options]`
  - generates a codebook JSON from an MHD code (auto-sized to the target count, `Blank-N` decoys fill spare capacity; round-confounder codewords are swapped onto blanks).
  - `--expression NAME_OR_PATH` (**optional**): balance total expression load across readout bits by trying `--iterations` assignments and keeping the highest-entropy one. Accepts the name of an annotation table registered in the dataset or a parquet/csv/tsv path (needs a `transcript_id`/`gene_id` column; `--expression-column` disambiguates the value column). Without it, assignment is a plain seeded shuffle (`--seed`).
  - `--existing-codebook prev.json` extends a previous panel with non-overlapping bits; `--n-bits`/`--offset` for manual control.

### Probe generation

- `mkprobes run-panel PATH CODEBOOK [GENE] [options]`
  - the batch driver: designs probes for every target in the codebook (`candidates -> screen -> construct`) across parallel workers; skips finished targets, applies `<codebook>.acceptable.json` (or `--allow-file`), records failures in `<codebook>.failed.txt`, and exits non-zero if any gene fails. `--list-failed`/`--list-failed-all` triage missing outputs. Production defaults: `--minimum 60 --maxoverlap 0 --restriction BamHI,KpnI --target-probes 48`.
- `mkprobes candidates PATH --gene GENE --output output/ [options]`
- `mkprobes screen DATA_DIR GENE [--minimum M] [--overlap L] [--maxoverlap L] [--restriction ...] [--overwrite]`
- `mkprobes construct PATH OUTPUT_PATH --gene GENE --codebook CODEBOOK.json [options]`

### Panel QC and provenance

- `mkprobes filter-genes OUTPUT_PATH --genes genes.txt --min-probes N [--out out.txt]`
- `mkprobes hash CODEBOOK.json`

## Manifest assembly

Manifest assembly and final export (see {doc}`../workflows/phase_5_manifest_assembly`):

```bash
mkprobes assemble manifest.json short 12     # triage under-provisioned genes (interactive off-target accept)
mkprobes assemble manifest.json gen          # assemble the orderable oligo pool (deterministic)
```

`gen` accepts `--rm-species '<taxon>'` / `--skip-repeatmasker` for non-model species and `--headerfooter` to override the vendored table.

The remaining files under `scripts/probegen/` are deprecated shims (`o_codebook.py` -> `mkprobes make-codebook`; `1_run_codebook*.py` -> `mkprobes run-panel`; `2_assemble_manifest.py` -> `mkprobes assemble`) plus exploratory notebooks (`simulate.py` in-silico validation, `foridt.py`/`adt.py` IDT-ordering examples whose logic already lives in the package).
