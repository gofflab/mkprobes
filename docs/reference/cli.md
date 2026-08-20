# CLI reference

This page summarizes the `mkprobes` command surface and where each command fits. The package lives in its own repository ([github.com/gofflab/mkprobes](https://github.com/gofflab/mkprobes)); run commands with `uv run mkprobes ...` from the repo root.

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

### Probe generation

- `mkprobes candidates PATH --gene GENE --output output/ [options]`
- `mkprobes screen DATA_DIR GENE [--minimum M] [--overlap L] [--maxoverlap L] [--restriction ...] [--overwrite]`
- `mkprobes construct PATH OUTPUT_PATH --gene GENE --codebook CODEBOOK.json [options]`

### Panel QC and provenance

- `mkprobes filter-genes OUTPUT_PATH --genes genes.txt --min-probes N [--out out.txt]`
- `mkprobes hash CODEBOOK.json`

## Probe-generation scripts

Manifest assembly and final export are driven by `scripts/probegen/2_assemble_manifest.py` (see {doc}`../workflows/phase_5_manifest_assembly`):

```bash
uv run python scripts/probegen/2_assemble_manifest.py manifest.json short 12
uv run python scripts/probegen/2_assemble_manifest.py manifest.json gen
```

Other scripts under `scripts/probegen/` (codebook generation, simulation) are lab-specific orchestration; for documentation-driven workflows, use the `mkprobes` CLI commands directly.
