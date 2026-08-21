# Getting Started

The `mkprobes` workflow, start to finish. Each step is one command; every
step validates the prerequisites of the next. Deep dives are linked per step.

```text
1. dataset    mkprobes prepare | ingest      reference (mouse/human) | any species
2. targets    mkprobes chkgenes / convert-to-transcripts
3. codebook   mkprobes make-codebook
4. probes     mkprobes run-panel             candidates -> screen -> construct, all targets in parallel
5. panel QC   mkprobes filter-genes
6. assembly   mkprobes assemble              -> orderable oligos
```

Before you start: installation complete ({doc}`installation`), `mkprobes --help`
works, and `bowtie2`/`jellyfish`/`gffread` are on `PATH`.

## Suggested project layout

```text
project/
├── data/
│   └── mouse/                      # dataset (prepare or ingest)
├── panel_a/
│   ├── genes.txt                   # your target list
│   ├── codebook.json               # from make-codebook
│   └── output/                     # per-target parquet outputs
```

## 1. Build the dataset (once per species)

Mouse or human reference:

```bash
mkprobes prepare data --species mouse --threads 16
```

Any other species, from a genome FASTA + GTF/GFF3:

```bash
mkprobes ingest data/myspecies --genome genome.fa.gz --gtf annotation.gtf --species myspecies
```

Run `ingest` with `--validate-only` first and read the report. Full walkthrough:
{doc}`workflows/solar_new_species`.

## 2. Validate targets

`panel_a/genes.txt`, one target per line (gene names for reference datasets;
gene or transcript IDs for custom ones):

```bash
mkprobes chkgenes data/mouse panel_a/genes.txt
mkprobes convert-to-transcripts data/mouse panel_a/genes.converted.txt   # -m longest for custom datasets
```

This normalizes names and resolves each gene to a transcript
(`genes.converted.tss.txt`). Details: {doc}`workflows/phase_2_codebook_design`.

## 3. Generate the codebook

```bash
mkprobes make-codebook data/mouse panel_a/genes.converted.tss.txt -o panel_a/codebook.json
```

Auto-sizes the code, assigns three readout bits per target, fills spare
capacity with blanks, and logs the codebook hash. Have expression data?
Add `--expression <table>` to balance load across bits — optional.
Details: {doc}`workflows/phase_2_codebook_design`.

## 4. Design probes (all targets)

```bash
mkprobes run-panel data/mouse panel_a/codebook.json
```

The codebook is the work list: every target runs
`candidates -> screen -> construct` in parallel (16 workers by default),
finished targets are skipped on re-runs, failures are isolated per gene and
collected into `codebook.failed.txt`, and a `codebook.acceptable.json`
allow-list (from off-target triage) is applied automatically. Re-run a single
target with `mkprobes run-panel data/mouse panel_a/codebook.json <gene>`.
Details and the underlying single-target commands:
{doc}`workflows/phase_3_candidate_screen_construct`.

## 5. Panel QC

```bash
mkprobes filter-genes panel_a/output --genes panel_a/genes.converted.tss.txt \
    --min-probes 48 --out panel_a/genes.pass.txt
```

Targets below the threshold are warned individually; rework or drop them
before assembly. Details: {doc}`workflows/phase_4_panel_qc_export`.

## 6. Assemble the orderable pool

```bash
mkprobes assemble panel_a/manifest.json gen
```

Emits the final oligo pool plus a provenance sidecar.
Details: {doc}`workflows/phase_5_manifest_assembly`.

## Expected outputs per target

`*_all.parquet` → `*_bowtie.parquet` → `*_crawled.parquet` →
`*_screened_ol*.parquet` → `*_final_BamHIKpnI_<bits>.parquet`

## Next steps

- New-species end-to-end runbook: {doc}`workflows/solar_new_species`
- How the assay and algorithms work: {doc}`workflows/assay_and_under_the_hood`
- Input/output flow map: {doc}`workflows/dataflow_map`
- Command details: {doc}`reference/cli`
