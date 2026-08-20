# Phase 1: Dataset preparation

Phase 1 produces the indexed reference data required for candidate search and off-target filtering when designing SOLAR (splint/padlock, STARmap-style) probesets.

## Why this phase matters

Probe design quality is bounded by reference quality. If this phase is wrong or incomplete, later phases can produce misleading off-target calls or low-yield candidates. Treat dataset prep as foundational infrastructure, not an optional pre-step.

## Paths supported

1. Reference mode (`human`/`mouse`) with `mkprobes prepare`.
2. Ingest mode (genome FASTA + GTF/GFF3, any species) with `mkprobes ingest`.
3. Generic mode (transcriptome FASTA) with `mkprobes create-dataset`.

## 1A. Reference mode (recommended for human/mouse)

```bash
uv run mkprobes prepare <DATA_ROOT> --species mouse --threads 16
```

Example:

```bash
uv run mkprobes prepare data --species human --threads 24
```

### What it does

- Downloads required reference files listed in species URL TSV.
- Creates `cdna_ncrna_trna.fasta`.
- Builds bowtie index (`txome`).
- Runs jellyfish for:
  - transcriptome 18-mers (`cdna18.jf`)
  - t/r/snoRNA exclusion kmers (`r_t_snorna15.jf`)
- Validates dataset by constructing `ReferenceDataset`.

### Key arguments explained

- `PATH`:
  - **Intention:** root folder where species subfolders and derived indices are stored.
  - **Use carefully because:** changing this path changes all downstream command inputs.
- `--species`:
  - **Intention:** selects curated reference source set (`human` or `mouse`).
  - **Use carefully because:** species mismatch causes invalid transcript mapping/candidate behavior.
- `--threads`:
  - **Intention:** parallelism for indexing/k-mer generation.
  - **Use carefully because:** oversizing can hurt shared-node performance; align with allocated CPU.

### Expected output location

`<DATA_ROOT>/<species>/`, e.g. `data/mouse/`.

## 1B. Ingest mode (any species, genome + annotation)

```bash
uv run mkprobes ingest <DATASET_DIR> --genome genome.fa.gz --gtf annotation.gtf --species <SPECIES_NAME>
```

`ingest` validates the annotation, extracts transcript (or CDS) sequences with `gffread`, builds the bowtie2 and k-mer indices, and writes the dataset manifest. Useful flags include `--extract transcripts|cds`, `--rrna-fasta`/`--trna-fasta` (blocklist inputs), `--blocklist-biotypes`, `--annotation-table NAME=PATH`, `--keep-genome`, `--fasta-key-regex`, `--strip-version/--no-strip-version` (default: no strip), `--validate-only`, and `--overwrite`.

This is the recommended path for non-traditional model species; the full runbook is {doc}`solar_new_species`.

Ingested datasets additionally contain:

- `solar_intake.yaml`: provenance manifest (input sha256s, tool versions, the literal command run, QC counts). The operator completes the stub fields.
- `validation_report.json`: annotation/sequence validation results.

## 1C. Generic/custom dataset mode (transcriptome FASTA)

```bash
uv run mkprobes create-dataset <DATASET_DIR> --fasta <REFERENCE_FASTA> --species <SPECIES_NAME>
```

Example:

```bash
uv run mkprobes create-dataset data/squid --fasta refs/squid_txome.fasta --species doryteuthis
```

This writes `dataset.json` and required index/k-mer files in the target folder.

### Key arguments explained

- `PATH`:
  - **Intention:** output folder for generic dataset files.
- `--fasta`:
  - **Intention:** sequence space to search during candidate generation.
  - **Use carefully because:** wrong FASTA means all downstream probes target the wrong universe.
- `--species`:
  - **Intention:** metadata label used for reproducibility/reporting.
  - **Use carefully because:** does not automatically validate taxonomy compatibility with FASTA.
- `--gtf`:
  - **Intention:** optional annotation used for gene/isoform relationships (e.g. sibling-isoform allowance in phase 3).
- `--blocklist-fasta`:
  - **Intention:** rRNA/tRNA sequences to build a probe blocklist k-mer database from.
- `--fasta-key-regex`:
  - **Intention:** regex for extracting record IDs from FASTA headers; persisted in `dataset.json`.
- `--strip-version/--no-strip-version`:
  - **Intention:** whether to strip trailing version suffixes from record IDs (default: strip); persisted in `dataset.json`.
- `--annotation NAME=PATH`:
  - **Intention:** register a lookup table (e.g. orthologs) for symbol resolution in transcript selection; tables must carry `transcript_id` and/or `gene_id` columns.

### Custom dataset manifest (`dataset.json`)

Generic/custom datasets support these fields (all newer fields are optional and back-compatible):

- `gtf_name` / `cache_name`: annotation file and its parsed cache.
- `blocklist_kmer_name`: rRNA/tRNA blocklist k-mer database (typically `blocklist15.jf`); enforced automatically in phase 3 when present.
- `genome_fasta_name`: retained genome FASTA (from `ingest --keep-genome`).
- `annotations`: registry of named lookup tables, e.g. `{"orthologs": "orthologs.tsv"}`.
- `fasta_key_regex`, `strip_version`: header-parsing settings persisted per dataset.

## Validation checks

Run basic checks after prep:

```bash
test -f data/mouse/cdna18.jf
test -f data/mouse/r_t_snorna15.jf
```

And (for reference mode):

```bash
test -f data/mouse/gencode.gtf.gz
test -f data/mouse/ensembl.gtf.gz
```

## HPC add-on guidance

- Run dataset prep once per species/version, then treat as immutable shared input.
- Use explicit versioned roots (example: `refs/mouse_release_110/`).
- Prefer fewer high-memory prep jobs over many duplicate prep jobs.
- Document provenance: the source URL manifest (`urls.tsv`) for reference mode, or `solar_intake.yaml` for ingested datasets.

## Failure modes

- Missing downloaded reference files:
  - verify outbound network access on setup node.
- `jellyfish`, `bowtie`, or `gffread` command not found:
  - install tools in the same environment used for prep (all available on bioconda).
- Species not supported in reference mode:
  - use `ingest` (genome + annotation) or generic `create-dataset` (transcriptome FASTA).
- Reusing stale dataset files:
  - re-run prep in a clean/versioned directory when changing genome/transcript release.
