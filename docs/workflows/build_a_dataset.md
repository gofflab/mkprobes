# Build a dataset

**Step 1 of the {doc}`../getting_started` workflow.**

A dataset is the sequence universe your probes are searched in and screened
against. Probe specificity is bounded by it: if the dataset is wrong or
incomplete, every later step produces confident-looking probes that bind the
wrong things. Build it once per species and release, then treat it as
immutable shared input.

There are three ways to build one. Pick by what you have:

| You have | Use | Section |
| --- | --- | --- |
| Nothing — you want mouse or human | `mkprobes prepare` | [A](#a-mouse-or-human) |
| A genome FASTA and a GTF/GFF3 annotation | `mkprobes ingest` | [B](#b-any-species-from-a-genome-and-annotation) |
| Transcript sequences already | `mkprobes create-dataset` | [C](#c-any-species-from-a-transcriptome-fasta) |

## A. Mouse or human

```bash
mkprobes prepare data --species mouse --threads 16
```

This downloads the curated reference set and indexes it. It writes into a
species subdirectory, so the above creates `data/mouse/`, and every later
command takes `data/mouse` — not `data`.

What it produces:

- `cdna_ncrna_trna.fasta`, the sequence universe
- a bowtie2 index (`txome`) for off-target search
- transcriptome 18-mers (`cdna18.jf`) and t/r/snoRNA exclusion 15-mers
  (`r_t_snorna15.jf`) from jellyfish
- `gencode.gtf.gz` and `ensembl.gtf.gz`

Set `--threads` to the CPUs you were actually allocated; oversizing it hurts
on a shared node.

Check it worked:

```bash
ls data/mouse/cdna18.jf data/mouse/r_t_snorna15.jf data/mouse/gencode.gtf.gz
```

## B. Any species, from a genome and annotation

```bash
mkprobes ingest data/myspecies \
    --genome genome.fa.gz --gtf annotation.gtf --species myspecies
```

`ingest` validates the annotation, extracts transcript sequences with
`gffread`, builds the bowtie2 and k-mer indices, and writes the dataset
manifest.

**Always run `--validate-only` first** and read the report. The failure it
catches most often is contig names that do not match between genome and
annotation (`chr1` vs `1` vs `scaffold_1`), which otherwise produces a
silently *empty* dataset rather than an error.

Flags you are likely to want:

- `--extract transcripts|cds` — spliced transcripts with UTRs (default,
  matching how the reference datasets are built) vs CDS only.
- `--rrna-fasta` / `--trna-fasta` — sequences for the probe blocklist, both
  repeatable. Without a blocklist, nothing stops probes landing in
  rRNA/tRNA-derived sequence.
- `--blocklist-biotypes rRNA,tRNA,snoRNA` — build the blocklist from the GTF's
  biotype column instead, when it has one.
- `--annotation-table NAME=PATH` — register a lookup table (orthologs,
  aliases, expression) so you can select targets by e.g. human ortholog
  symbol. Needs a `transcript_id` and/or `gene_id` column. Repeatable.
- `--keep-genome` — copy the genome into the dataset directory. Off by
  default because genomes are large; the sha256 is recorded either way.
- `--strip-version` — off by default here, and should stay off for de novo
  annotations: stripping `.N` suffixes merges StringTie isoforms
  (`STRG.1.1` and `STRG.1.2` become one).

Ingested datasets carry two extra files: `validation_report.json`, and
`solar_intake.yaml`, a provenance manifest with input sha256s, tool versions,
the literal command run and QC counts. **Fill in its stub fields** (assembly
source, annotation method, data owner) — that file is how another lab member
reproduces your dataset.

This is the path for non-traditional model species, and it has a full worked
runbook: {doc}`solar_new_species`.

## C. Any species, from a transcriptome FASTA

Use this only when you already have transcript sequences and no annotation to
extract them from. Prefer `ingest` when you have a genome and a GTF.

```bash
mkprobes create-dataset data/squid --fasta refs/squid_txome.fasta --species doryteuthis
```

Flags:

- `--gtf` — an annotation matching the FASTA. Optional, but supplying it is
  what enables gene-name lookups and sibling-isoform allowance during
  screening, so supply it if you have one.
- `--blocklist-fasta` — rRNA/tRNA sequences for the blocklist. Repeatable.
- `--annotation NAME=PATH` — register a lookup table (note: `--annotation`
  here, `--annotation-table` on `ingest`). Repeatable.
- `--fasta-key-regex` — how to extract the record ID from FASTA headers.
- `--strip-version` / `--no-strip-version` — on by default here, unlike
  `ingest`. Use `--no-strip-version` for de novo annotations whose IDs embed
  meaningful dots (`STRG.1.1`, `g1.t1`).

`--species` is metadata only. It is not checked against the FASTA, and
`human` and `mouse` are reserved for reference datasets.

## Reference or custom: why it matters

The two kinds of dataset are screened differently, and `mkprobes` decides
which you have from the directory contents. Reference datasets screen
candidates against pseudogenes and treat every Ensembl isoform of a target as
an acceptable binder; custom datasets do neither. Keep the two in separate
directories — a directory holding both is rejected rather than guessed at.

The exact rules, and the fields inside `dataset.json`, are in
[Which kind of dataset a command loads](../reference/file_formats.md#which-kind-of-dataset-a-command-loads).

## On a cluster

- Build once per species and version, into an explicit versioned root
  (`refs/mouse_release_110/`), then treat it as read-only.
- Prefer one high-memory prep job over several duplicate ones.
- Record provenance: `urls.tsv` for reference mode, `solar_intake.yaml` for
  ingested datasets.

## When it goes wrong

- **Downloads fail** — check outbound network access from the node you are on.
- **`jellyfish`, `bowtie2` or `gffread`: command not found** — install them
  into the environment you are running from; all three are on bioconda. See
  {doc}`../before_you_start`.
- **Species not supported by `prepare`** — that command is mouse and human
  only. Use `ingest` or `create-dataset`.
- **Stale files after a release change** — build into a clean, versioned
  directory rather than over the top of the old one.

---

Next: {doc}`choose_your_targets`.
