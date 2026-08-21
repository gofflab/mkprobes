# SOLAR probesets for a new species

This is the end-to-end runbook for creating a **SOLAR** probeset (the lab's
splint/padlock, STARmap-style combinatorial FISH assay) for a
**non-traditional model species**, starting from nothing but a genome FASTA
and an annotation (GTF or GFF3). Every command below was validated live on
*Octopus chierchiae* (StringTie-merged annotation, 1.3 GB genome, 53,045
transcripts); the example outputs shown are from that run.

The pipeline is annotation-robust by design: it handles de novo annotations
from StringTie, AUGUSTUS/BRAKER, MAKER, Trinity-derived GTFs, and
NCBI-Gnomon-style files, including IDs containing dots and underscores
(`STRG.1.1`, `g1.t1`, `TRINITY_DN123_c0_g1_i1`, `rna-XM_012345.1`).

## Overview

```text
genome.fa(.gz) + annotation.gtf/.gff3(.gz)
        │  mkprobes ingest          validate → gffread → indices → manifest
        ▼
dataset directory (dataset.json, transcripts.fasta, bowtie2 + k-mer indices)
        │  mkprobes transcripts     pick target transcripts (--longest)
        ▼
        │  mkprobes make-codebook   assign readout bits (optionally expression-informed)
        ▼
        │  mkprobes run-panel       candidates → screen → construct, all targets in parallel
        ▼
        │  mkprobes filter-genes    panel QC
        ▼
mkprobes assemble gen   →  orderable oligo pool
```

## 0. Prerequisites

- Python ≥ 3.12, this package installed and its environment activated
  (see [Installation](../installation.md) — uv or conda/mamba).
- External tools on `PATH`: **gffread**, **bowtie2**, **jellyfish**
  (all on bioconda; `mkprobes ingest` checks them up front and prints exact
  install hints for anything missing). **RepeatMasker** is optional and only
  used at final assembly.

Inputs you need:

| Input | Required | Notes |
|---|---|---|
| Genome FASTA | yes | `.gz` fine. Contig names must match the GTF `seqname` column — validated. |
| Annotation | yes | GTF or GFF3, `.gz` fine. GFF3 is converted with `gffread -T` automatically. |
| Species name | yes | Free-form metadata (not `human`/`mouse` — those are reserved for reference datasets). |
| rRNA/tRNA FASTA(s) | recommended | For the probe blocklist. Sources: [SILVA](https://www.arb-silva.de/) (rRNA), [GtRNAdb](http://gtrnadb.ucsc.edu/) or `tRNAscan-SE` output (tRNA), or `--blocklist-biotypes` if your GTF carries biotype attributes. **Without a blocklist, nothing stops probes from landing in rRNA/tRNA-derived sequence.** |
| Annotation tables | optional | Ortholog/alias/expression tables (parquet/csv/tsv with a `transcript_id` and/or `gene_id` column). Lets you select targets by e.g. human ortholog symbol. |

## 1. Ingest: genome + annotation → dataset

Always run validation first and read the report:

```bash
mkprobes ingest data/ochierchiae \
    --genome refs/Ochierchiae_genome.fa.gz \
    --gtf refs/Ochier_stringtie_merged.gtf \
    --species octopus_chierchiae \
    --validate-only
```

The validation report (also saved as `validation_report.json`) checks, with
named error codes and a fix for each:

- **`SEQNAME_MISMATCH`** — GTF seqnames absent from the genome (the classic
  `chr1` vs `1` vs `scaffold_1` mismatch, which otherwise yields a silently
  *empty* gffread output).
- **`DUPLICATE_TRANSCRIPT_ID`**, **`REQUIRED_ATTR_MISSING`** (gene_id /
  transcript_id), **`NO_EXON_ROWS`**, **`COORDINATES_INVERTED`**,
  **`ID_FORBIDDEN_CHARS`** — hard errors; the dataset is not built.
- **`STRAND_MISSING`**, **`TRANSCRIPT_WITHOUT_EXONS`**, **`ID_COLON`**,
  **`GENE_NAME_FALLBACK`** — warnings. (The octopus run had 21,500 unstranded
  rows; gffread handled them identically to the lab's reference extraction.)

It also reports gene/transcript/isoform counts and which **gene-name source**
applies. De novo GTFs usually have no `gene_name`; names then fall back
`gene_name ← Name ← gene ← gene_id`, meaning targets are selected by ID
unless you register an ortholog table (below).

Then run the full ingest:

```bash
mkprobes ingest data/ochierchiae \
    --genome refs/Ochierchiae_genome.fa.gz \
    --gtf refs/Ochier_stringtie_merged.gtf \
    --species octopus_chierchiae \
    --rrna-fasta refs/octopus_rrna.fasta \
    --trna-fasta refs/octopus_trna.fasta \
    --annotation-table orthologs=refs/human_orthologs.tsv
```

Useful flags:

- `--extract transcripts|cds` — gffread `-w` (spliced transcripts with UTRs,
  default, matches how reference cDNA datasets are built) vs `-x` (CDS only).
- `--blocklist-biotypes rRNA,tRNA,snoRNA` — auto-extract blocklist sequences
  by the GTF's biotype column, when it has one (StringTie GTFs don't).
- `--keep-genome` — copy the genome into the dataset dir (needed later only
  for genome-mode simulation; off by default because genomes are large;
  provenance sha256 is recorded either way).
- `--strip-version` — off by default and should stay off for de novo
  annotations: gffread FASTA headers match the GTF `transcript_id` verbatim,
  and stripping `.N` suffixes would merge StringTie isoforms
  (`STRG.1.1`/`STRG.1.2`).

What lands in the dataset directory:

```text
data/ochierchiae/
├── dataset.json               # machine-readable definition (see file_formats)
├── solar_intake.yaml          # provenance manifest — COMPLETE THE STUBS
├── validation_report.json
├── annotation.gtf             # normalized plain-text GTF
├── transcripts.fasta          # gffread output (53,045 records for octopus)
├── transcripts.parquet        # parsed GTF cache
├── transcripts.{1..4,rev.1,rev.2}.bt2   # bowtie2 index
├── transcripts.jf             # 18-mer counts
└── blocklist15.jf             # 15-mer rRNA/tRNA blocklist (if provided)
```

`solar_intake.yaml` auto-fills input sha256s, tool versions, the literal
command, and QC counts. **Fill in the provenance stubs** (assembly source,
annotation method, data owner) — the manifest is the record another lab
member uses to reproduce your dataset.

Sanity check that everything round-trips (ingest already did this, but it's
cheap to repeat after any manual change):

```bash
mkprobes transcripts data/ochierchiae --gene <any_gene_id> --longest
```

If you also have an independently produced transcriptome FASTA, compare it
against `transcripts.fasta` — for octopus, gffread's output was 100.00%
identical to the lab's reference transcriptome (all 53,045 records).

## 2. Pick target transcripts

For custom datasets, transcript selection is offline (no Ensembl):

```bash
# one gene -> its longest isoform (default for custom datasets)
mkprobes transcripts data/ochierchiae --gene Och.576 --longest

# a file of genes/IDs -> genes.tss.txt
mkprobes convert-to-transcripts data/ochierchiae genes.txt -m longest
```

- `--longest` picks, per gene, the isoform with the longest **sequence**
  (measured from the FASTA, so introns don't distort the choice).
- `--all` returns every isoform.
- Tokens may be transcript IDs (`Och.576.10`, passed through), gene IDs
  (`Och.576`, expanded), or — if you registered an ortholog/alias table —
  external symbols (`SHANK3`, case-insensitive).
- Unresolvable tokens fail loudly with close-match suggestions.

`mkprobes chkgenes data/ochierchiae genes.txt` performs the same resolution
as a standalone validation step and writes `genes.converted.txt`.

## 3. Generate the codebook

The codebook maps each target to three distinct readout bits
(`{"Och.687.1": [1, 2, 3], ...}`):

```bash
mkprobes make-codebook data/ochierchiae genes.tss.txt -o codebook.json
```

This auto-sizes an MHD code, assigns codewords (seeded), fills spare
capacity with `Blank-N` decoys, and logs the codebook hash. **Optionally**,
if expression data exists (a table registered at ingest, or any
parquet/csv/tsv on disk), add `--expression <name-or-path>` to balance
total expression load across readout bits — without it, the plain seeded
assignment is used. Details: [Phase 2](phase_2_codebook_design.md).

## 4. Design probes for the whole panel

One command runs `candidates → screen → construct` for every target in the
codebook, in parallel:

```bash
mkprobes run-panel data/ochierchiae codebook.json
```

Finished targets are skipped on re-runs; failures are isolated per gene
(collected in `codebook.failed.txt` — triage with `--list-failed-all`); an
off-target allow-list at `codebook.acceptable.json` is applied automatically.
Re-run one target: `mkprobes run-panel data/ochierchiae codebook.json Och.687.1`.

For debugging a single target step-by-step, the underlying commands are:

```bash
mkprobes candidates data/ochierchiae -g Och.687.1 -o output/
mkprobes screen output/ Och.687.1 --minimum 60 --maxoverlap 0 --restriction BamHI,KpnI
mkprobes construct data/ochierchiae output/ -g Och.687.1 \
    -c codebook.json --restriction BamHI --restriction KpnI
```

Notes specific to custom datasets:

- **Sibling isoforms are auto-allowed**: probes for `Och.576.10` are not
  penalized for binding `Och.576.1`–`.9` (the other isoforms of the same
  gene, looked up from the GTF). This mirrors the reference-dataset behavior;
  without it, any multi-isoform gene would yield zero probes.
- `--allow`/`--disallow` take **transcript IDs** (FASTA record IDs), not gene
  names. Use them for known cross-reactive homologs surfaced in the
  "Most common binders" log table or `<target>_offtarget_counts.csv`.
- The rRNA/tRNA blocklist is enforced automatically when the dataset has one
  (a single warning is printed when it doesn't).

Octopus reference numbers (17.9 kb / 16.6 kb / 16.0 kb targets): 1,700–2,800
candidates each → 57–73 screened pairs → 54–69 constructed probes.

## 5. Panel QC

```bash
mkprobes filter-genes output/ --genes genes.tss.txt \
    --min-probes 48 --out genes.pass.txt
```

Genes below `--min-probes` are warned about individually; `genes.pass.txt`
gets everything at or above the threshold. For low-yield targets: try a
different isoform, loosen `screen --minimum/--maxoverlap`, or accept verified
homologous off-targets via `--allow`.

## 6. Final assembly → orderable oligos

Write a manifest (list of probesets) next to your codebook:

```json
[
  {
    "name": "och_panel_v1",
    "species": "octopus_chierchiae",
    "codebook": "codebook.json",
    "bcidx": 0,
    "n_probes": 16
  }
]
```

Then assemble:

```bash
# non-model species: either give RepeatMasker a supported taxon...
mkprobes assemble manifest.json gen --rm-species mollusca
# ...or skip it explicitly
mkprobes assemble manifest.json gen --skip-repeatmasker
```

Outputs under `generated/`: `<name>.parquet`, `<name>_pad.fasta`,
`<name>_splint.fasta`, the orderable pool `<name>_final.txt`
(one oligo per line), and `<name>.provenance.json` (timestamp, mkprobes
version, codebook hash, probeset config, RepeatMasker status). Every pair is
assert-checked for splint/padlock geometry and 139–150 nt padlock length
during assembly.

The interactive off-target triage (`mkprobes assemble manifest.json short <N>`)
reviews low-count genes and writes accepted off-targets to
`<codebook>.acceptable.json`, which `mkprobes run-panel` picks up
automatically on the next run. Assembly output is deterministic: the same
inputs always produce the same oligo pool.

## Troubleshooting quick hits

- **gffread output empty / `SEQNAME_MISMATCH`** — the annotation was not
  built for this assembly, or contigs were renamed. Fix the names; don't mix
  assemblies.
- **"looks like GFF3"** — pass the GFF3 straight to `mkprobes ingest`
  (it converts), or `gffread in.gff3 -T -o out.gtf`.
- **`Could not resolve 'X'`** in transcript selection — you used a symbol the
  annotation doesn't know. Use the dataset's own IDs, or register an
  ortholog/alias table at ingest.
- **Zero probes for a multi-isoform gene** — should not happen (siblings are
  auto-allowed); if it does, inspect `<target>_offtarget_counts.csv` for a
  homolog and `--allow` it after verifying.
- **StringTie isoforms collapse into one ID** — the dataset was built with
  `--strip-version`. Rebuild with the default (`--no-strip-version`).
