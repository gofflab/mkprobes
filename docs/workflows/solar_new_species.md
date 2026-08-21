# SOLAR probesets for a new species

Everything that is *different* when your species is not mouse or human.

This page is a companion to {doc}`../getting_started`, not a replacement for
it. The workflow is the same six steps in the same order; read that page for
the steps themselves and come back here for the deltas. Every command below
was validated live on *Octopus chierchiae* — a StringTie-merged annotation,
1.3 GB genome, 53,045 transcripts — and the numbers quoted are from that run.

The pipeline is annotation-robust by design: StringTie, AUGUSTUS/BRAKER,
MAKER, Trinity-derived GTFs and NCBI-Gnomon-style files all work, including
IDs containing dots and underscores (`STRG.1.1`, `g1.t1`,
`TRINITY_DN123_c0_g1_i1`, `rna-XM_012345.1`).

## What changes

| Step | Reference species | Your species |
| --- | --- | --- |
| 1. dataset | `mkprobes prepare` downloads it | `mkprobes ingest` builds it from your genome + GTF |
| 2. targets | gene names, checked against Ensembl | dataset IDs, resolved offline; `-m longest` |
| 3. codebook | unchanged | unchanged |
| 4. probes | unchanged | `--allow`/`--disallow` take transcript IDs |
| 5. panel QC | unchanged | unchanged |
| 6. assembly | RepeatMasker taxon inferred | `--rm-species` or `--skip-repeatmasker` |

## What you need before you start

| Input | Required | Notes |
| --- | --- | --- |
| Genome FASTA | yes | `.gz` fine. Contig names must match the GTF `seqname` column — validated. |
| Annotation | yes | GTF or GFF3, `.gz` fine. GFF3 is converted with `gffread -T` automatically. |
| Species name | yes | Free-form metadata. Not `human`/`mouse` — those are reserved for reference datasets. |
| rRNA/tRNA FASTA(s) | recommended | For the probe blocklist. [SILVA](https://www.arb-silva.de/) for rRNA, [GtRNAdb](http://gtrnadb.ucsc.edu/) or `tRNAscan-SE` output for tRNA. **Without a blocklist, nothing stops probes landing in rRNA/tRNA-derived sequence.** |
| Annotation tables | optional | Ortholog/alias/expression tables (parquet/csv/tsv with `transcript_id` and/or `gene_id`). Lets you pick targets by human ortholog symbol. |

Disk, memory and time expectations: {doc}`../before_you_start`.

## 1. Ingest, and read the validation report

Always validate first:

```bash
mkprobes ingest data/ochierchiae \
    --genome refs/Ochierchiae_genome.fa.gz \
    --gtf refs/Ochier_stringtie_merged.gtf \
    --species octopus_chierchiae \
    --validate-only
```

The report (also written to `validation_report.json`) uses named codes, each
with a fix:

- **`SEQNAME_MISMATCH`** — GTF seqnames absent from the genome: the classic
  `chr1` vs `1` vs `scaffold_1` mismatch. This is the one that matters most,
  because without the check it produces a silently *empty* gffread output
  rather than an error.
- **`DUPLICATE_TRANSCRIPT_ID`**, **`REQUIRED_ATTR_MISSING`**, **`NO_EXON_ROWS`**,
  **`COORDINATES_INVERTED`**, **`ID_FORBIDDEN_CHARS`** — hard errors; nothing
  is built.
- **`STRAND_MISSING`**, **`TRANSCRIPT_WITHOUT_EXONS`**, **`ID_COLON`**,
  **`GENE_NAME_FALLBACK`** — warnings. (The octopus run had 21,500 unstranded
  rows; gffread handled them identically to the lab's reference extraction.)

It also reports which **gene-name source** applies. De novo GTFs usually carry
no `gene_name`, so names fall back `gene_name ← Name ← gene ← gene_id` —
meaning you select targets by ID unless you register an ortholog table.

Then run it for real:

```bash
mkprobes ingest data/ochierchiae \
    --genome refs/Ochierchiae_genome.fa.gz \
    --gtf refs/Ochier_stringtie_merged.gtf \
    --species octopus_chierchiae \
    --rrna-fasta refs/octopus_rrna.fasta \
    --trna-fasta refs/octopus_trna.fasta \
    --annotation-table orthologs=refs/human_orthologs.tsv
```

Flag-by-flag detail is in {doc}`build_a_dataset`. The one worth repeating:
**leave `--strip-version` off** (its default here). Stripping `.N` suffixes
merges StringTie isoforms — `STRG.1.1` and `STRG.1.2` collapse into one.

What lands in the dataset directory:

```text
data/ochierchiae/
├── dataset.json               # machine-readable definition
├── solar_intake.yaml          # provenance manifest — COMPLETE THE STUBS
├── validation_report.json
├── annotation.gtf             # normalized plain-text GTF
├── transcripts.fasta          # gffread output (53,045 records for octopus)
├── transcripts.parquet        # parsed GTF cache
├── transcripts.{1..4,rev.1,rev.2}.bt2   # bowtie2 index
├── transcripts.jf             # 18-mer counts
└── blocklist15.jf             # 15-mer rRNA/tRNA blocklist (if provided)
```

**Fill in the `solar_intake.yaml` stubs** (assembly source, annotation method,
data owner). The sha256s, tool versions, literal command and QC counts are
filled in for you; the stubs are what another lab member needs to reproduce
your dataset.

Sanity check that it round-trips:

```bash
mkprobes transcripts data/ochierchiae --gene <any_gene_id> --longest
```

If you have an independently produced transcriptome FASTA, compare it against
`transcripts.fasta`. For octopus, gffread's output was 100.00% identical to
the lab's reference transcriptome across all 53,045 records.

## 2. Targets, offline

No Ensembl, no network:

```bash
mkprobes convert-to-transcripts data/ochierchiae genes.txt -m longest
```

Pass `-m longest` explicitly. Tokens may be transcript IDs (`Och.576.10`,
passed through), gene IDs (`Och.576`, expanded to its isoforms), or — if you
registered an ortholog table — external symbols (`SHANK3`, case-insensitive).
Unresolvable tokens fail with close-match suggestions.

The rest of {doc}`choose_your_targets` applies unchanged.

## 3. Codebook

Nothing species-specific. See {doc}`design_the_codebook`.

## 4. Probes

```bash
mkprobes run-panel data/ochierchiae codebook.json
```

Three things behave differently on a custom dataset:

- **Sibling isoforms are auto-allowed.** Probes for `Och.576.10` are not
  penalized for binding `Och.576.1`–`.9`, looked up from the GTF. This mirrors
  reference-dataset behaviour; without it, every multi-isoform gene would
  yield zero probes.
- **`--allow`/`--disallow` take transcript IDs**, not gene names — the FASTA
  record IDs. Get them from the "Most common binders" log table or
  `<target>_offtarget_counts.csv`.
- **The rRNA/tRNA blocklist is applied automatically** when the dataset has
  one. A single warning is printed when it does not.

Octopus reference numbers, for 17.9 / 16.6 / 16.0 kb targets: 1,700–2,800
candidates each → 57–73 screened pairs → 54–69 constructed probes.

## 5. Panel QC

Unchanged; see {doc}`qc_your_panel`. For low-yield targets the non-model
options are the usual ones: try a different isoform, loosen
`--minimum`/`--maxoverlap`, or accept verified homologous off-targets.

## 6. Assembly

The only difference is RepeatMasker, which has no built-in taxon mapping for
your species:

```bash
mkprobes assemble manifest.json gen --rm-species mollusca   # a taxon its library knows
mkprobes assemble manifest.json gen --skip-repeatmasker     # or skip it explicitly
```

Everything else — manifest fields, `short` triage, outputs — is
{doc}`order_your_oligos`.

## Troubleshooting, species-specific

- **gffread output empty / `SEQNAME_MISMATCH`** — the annotation was not built
  for this assembly, or contigs were renamed. Fix the names; do not mix
  assemblies.
- **"looks like GFF3"** — pass the GFF3 straight to `mkprobes ingest`, which
  converts it, or convert manually with `gffread in.gff3 -T -o out.gtf`.
- **`Could not resolve 'X'`** — you used a symbol the annotation does not
  know. Use the dataset's own IDs, or register an ortholog table at ingest.
- **StringTie isoforms collapsed into one ID** — the dataset was built with
  `--strip-version`. Rebuild with the default `--no-strip-version`.
- **Zero probes for a multi-isoform gene** — should not happen. If it does,
  inspect `<target>_offtarget_counts.csv` for a homolog and `--allow` it after
  verifying.
