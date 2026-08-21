# Before you start

The commands in this documentation assume three things are already true: the
right helper programs are installed, the right reference files are downloaded,
and the machine you are running on is big enough. None of those are obvious if
you have not done this before, and all three are places where a first attempt
typically stalls. This page covers them.

Software installation itself is in {doc}`installation`. This page is about
what that software needs to work with, and what it will cost you.

## The four external programs

`mkprobes` is a Python package, but it does not do the heavy sequence work
itself. It shells out to four established bioinformatics programs, which must
be on your `PATH` — meaning your shell can find them by name when you type
them. Three of the four are installed for you if you use the conda/mamba
environment described in {doc}`installation`.

`bowtie2` — the aligner
: Answers the question "where else in this transcriptome does this 20-base
  sequence stick?" for every candidate probe. This is the specificity check
  and the slowest step in probe design. It is used twice: `bowtie2-build`
  makes a searchable index of your transcriptome once when the dataset is
  built, and `bowtie2` itself queries that index every time a gene is designed.

  **It must be `bowtie2`, not `bowtie`.** These are two different programs.
  The original `bowtie` (sometimes called bowtie1) uses a different, mutually
  incompatible index format and does not support the local-alignment mode this
  pipeline relies on. If a page elsewhere says "bowtie", read it as `bowtie2`.
  Check with `bowtie2 --version`, which should print a version starting with 2.

`jellyfish` — the k-mer counter
: Builds fast lookup tables of every 18-letter (and, for the blocklist,
  15-letter) substring in a set of sequences, so the pipeline can ask "does
  this probe share an exact stretch with anything I want to avoid?" without
  running an alignment. Used when a dataset is built; not used during probe
  design itself. On bioconda the package is called `kmer-jellyfish` but
  installs a binary named `jellyfish`.

`gffread` — the annotation-to-sequence extractor
: Takes your genome FASTA plus your GTF/GFF3 annotation and writes out the
  actual spliced transcript sequences — that is, it looks up each transcript's
  exon coordinates, pulls those stretches out of the genome, and glues them
  together. It also converts GFF3 to GTF when needed. Required only when
  ingesting a new species from a genome; not needed if you already have a
  transcriptome FASTA, and not needed during probe design.

`RepeatMasker` — optional, final assembly only
: Masks repetitive sequence in the finished splint and padlock oligos, after
  which probes containing masked bases are dropped. It is invoked **only** at
  the final assembly step, and only when a RepeatMasker taxon is available for
  your species.

  RepeatMasker is **commented out** in the repository's `environment.yml`, so
  a standard conda install does not give it to you. That is deliberate: it
  pulls in large repeat libraries and most of the workflow does not need it.
  Install it separately when you get to assembly, or skip it explicitly. For a
  non-model species with no RepeatMasker library entry, skipping is the normal
  outcome and the assembly step will tell you so.

Which step needs which:

| Step | bowtie2 | jellyfish | gffread | RepeatMasker |
| --- | :---: | :---: | :---: | :---: |
| Building a reference dataset (human/mouse download) | yes | yes | – | – |
| Ingesting a new species (genome + annotation) | yes | yes | yes | – |
| Generating a codebook | – | – | – | – |
| Designing probes | yes | – | – | – |
| Panel QC | – | – | – | – |
| Final assembly | – | – | – | optional |

## Reference files you need to obtain

For human and mouse, you do not need to download anything by hand — the
reference-preparation step fetches everything itself: a GENCODE annotation, an
Ensembl annotation, cDNA and ncRNA sequence files, a tRNA archive, and an
APPRIS principal-isoform table.

For **any other species**, you supply the files. There are two required and
two strongly recommended.

### Required: a genome FASTA and a matching annotation

**Genome FASTA.** The DNA sequence of the assembly, one record per chromosome
or scaffold. May be gzipped; it is decompressed to a temporary file during
ingestion because `gffread` needs random access to it.

**Annotation, in GTF or GFF3 format.** This is the file that says which parts
of the genome are exons of which transcript of which gene. Without it, a
genome is just letters. A GTF line is one feature — a gene, a transcript, an
exon — with a sequence name, start and end coordinates, a strand, and a final
column of attributes such as `gene_id "..."; transcript_id "...";`. GFF3 is
the same idea with `key=value` attributes; `mkprobes` detects which one you
gave it and converts GFF3 to GTF automatically.

The single most important property is that **the annotation must have been
built for the exact assembly of the genome FASTA you downloaded**. If it was
not, the sequence names will not match — the classic symptom being an
annotation that says `chr1` against a genome whose records are named `1`, or a
scaffold-level annotation against a chromosome-level assembly. The result is
not an error; it is silently empty output. The ingestion step checks this for
you and reports it as `SEQNAME_MISMATCH`, which is why it is worth running the
validation-only mode first and reading the report before committing to a full
build.

Two other things the validator will tell you, both worth knowing in advance:

- **Every transcript needs exon rows.** Transcripts without them are dropped,
  because `gffread` has nothing to extract.
- **Gene names are optional.** Many de novo annotations (StringTie, AUGUSTUS,
  MAKER) have no gene-name attribute at all, only IDs. That is fine, but it
  means you must select targets by ID rather than by a familiar gene symbol —
  the report flags this as `GENE_NAME_FALLBACK`. If you have an ortholog
  table mapping your IDs to, say, human symbols, you can register it with the
  dataset and select targets through it.

### Strongly recommended: rRNA and tRNA sequences for the blocklist

Supply FASTA files of your species' ribosomal RNA and transfer RNA sequences.
Their 15-letter subsequences are counted into a **blocklist**, and any
candidate probe that shares even one 15-mer with that set is thrown out before
it is ever considered.

This matters more than it sounds. Ribosomal RNA is not a minor contaminant of
a cell — it is the overwhelming majority of the RNA present, often well over
80% of total RNA. A probe that happens to have a short exact stretch in common
with rRNA will find an enormous number of targets in every cell, in every
compartment, and it will do so brightly. The failure mode is not subtle
background: it is a probe that lights up the whole tissue and quietly ruins
whichever bit it was assigned to. Because rRNA and tRNA are also highly
structured and highly conserved, the risk does not go away just because your
species is unusual.

If you build a dataset without a blocklist, the pipeline does not fail — it
logs `No rRNA/tRNA blocklist k-mers in this dataset; blocklist filtering is a
no-op` on every single run. If you see that line, you are designing without
this protection.

A partial alternative, if your GTF has a biotype column, is to let the
ingestion step pull rRNA/tRNA/snoRNA transcripts out of your own annotation
and build the blocklist from those. Many de novo annotations have no biotype
column, in which case you need external FASTAs.

### Where these files come from

Without inventing specific URLs, the usual sources are:

- **Genome assemblies and annotations**: NCBI (Genome / RefSeq / Datasets) and
  Ensembl, including Ensembl's clade-specific sites for non-vertebrates. For
  many non-model organisms the assembly comes from a consortium or a single
  publication's data repository rather than from either.
- **Your own de novo annotation**: perfectly acceptable, and common for
  non-model species. StringTie, MAKER, AUGUSTUS/BRAKER and NCBI Gnomon outputs
  all work; the ingestion validator is written with them in mind.
- **rRNA sequences**: SILVA is the standard curated rRNA database. Failing
  that, rRNA genes from your own annotation, or rRNA from a close relative,
  are better than nothing.
- **tRNA sequences**: GtRNAdb hosts tRNA gene sets for a wide range of
  genomes. tRNAscan-SE run over your own assembly is the fallback when your
  species is not in it.

Whatever you download, keep a note of the release or version. The provenance
manifest written at ingestion records a checksum of each input file and the
exact command you ran, but the fields for assembly name, source, and release
date are stubs for you to fill in. Fill them in. Six months later that file is
the only record of what your panel was designed against.

## What it costs: disk, memory, and time

The figures below are measured from this repository, on a laptop, using a
53,045-transcript *Octopus chierchiae* dataset built from a StringTie
annotation (168 MB of extracted transcript sequence). They are anchors, not
promises — the honest scaling rule is given alongside each one, and your
numbers will differ with transcriptome size, transcript length, and how many
cores you actually have.

| Step | Disk | Memory | Time |
| --- | --- | --- | --- |
| Reference download (human/mouse) | the six downloaded files, plus a concatenated transcriptome FASTA and its index | dominated by k-mer counting, see below | dominated by download speed |
| Ingest a new species | ≈ 3× the extracted transcript FASTA, plus temporary space for the decompressed genome | scales with transcriptome size | ≈ 9 min for 53k transcripts / 168 MB |
| Generate a codebook | kilobytes | negligible | seconds |
| Design probes | ≈ 1–2.5 MB per target | ≈ 0.5 GB **per parallel worker** | ≈ 5 s per target, per worker |
| Panel QC | none | negligible | seconds |
| Final assembly | a few MB per panel | modest, unless RepeatMasker runs | minutes, plus RepeatMasker if enabled |

### Disk

A built dataset is roughly **three times the size of the extracted transcript
FASTA**. For the octopus dataset: 168 MB of transcript sequence became 482 MB
on disk, made up of the FASTA itself, a 245 MB bowtie2 index (six `.bt2`
files), a 33 MB k-mer table, and a small annotation cache. Scale that ratio to
your own transcriptome.

Ingestion additionally needs **temporary free space equal to the uncompressed
genome**, because a gzipped genome is decompressed before `gffread` can use
it. That temporary copy is removed afterwards unless you asked to keep the
genome inside the dataset.

Design output is around **1–2.5 MB per target**, of which the great majority is
the `*_all.parquet` file holding every alignment of every candidate. A
300-gene panel therefore lands under a gigabyte, and you can reclaim most of
it by deleting the `_all` files once each gene is finished.

### Memory

Memory is the requirement most likely to bite you, and the reason is
parallelism rather than any single large object.

Probe design runs **16 worker processes by default**, one target per worker.
Each worker independently loads the dataset — including the full k-mer table
as an in-memory set. For the octopus dataset that k-mer table holds 1.5
million entries, and a freshly loaded dataset measures about **0.45 GB of
resident memory per process**. Sixteen of those is roughly 7 GB before any
alignment has happened. On top of that, each worker launches its own `bowtie2`,
which loads the transcriptome index (the `.1.bt2` and `.2.bt2` files —
about 100 MB for this dataset).

The practical rules:

- Peak memory during probe design ≈ (workers) × (dataset load + bowtie2 index).
  Estimate the per-worker cost from your own dataset's k-mer file and `.bt2`
  file sizes.
- **If you run out of memory, reduce the worker count.** That is the single
  effective lever. It costs wall-clock time and nothing else.
- The same argument applies to CPU. Each worker's `bowtie2` is invoked asking
  for 16 threads, so 16 workers request 256 threads between them. On a laptop
  or a small scheduler allocation, that oversubscription makes things slower,
  not faster. Match the worker count to the cores you actually have.

Dataset **building** has a separate memory peak: k-mer counting. When
preparing a human or mouse reference, `jellyfish` is asked for a 10-billion-slot
initial hash table (`-s 10G`), which is by far the largest single allocation
anywhere in the pipeline and sets the memory ceiling for that step. If any step
runs out of memory, expect it to be this one, and run it somewhere with plenty
of RAM rather than on a laptop. Custom datasets built by ingestion do **not**
use that fixed size — they size the hash from the FASTA (a quarter of its byte
size, with a 10-million floor), which is far smaller and adapts to your input.

### Time

Per-target design is fast: on the octopus dataset, one 16.6 kb transcript went
from start to finished construct in **5.4 seconds** — about 0.8 s to enumerate
candidates, 2.3 s for the `bowtie2` alignment (which produced 171,659
alignments), and roughly a second for screening and construction. Multiply by
your gene count and divide by your worker count for a first estimate, then
allow for the fact that longer transcripts and larger transcriptomes both cost
more alignment time.

Dataset building is the slow part and is done once per species. Ingesting the
octopus genome and annotation — validation, transcript extraction, bowtie2
index, k-mer counting — took about **9 minutes** end to end. A larger, more
fragmented transcriptome will take proportionally longer; the bowtie2 index
build is usually the dominant term.

Reference preparation for human or mouse is dominated by download time for six
files from three different servers, and is worth doing once into a shared
location rather than repeating per project.

Final assembly is minutes — unless RepeatMasker is enabled, in which case
RepeatMasker's runtime dominates everything else in the workflow.

## A sensible order of operations

1. Install, activate the environment, and confirm all three core tools answer
   `--version`. See {doc}`installation`.
2. Obtain your genome FASTA, its matching annotation, and rRNA/tRNA FASTAs if
   you can get them.
3. Run the dataset ingestion in **validation-only** mode first and read the
   report. Fix sequence-name mismatches and missing exon rows *before*
   building anything — a build on a mismatched pair produces an empty or
   silently truncated transcriptome.
4. Build the dataset, once, somewhere shared and backed up.
5. Design one or two genes first and look at the output before launching the
   whole panel. Check the alignment summary file for unexpected top binders,
   and check that the strict filter tiers are producing most of the probes.
6. Then run the panel.

## Where to go next

- What the assay is doing and why: {doc}`what_is_solar`
- The vocabulary, defined: {doc}`glossary`
- The commands, in order: {doc}`getting_started`
- The full new-species runbook: {doc}`workflows/solar_new_species`
- What the output files contain: {doc}`reference/columns`
