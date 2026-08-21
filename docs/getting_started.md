# Getting started

This is the walkthrough: the whole `mkprobes` workflow, start to finish, on a
mouse panel. Every other workflow page assumes you arrived from here.

Design is six steps, each one command, after a one-command setup:

```text
0. project    mkprobes init          scaffold genes.txt + manifest.json

1. dataset    mkprobes prepare       (mouse/human) or mkprobes ingest (any species)
2. targets    mkprobes chkgenes  ->  mkprobes convert-to-transcripts
3. codebook   mkprobes make-codebook
4. probes     mkprobes run-panel     candidates -> screen -> construct, all targets
5. panel QC   mkprobes filter-genes
6. assembly   mkprobes assemble      -> orderable oligos
```

This is the same six-step order `mkprobes --help` prints, and the same order
`mkprobes init` writes into your project's `README.md`.

Before you begin: installation is complete ({doc}`installation`), `mkprobes --help`
works, and `bowtie2`, `jellyfish` and `gffread` are on your `PATH`. If any of
that is unfamiliar, read {doc}`before_you_start` first — it covers what to
download, how much disk and time each step needs, and what the external
programs are for. New to the assay itself? {doc}`what_is_solar` explains what
these probes physically are before you design any.

If a term here is unfamiliar, it is almost certainly in the {doc}`glossary`.

## 0. Create the project

```bash
mkprobes init panel_a --species mouse --dataset ../data/mouse
```

This writes a directory that is already valid:

| File | What it is |
| --- | --- |
| `genes.txt` | Your targets, one per line. Edit this first. |
| `manifest.json` | Describes the panel for the assembly step, with every field commented. |
| `README.md` | The same commands as below, filled in with your paths. |

Start here rather than hand-writing `manifest.json`. The manifest has two
fields that are easy to get wrong — `bcidx` indexes rows of an internal
table, and `n_probes` accepts two magic strings — and `init` fills both with
working values. Pass `--force` to overwrite an existing project.

Now edit `panel_a/genes.txt`. One target per line; blank lines are ignored,
and everything after a `#` is a comment, so you can record why each target is
in the panel:

```text
Sox2
Pax6      # dorsal telencephalon marker
# Gad1    - dropped, too low in this region
```

Naming a target twice is an error, not a warning: a duplicate would claim a
second set of readout bits and quietly corrupt the codebook.

## 1. Build the dataset (once per species)

A dataset is the sequence universe probes are searched in and screened
against. Build it once and treat it as shared, immutable input.

Mouse or human — downloaded and indexed for you:

```bash
mkprobes prepare data --species mouse --threads 16
```

This writes `data/mouse/`. Note the species subdirectory: later commands take
`data/mouse`, not `data`.

Any other species, from a genome FASTA plus an annotation:

```bash
mkprobes ingest data/myspecies --genome genome.fa.gz --gtf annotation.gtf --species myspecies
```

Run `ingest` with `--validate-only` first and read the report — it catches the
mismatches that otherwise produce a silently empty dataset.

Details: {doc}`workflows/build_a_dataset`. Non-model species from end to end:
{doc}`workflows/solar_new_species`.

## 2. Resolve your targets

Two commands: check that the names exist, then pick one transcript per gene.

:::{tip}
Have expression data and room left in the panel? `mkprobes suggest-targets` can
propose genes that add information your current list does not already carry, and
score how much of the biology the panel captures. Optional, and it reads your
expression data rather than the dataset — see
{doc}`workflows/choose_your_targets`.
:::

```bash
mkprobes chkgenes data/mouse panel_a/genes.txt
mkprobes convert-to-transcripts data/mouse panel_a/genes.converted.txt
```

Each command names its output after its input, so the files chain:

```text
genes.txt  --chkgenes-->  genes.converted.txt  --convert-to-transcripts-->  genes.converted.tss.txt
```

`genes.converted.tss.txt` is the target list every later step uses. On a
custom dataset add `-m longest` to `convert-to-transcripts`; the default
`canonical` mode needs Ensembl and falls back to `longest` anyway.

Details: {doc}`workflows/choose_your_targets`.

## 3. Generate the codebook

```bash
mkprobes make-codebook data/mouse panel_a/genes.converted.tss.txt -o panel_a/codebook.json
```

This sizes the code from your target count, assigns three readout bits per
target, fills spare capacity with `Blank-N` decoys, and logs a hash of the
result for provenance.

Pass `-o` explicitly, as above. Without it the output name is derived from the
input list — `genes.converted.tss.txt` becomes `genes.converted.tss.codebook.json`,
which is not what most people expect.

Have per-target expression data? `--expression <table>` balances signal across
readout bits. It is optional and the panel works without it.

Details: {doc}`workflows/design_the_codebook`.

## 4. Design probes for every target

```bash
mkprobes run-panel data/mouse panel_a/codebook.json
```

This is the long step — hours for a real panel. The codebook is the work list:
every target in it runs `candidates -> screen -> construct` in parallel (16
workers by default).

It is safe to interrupt and re-run. Finished targets are skipped, failures are
isolated per gene rather than killing the run, and the names of failed genes
are collected in `codebook.failed.txt`. Re-run one target with
`mkprobes run-panel data/mouse panel_a/codebook.json Sox2`.

Details, and the single-target commands underneath it:
{doc}`workflows/design_probes`.

## 5. Check the panel

```bash
mkprobes filter-genes panel_a/output --genes panel_a/genes.converted.tss.txt \
    --min-probes 48 --out panel_a/genes.pass.txt
```

This counts what each target would actually contribute to an oligo order and
warns about each thin one individually. Targets that never produced probes at
all are reported separately — those are `run-panel` failures, not thin panels.

Details: {doc}`workflows/qc_your_panel`.

## 6. Assemble the orderable pool

Triage the thin targets first, then build:

```bash
mkprobes check-manifest panel_a/manifest.json
mkprobes assemble panel_a/manifest.json short 12
mkprobes assemble panel_a/manifest.json gen
```

`short` walks you through low-count targets and records any off-targets you
decide to accept in `codebook.acceptable.json`, which `run-panel` picks up on
the next run. `gen` writes the pool. Assembly is deterministic: the same
inputs always produce the same oligos.

Note the argument order — the manifest comes *before* `gen` or `short`,
because both share it.

Details: {doc}`workflows/order_your_oligos`.

## What you end up with

Per target, in `panel_a/output/`:

```text
*_all.parquet -> *_bowtie.parquet -> *_crawled.parquet
    -> *_screened_ol*.parquet -> *_final_BamHIKpnI_<bits>.parquet
```

Per panel, in `panel_a/generated/`: the parquet, the splint and padlock
FASTAs, `<panel>_final.txt` (the orderable pool, one oligo per line), and
`<panel>.provenance.json`.

Every parquet also carries its own provenance internally. To find out how any
output file was made:

```bash
mkprobes provenance panel_a/output/Sox2_final_BamHIKpnI_1,2,3.parquet
```

What the columns in those files mean: {doc}`reference/columns`.

## When something goes wrong

Commands report failures as a single actionable line rather than a Python
traceback. To see the traceback, put `--debug` before the command name:

```bash
mkprobes --debug run-panel data/mouse panel_a/codebook.json
```

Common problems and their fixes: {doc}`troubleshooting`.

## Next steps

- Non-model species, end to end: {doc}`workflows/solar_new_species`
- How the assay and the algorithms work: {doc}`workflows/assay_and_under_the_hood`
- Input/output flow map: {doc}`workflows/dataflow_map`
- Every command and flag: {doc}`reference/cli`
