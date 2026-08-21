# Troubleshooting

## First: see the actual error

Commands report failures as a single actionable line rather than a Python
traceback, because the traceback is rarely the useful part. When you need the
full one, put `--debug` **before** the command name:

```bash
mkprobes --debug run-panel data/mouse panel_a/codebook.json
```

It is an option on `mkprobes` itself, so `mkprobes run-panel --debug` will not
work. The traceback also always reaches the log file, whether or not you pass
the flag.

## Installation and environment

- **`mkprobes: command not found`** — activate your environment first:
  `source .venv/bin/activate` for uv, `conda activate mkprobes` for conda. If
  it is stale, re-run `uv sync` or `mamba env update -f environment.yml`.
- **`ModuleNotFoundError`** — same cause, same fix.
- **External tools fail** — confirm `bowtie2`, `jellyfish` and `gffread` are on
  `PATH`, on compute nodes as well as the login node. All three are on
  bioconda. `RepeatMasker` is optional and only needed at final assembly.

## Datasets

- **`is not a recognized probe dataset`** — the path you gave is not a dataset
  directory. Note that `mkprobes prepare data --species mouse` creates
  `data/mouse`, so later commands take `data/mouse`, not `data`.
- **Empty dataset after `ingest`** — almost always contig names that differ
  between genome and annotation. Run `mkprobes ingest ... --validate-only` and
  read the `SEQNAME_MISMATCH` entry.
- **A directory holding both a reference build and a `dataset.json`** — this is
  rejected rather than guessed at, because the two are screened differently.
  See [Which kind of dataset a command loads](reference/file_formats.md#which-kind-of-dataset-a-command-loads).

## Targets and codebook

- **`Transcript not found in ensembl`** — you passed a transcript where a gene
  was expected, or you are on a custom dataset where the Ensembl-backed modes
  do not apply. Use the offline modes: `--longest` or `--all`.
- **`Could not resolve 'X'`** — the annotation does not know that symbol. Use
  the dataset's own IDs, or register an ortholog table when building the
  dataset.
- **`lists N target(s) more than once`** — remove the repeats. A duplicate
  would take a second set of readout bits and corrupt the codebook, so it is
  refused rather than warned about.
- **Codebook landed somewhere unexpected** — `make-codebook` names its output
  after its input unless you pass `-o`. Pass `-o codebook.json`.

## Probe design

- **Very low probe counts** — inspect the target's `_crawled.stats.json`, or
  run `mkprobes run-panel ... --list-failed-all` for the off-target picture.
  One dominant cross-reactive binder is a different problem from diffuse loss.
- **Outputs skipped unexpectedly** — finished targets are skipped by design.
  Re-run with `--overwrite`, or name a single gene as the third argument to
  `run-panel` to force just that one.
- **`--restriction` rejected** — SOLAR chemistry fixes the pair to BamHI +
  KpnI. Drop the option; the default is correct. See
  {doc}`workflows/design_probes`.
- **Some targets never finished** — read `codebook.failed.txt` beside the
  codebook, and each target's `output/<gene>.log`.

## Assembly

- **Manifest rejected** — run `mkprobes check-manifest manifest.json`; it names
  the field and the fix. `codebook` paths resolve relative to the manifest
  file, not to your shell's working directory.
- **`bcidx` out of range** — each index consumes two rows of the header/footer
  table; the error names the maximum.
- **`RepeatMasker` not found** — install it, or pass `--skip-repeatmasker`. For
  non-model species use `--rm-species <taxon>`.
- **Missing `_final_` parquet files** — `run-panel` did not complete for every
  target in the codebook.

## Working out what produced a file

```bash
mkprobes provenance panel_a/output/Sox2_final_BamHIKpnI_1,2,3.parquet
```

This prints the mkprobes version, timestamp, command line, dataset and
parameters embedded in any output parquet. Files written before provenance was
recorded, and files not written by `mkprobes`, have none — the command says so
rather than guessing.
