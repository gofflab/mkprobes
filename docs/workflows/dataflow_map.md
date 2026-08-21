# Dataflow map

A compact map of what goes into each step and what comes out. Use it to work
out where a missing file should have come from.

The steps themselves are in {doc}`../getting_started`.

## Step by step

| Step | Command(s) | Takes | Produces |
| --- | --- | --- | --- |
| 0. project | `init` | a directory name | `genes.txt`, `manifest.json`, `README.md` |
| 1. dataset | `prepare`, `ingest`, or `create-dataset` | species choice, or genome + GTF, or a FASTA | indexed dataset directory (bowtie2 index, `.jf` k-mer files, annotation) |
| 2. targets | `chkgenes`, then `convert-to-transcripts` | target list + dataset | `genes.converted.txt`, then `genes.converted.tss.txt` |
| 3. codebook | `make-codebook` | transcript list (+ optional expression table) | `codebook.json` + a logged hash |
| 4. probes | `run-panel` (or `candidates`, `screen`, `construct`) | dataset + codebook | per-target parquet chain, ending in `_final_` |
| 5. panel QC | `filter-genes` | output directory + target list | `genes.pass.txt`, warnings per thin target |
| 6. assembly | `check-manifest`, then `assemble short` / `assemble gen` | `manifest.json` + `_final_` parquet files | `generated/`: pool, FASTAs, provenance |

## Filenames chain from their inputs

This is the part that trips people up: most commands name their output after
their input rather than using a fixed name.

```text
genes.txt
  --chkgenes-->                genes.converted.txt
  --convert-to-transcripts-->  genes.converted.tss.txt
  --make-codebook-->           genes.converted.tss.codebook.json   (unless you pass -o)
```

Pass `-o codebook.json` to `make-codebook` to break that chain, which is what
{doc}`../getting_started` and the `mkprobes init` template both do.

Per target `T`, inside `output/`:

```text
T_all.parquet          every candidate position
T_bowtie.parquet       raw alignments
T_crawled.parquet      candidates + off-target context   (+ T_crawled.stats.json)
T_screened_ol<N>.parquet   selected probe pairs          (+ .stats.json)
T_final_BamHIKpnI_<bits>.parquet   encoded constructs
```

`_final_` is the one that counts as a result: `filter-genes` counts it, and
`assemble` reads it. The `<bits>` in the name come from the codebook, and
`BamHIKpnI` is fixed by the assay chemistry.

Decoding these filenames in full — including why the `ol` number is a gap and
not an overlap — is in {doc}`../reference/columns`.

Panel-level:

```text
codebook.json               from make-codebook
codebook.failed.txt         targets run-panel could not finish
codebook.acceptable.json    off-targets you accepted during `assemble short`
genes.pass.txt              from filter-genes
generated/<name>_final.txt  the orderable pool
```

## Checkpoints

1. **After step 1** — the dataset's index and `.jf` files exist and are readable.
2. **After step 2** — `genes.converted.tss.txt` has one line per target you expected.
3. **After step 3** — the codebook covers exactly your target list; hash recorded.
4. **During step 4** — a `_final_` parquet appears per target; `codebook.failed.txt` is empty.
5. **After step 5** — the pass list is as long as your panel.
6. **After step 6** — `generated/` holds the pool and `<panel>.provenance.json`.

## If a file is missing

| Missing | Look at |
| --- | --- |
| `*_crawled.parquet` | the `candidates` step — dataset path and target name |
| `*_screened_ol*.parquet` | `screen` constraints: `--minimum`, `--maxoverlap` |
| `*_final_*.parquet` | the target is in the codebook, and screened input exists |
| many targets in `codebook.failed.txt` | `run-panel --list-failed-all` for the off-target picture |
| everything, after a re-run | outputs are skipped when present; use `--overwrite` |

Any output parquet will tell you how it was made: `mkprobes provenance <file>`.
