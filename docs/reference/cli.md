# CLI reference

Every command, every flag, generated directly from the code at build time.
Nothing on this page can drift out of date; if a flag is here, it exists.

This is a lookup table, not a tutorial. If you are designing a panel for the
first time, read {doc}`../getting_started` instead — it puts these commands in
order and explains what each one is for.

## The order to run them in

| Step | Command | How-to |
| --- | --- | --- |
| 0. project | `init` | {doc}`../getting_started` |
| 1. dataset | `prepare` (mouse/human), `ingest` or `create-dataset` (any species) | {doc}`../workflows/build_a_dataset` |
| 2. targets | `chkgenes`, then `convert-to-transcripts` | {doc}`../workflows/choose_your_targets` |
| 3. codebook | `make-codebook` | {doc}`../workflows/design_the_codebook` |
| 4. probes | `run-panel` (wraps `candidates`, `screen`, `construct`) | {doc}`../workflows/design_probes` |
| 5. panel QC | `filter-genes` | {doc}`../workflows/qc_your_panel` |
| 6. assembly | `check-manifest`, then `assemble short` / `assemble gen` | {doc}`../workflows/order_your_oligos` |

Supporting commands, usable at any point: `provenance` (how was this file
made?), `hash` (codebook identity), `transcripts` (one-off transcript lookup).

## Things worth knowing before you read the list

**`--debug` goes before the command name.** It is an option on `mkprobes`
itself, not on individual commands, so it is `mkprobes --debug run-panel ...`,
never `mkprobes run-panel --debug`. Without it, a failure is reported as one
actionable line; with it, you get the full Python traceback.

**`--restriction` is not a free choice.** It appears on `screen`, `construct`
and `run-panel`, but SOLAR chemistry fixes the pair to **BamHI + KpnI**. The
header/footer sequences carry those two sites and final assembly excises the
probe with a KpnI/BamHI double digest, so any other pair yields probes that
nothing downstream can cut out. Anything else is refused up front, with an
explanation. In practice: leave the option alone.

The spelling does differ between commands, which is worth knowing when you
copy a command line around:

- `screen` and `run-panel` take one comma-separated value: `--restriction BamHI,KpnI`
- `construct` takes a repeatable option: `--restriction BamHI --restriction KpnI`

**`assemble` takes its manifest before the subcommand.** The manifest is an
argument on the group, because `gen` and `short` share it:

```bash
mkprobes assemble panel_a/manifest.json gen      # correct
mkprobes assemble gen panel_a/manifest.json      # wrong
```

The same applies to `--headerfooter`, which belongs to the `assemble` group
rather than to `gen`. A consequence: `mkprobes assemble ... gen --help` still
has to parse and validate the manifest first, so it needs a real one.

**Which dataset a command loads is inferred from the directory.** Reference
datasets (`prepare`) and custom ones (`ingest` / `create-dataset`) are screened
differently, and the difference is not cosmetic. See
[Which kind of dataset a command loads](file_formats.md#which-kind-of-dataset-a-command-loads).

**Target lists tolerate comments.** Anywhere a command takes a `GENES` file,
blank lines are skipped and everything after a `#` is a comment, so you can
record why a target is in the panel. Listing a target twice is an error.

**Every parquet output records how it was made.** `mkprobes provenance <file>`
prints the version, timestamp, command line, dataset and parameters embedded
in any output parquet. The `.stats.json` sidecars carry the same record under
a `provenance` key, and `assemble` writes a `<panel>.provenance.json` beside
the pool.

## Commands

```{eval-rst}
.. click:: mkprobes.cli:main
   :prog: mkprobes
   :nested: full
```

## Deprecated script shims

The files under `scripts/probegen/` predate the package and are kept only as
shims:

| Script | Replaced by |
| --- | --- |
| `o_codebook.py` | `mkprobes make-codebook` |
| `1_run_codebook*.py` | `mkprobes run-panel` |
| `2_assemble_manifest.py` | `mkprobes assemble` |

The rest of that directory is exploratory notebooks: `simulate.py` (in-silico
validation) and `foridt.py` / `adt.py` (IDT ordering examples whose logic
already lives in the package).
