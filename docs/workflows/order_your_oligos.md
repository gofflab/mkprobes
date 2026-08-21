# Order your oligos

**Step 6 of the {doc}`../getting_started` workflow.**

The last step turns per-target designs into the actual list of DNA sequences
you send to a vendor. Output is deterministic: the same inputs always produce
the same pool.

## The manifest

Assembly is driven by a `manifest.json` describing the panel. `mkprobes init`
writes a valid one with every field commented, which is the easiest way to get
a correct one:

```json
[
  {
    "name": "panel_a",
    "species": "mouse",
    "codebook": "codebook.json",
    "bcidx": 0,
    "n_probes": 24
  }
]
```

| Field | What it means |
| --- | --- |
| `name` | names the output files in `generated/` |
| `species` | picks the RepeatMasker taxon; any name is accepted |
| `codebook` | path to the codebook, **relative to the manifest** |
| `bcidx` | which header/footer pair to build against |
| `n_probes` | maximum probes per target: a number, or `"high"` (34) or `"low"` (16). Omit to let the species decide |

`bcidx` is the field people get wrong. Each index consumes two rows of an
internal header/footer table, so the valid range is bounded — use a different
index for each panel you intend to pool together, and `check-manifest` will
tell you the maximum if you exceed it.

## Check it first

```bash
mkprobes check-manifest panel_a/manifest.json
```

This validates the whole manifest — schema, `bcidx` range, that each named
codebook actually exists, that panel names are unique — and prints a summary
line per probe set. It takes a second. Run it before `gen`, which otherwise
spends hours proving the same point.

(`assemble` runs the same validation itself, so a bad manifest fails fast
either way. `check-manifest` just lets you check without starting anything.)

One consequence: because it checks that the codebook exists, `check-manifest`
will fail on a project straight out of `mkprobes init` — the codebook does not
exist until step 3. That is correct behaviour, not a broken manifest.

## Triage the thin targets

```bash
mkprobes assemble panel_a/manifest.json short 12
```

This finds targets with fewer than 12 probes and, where off-target tables
exist, walks you through their binders interactively so you can accept the
ones you do not mind labelling. Accepted off-targets are written to
`<codebook>.acceptable.json`.

That file is not applied retroactively — `mkprobes run-panel` picks it up on
its **next** run and re-designs the affected targets. So the loop is:

```text
assemble ... short N   ->   run-panel (re-designs affected targets)   ->   assemble ... gen
```

## Build the pool

```bash
mkprobes assemble panel_a/manifest.json gen
```

For non-model species, RepeatMasker has no built-in mapping, so either give it
a taxon it knows or skip it explicitly:

```bash
mkprobes assemble panel_a/manifest.json gen --rm-species mollusca
mkprobes assemble panel_a/manifest.json gen --skip-repeatmasker
```

RepeatMasker is an optional dependency, used only at this step.

## Note the argument order

The manifest comes **before** the subcommand, because `gen` and `short` share
it:

```bash
mkprobes assemble panel_a/manifest.json gen      # correct
mkprobes assemble gen panel_a/manifest.json      # wrong
```

The same is true of `--headerfooter`, which overrides the vendored
splint/padlock header/footer table. It belongs to `assemble`, not to `gen`:

```bash
mkprobes assemble --headerfooter custom.csv panel_a/manifest.json gen
```

You will almost never need it; the table that ships with the package encodes
the assay chemistry.

## What assembly does

For each probe set in the manifest:

1. Loads the codebook and the per-target `_final_` parquet files.
2. Selects the top probes per target (`n_probes`, or the species default).
3. Optionally runs RepeatMasker.
4. Builds the final splint and padlock constructs, stitching on header and
   footer.
5. Asserts the geometry — every pair is checked for splint/padlock structure
   and 139–150 nt padlock length. A failure here is a bug, not a tuning knob.
6. Writes the outputs and a provenance record.

## What you get

Under `panel_a/generated/`:

| File | What it is |
| --- | --- |
| `<name>_final.txt` | **the orderable pool** — one oligo per line |
| `<name>.parquet` | the same content with all columns |
| `<name>_pad.fasta` | padlock sequences |
| `<name>_splint.fasta` | splint sequences |
| `<name>.provenance.json` | timestamp, version, codebook hash, probe set config, RepeatMasker status |
| `_allout<timestamp>.txt` | cumulative pool across runs |

`<name>_final.txt` is the file you send to the vendor. What you are physically
ordering, and why there are two oligos per site, is explained in
{doc}`../what_is_solar`.

## On a cluster

- Run `short` first as a cheap gate.
- Run `gen` only after every acceptable-off-target decision is frozen.
- Archive the manifest, codebook, `.acceptable.json` and everything in
  `generated/` together as one release bundle.

## When it goes wrong

- **Missing `_final_` parquet files** — `run-panel` did not finish for every
  target in the codebook. Check `codebook.failed.txt`.
- **RepeatMasker not found** — install it, or pass `--skip-repeatmasker`.
- **Manifest path problems** — `codebook` is resolved relative to the manifest
  file, not to your shell's working directory. `check-manifest` reports this
  precisely.
- **`bcidx` out of range** — the error names the maximum. Each index uses two
  rows of the header/footer table.

---

That is the whole workflow. To confirm how any output was produced, at any
point: `mkprobes provenance <file>`.
