# Choose your targets

**Step 2 of the {doc}`../getting_started` workflow.**

Two commands turn the list of genes you care about into the transcript-level
target list everything downstream uses. Doing this before the codebook exists
is deliberate: a name that turns out to be wrong is free to fix now and
expensive to fix after probes have been designed against it.

## Writing the target list

One target per line. Blank lines are ignored, everything after a `#` is a
comment, and inline comments work, so the file can carry its own rationale:

```text
# Panel A - dorsal telencephalon, v2
Sox2
Pax6      # dorsal telencephalon marker
Eomes
# Gad1    - dropped, too low in this region
```

For a reference dataset, use gene names. For a custom dataset, use whatever
IDs its annotation actually carries — gene IDs, transcript IDs, or symbols
resolvable through a registered ortholog table.

Naming a target twice is a hard error, not a warning: the duplicate would
claim a second set of readout bits and silently corrupt the codebook. Order is
preserved, because it feeds bit assignment.

## 1. Check the names

```bash
mkprobes chkgenes data/mouse panel_a/genes.txt
```

This resolves every name against the dataset and writes
`panel_a/genes.converted.txt` — the same list, normalized. On reference
datasets it checks against Ensembl and writes `genes.mapping.json` recording
any renames it applied. On custom datasets the check is entirely offline,
against the dataset's own annotation plus any registered alias or ortholog
tables.

Unresolvable names fail loudly, with close-match suggestions. Fix them in
`genes.txt` and re-run before continuing.

## 2. Resolve to transcripts

```bash
mkprobes convert-to-transcripts data/mouse panel_a/genes.converted.txt
```

This picks one transcript per gene and writes
`panel_a/genes.converted.tss.txt`. That file is the target list every later
step takes.

Both commands name their output after their input, so the chain is:

```text
genes.txt
  --chkgenes-->               genes.converted.txt
  --convert-to-transcripts--> genes.converted.tss.txt
```

If you skip `chkgenes` and run `convert-to-transcripts` on `genes.txt`
directly you get `genes.tss.txt` instead. Either is fine, as long as you pass
the right filename to the next step.

### Picking the selection mode

`-m/--mode` sets the policy. Fix it before you design probes and do not change
it mid-panel; switching modes changes which sequence probes were designed
against.

| Mode | What it picks | Where it works |
| --- | --- | --- |
| `canonical` | the canonical transcript (default) | reference; falls back to `longest` on custom |
| `longest` | per gene, the isoform with the longest sequence | anywhere, offline |
| `all` | every isoform | anywhere, offline |
| `gencode`, `ensembl`, `appris`, `apprisalt` | annotation-specific sets | reference only (human/mouse) |

On a custom dataset, pass `-m longest` explicitly. `canonical` will fall back
to it anyway, and being explicit keeps the command self-documenting.

`longest` measures the **sequence**, read from the FASTA, so introns do not
distort the choice.

## Looking up one transcript

`mkprobes transcripts` answers the same question for a single gene without
writing a list — useful for spot-checks:

```bash
mkprobes transcripts data/mouse --gene Sox2 --canonical
mkprobes transcripts data/myspecies --gene Och.576 --longest
```

Here the modes are flags (`--longest`) rather than `-m longest`.

## When it goes wrong

- **`Transcript not found in ensembl`** — the name is a transcript where a
  gene was expected, or the dataset is custom and the reference-only modes do
  not apply. Use `--longest` / `--all`.
- **`Could not resolve 'X'`** — the annotation does not know that symbol. Use
  the dataset's own IDs, or register an ortholog table when you build the
  dataset (`--annotation-table orthologs=...`).
- **`lists N target(s) more than once`** — the target list repeats a name.
  Every command that reads a target list rejects this before doing any work;
  remove the repeats and re-run.

---

Next: {doc}`design_the_codebook`.
