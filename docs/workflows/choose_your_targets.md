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
resolvable through a registered annotation table.

### Writing targets as gene names on a custom dataset

De novo annotations usually have no gene names at all: `mkprobes ingest`
reports `GENE_NAME_FALLBACK` and names fall back to IDs like `Och.958.1`,
which nobody wants to write a panel in. If you have a table mapping those IDs
to names — an ortholog assignment, a curated symbol list — register it and say
which column holds the names you want to use:

```bash
mkprobes ingest data/myspecies --genome genome.fa --gtf annotation.gtf \
    --species myspecies \
    --annotation-table annot=my_annotation.tsv \
    --gene-name-column Hsapiens_gene_name
```

Targets can then be written as those names, and `chkgenes` resolves them to
transcript IDs. Two details worth knowing:

- **Only that column is searched.** Without `--gene-name-column`, lookup scans
  every text column of every registered table. That is fine on a small table
  and slow on a wide one — and it will happily match a protein sequence or an
  embedding column if something in it looks like your gene name.
- **Comma-separated cells count as one name each.** Ortholog tables routinely
  map one transcript to several symbols (`UBE2A,UBE2B`), so each entry is
  matched separately rather than requiring the whole cell to match.

The column must exist in one of the registered tables; ingest checks at build
time and lists the available columns if it does not.

Naming a target twice is a hard error, not a warning: the duplicate would
claim a second set of readout bits and silently corrupt the codebook. Order is
preserved, because it feeds bit assignment.

## Optional: let expression data fill out the panel

Skip this section if you already know every gene you want.

A panel has room for a fixed number of genes, and a hand-picked list tends to
cluster: markers chosen for the same cell type report much the same thing, so
the panel measures one axis of the biology several times over and misses
others. If you have expression data for the tissue — your own, or a published
atlas — `suggest-targets` proposes genes that carry information your current
choices do not.

```bash
mkprobes suggest-targets atlas.h5ad --add 40 --have genes.txt -o genes.txt
```

It regresses every candidate gene against the ones you already hold, then picks
the genes that best span whatever variation is left over. It also reports how
much of the data's structure the panel captures, before and after, so the
suggestion comes with a number rather than only a list:

```text
Panel of 43 captures 38.5% of the variance in the top 30 PCs
  (up from 6.8% with your 3 alone).
```

The output is an ordinary target list with your genes first, in order, followed
by the suggestions. **Read it and edit it.** These are candidates ranked by a
statistical criterion, not by whether they make biological sense, are expressed
highly enough to detect, or matter to your question.

:::{important}
Filter your expression data to informative genes first — scanpy's
`highly_variable_genes` is the usual route.

A gene that correlates with nothing looks maximally independent to this method,
so unfiltered data makes it prefer genes that are merely noisy over genes that
report real biology. The effect is not subtle: on test data where 8 latent
programmes were present, adding a few hundred unstructured genes changed the
result from covering most programmes to picking nothing but noise.

The command checks for this and warns when the suggestions look like noise, but
the check is a backstop, not a substitute for filtering.
:::

Two options are worth knowing. `--layer` chooses which expression matrix to use,
if your file has more than one — normalized, log-transformed values usually work
better than raw counts. `--n-components` sets how many dimensions of the
leftover variation to select against; it defaults to `--add + 20` (minimum 50),
and lowering it helps when your data has only a few distinct programmes. The
choice materially changes which genes come out, so it is worth trying a couple
of values and comparing the variance-capture numbers.

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
