# Design probes

**Step 4 of the {doc}`../getting_started` workflow.**

This is where design intent becomes actual sequence, and it is the long step —
hours for a real panel. One command does the whole panel; the three commands
underneath it exist for debugging a single target.

## The one command

```bash
mkprobes run-panel data/mouse panel_a/codebook.json
```

The codebook is the work list: every target in it runs
`candidates -> screen -> construct` across parallel workers.

It is designed to be interrupted and re-run. Specifically:

- Targets with a final output already present are **skipped**. `--overwrite`
  redoes them.
- A failure in one gene does not kill the run. Failures are logged to
  `output/<gene>.log`, collected in `codebook.failed.txt`, and the command
  exits non-zero so a batch script notices.
- `codebook.acceptable.json` — the off-target allow-list produced by
  `mkprobes assemble ... short` — is picked up automatically, and forces the
  affected genes to be re-screened and re-constructed.
- Re-run a single target by naming it:
  `mkprobes run-panel data/mouse panel_a/codebook.json Sox2`. This forces
  overwrite for that target only.

Useful flags:

| Flag | Default | What it does |
| --- | --- | --- |
| `-j, --workers` | 16 | parallel worker processes |
| `-o, --output` | `output/` beside the codebook | where per-target files go |
| `--minimum` | 60 | minimum probes per gene at the screen stage |
| `--maxoverlap` | 0 | how much probe overlap to allow when reaching `--minimum` |
| `--target-probes` | 48 | maximum probes per gene at the construct stage |
| `--allow-file` | `<codebook>.acceptable.json` | per-gene acceptable off-targets |
| `--list-failed` | — | list targets with no final output, then exit |
| `--list-failed-all` | — | the same, plus each one's top off-target counts |

Match `--workers` to the CPUs you were actually allocated.

## Triaging what failed

```bash
mkprobes run-panel data/mouse panel_a/codebook.json --list-failed-all
```

This does no work; it reports which targets have no final output and what
their most common off-target binders were. That table is what tells you
whether a target failed because of one specific cross-reactive homolog (fixable
with `--allow`) or because the transcript is simply too short (not fixable —
pick a different isoform).

## The three commands underneath

Use these when tuning a single target. Running one gene at a time makes the
effect of each parameter obvious and wastes far less compute than re-running a
panel.

### 1. Candidates — the broad search

```bash
mkprobes candidates data/mouse --gene Sox2 --output panel_a/output
```

Enumerates every possible probe region on the target and works out what else
in the transcriptome each one would bind.

- `--allow gene1,gene2` — accept these as binders instead of rejecting the
  probe. Use after verifying the off-target is a genuine homolog you do not
  mind labelling.
- `--disallow gene1,gene2` — explicitly reject these.
- `--ignore-revcomp` — disables reverse-complement matching. Only with a clear
  biological reason.
- `--pseudogene-limit` — how many pseudogene hits to tolerate.
- `--overwrite` — redo an existing output.

On custom datasets: `--allow`/`--disallow` take **transcript IDs** (the FASTA
record IDs), not gene names. Sibling isoforms of the target's own gene are
allowed automatically from the GTF — without that, any multi-isoform gene
would yield zero probes. The rRNA/tRNA blocklist is enforced automatically
when the dataset carries one.

### 2. Screen — filter and select

```bash
mkprobes screen panel_a/output Sox2 --minimum 60 --maxoverlap 20 --overwrite
```

`OUTPUT_PATH` is the directory `candidates` wrote to — this step works on
the files it left there.

- `--minimum` — probes to aim for. Drives the adaptive overlap search.
- `-l, --overlap` — a fixed overlap. `--minimum` overrides it.
- `--maxoverlap` — how far the search may go to reach `--minimum`.
- `--restriction` — comma-separated. See the warning below.
- `--fpkm-path` — expression table used for weighting. Also accepted as
  `--fpkm_path`, the older spelling.

### 3. Construct — attach the readouts

```bash
mkprobes construct data/mouse panel_a/output --gene Sox2 --codebook panel_a/codebook.json
```

Both `--gene` and `--codebook` are **required**. This reads that target's
screened probes from the output directory and writes
`<target>_final_<enzymes>_<bits>.parquet` beside them.

- `-N, --target-probes` — maximum probes per gene (default 72 here; `run-panel`
  passes 48). Also accepted as `--target_probes`, the older spelling.
- `--restriction` — see below.

## About `--restriction`

It appears on all three commands, but it is **not a free choice**. SOLAR
chemistry fixes the pair to **BamHI + KpnI**: the header/footer sequences
carry those two sites, and final assembly excises the probe with a KpnI/BamHI
double digest. A different pair produces probes that nothing downstream can
cut out. Anything other than that pair is now refused up front, with an
explanation, rather than after the panel has been computed.

So: leave it alone. The default is already correct.

If you do write it out, the spelling differs between commands:

```bash
mkprobes screen    ... --restriction BamHI,KpnI                 # comma-separated
mkprobes run-panel ... --restriction BamHI,KpnI                 # comma-separated
mkprobes construct ... --restriction BamHI --restriction KpnI   # repeatable
```

## What lands in the output directory

Per target, in order:

```text
<target>_all.parquet          every candidate position
<target>_bowtie.parquet       raw alignment records
<target>_crawled.parquet      candidates with off-target context
<target>_crawled.stats.json
<target>_screened_ol*.parquet the selected probe pairs
<target>_screened_ol*.stats.json
<target>_final_BamHIKpnI_<bits>.parquet   the encoded constructs
```

Every parquet carries an embedded provenance record — version, timestamp,
command line, dataset and parameters. The `.stats.json` sidecars carry the
same record under a `provenance` key. Read it from any file with:

```bash
mkprobes provenance panel_a/output/Sox2_final_BamHIKpnI_1,2,3.parquet
```

What each column means, including the two genuinely surprising ones:
{doc}`../reference/columns`.

## Tuning loop

1. Run one target with the defaults.
2. Read its `_crawled.stats.json` for dominant off-target binders.
3. Decide the parameter policy — and fix it.
4. Run the whole panel with `run-panel`.

Changing parameters partway through a panel means targets were designed under
different rules. Prefer re-running the whole panel with `--overwrite`.

## On a cluster

- Keep `output/` on local scratch while running, then archive the final
  parquet files and logs.
- Match `--workers` to the job's CPU allocation.
- Shard by target list across jobs if one panel is too big for one job.

## When it goes wrong

- **Zero probes after filtering** — check the transcript is the one you meant
  (a short isoform gives few candidates), then look at the off-target table.
- **Outputs skipped that you wanted redone** — `--overwrite`, or name the
  single gene as the third argument to `run-panel`.
- **Sparse probe counts on many targets** — loosen `--minimum` /
  `--maxoverlap`, or accept verified homologs with `--allow`.
- **Zero probes for a multi-isoform gene on a custom dataset** — should not
  happen, siblings are auto-allowed. If it does, inspect
  `<target>_offtarget_counts.csv` for a homolog and `--allow` it.

---

Next: {doc}`qc_your_panel`.
