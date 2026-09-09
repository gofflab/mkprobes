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
| `--manifest` | `manifest.json` beside the codebook | where the panel's design settings are read from |
| `--minimum` | 60 | probes per gene the screen aims for |
| `--maxoverlap` | 0 | how far probes may overlap (nt, multiples of 5) to reach `--minimum` |
| `--tm-range` | `54,68` | crawler Tm window, °C at the design formamide |
| `--length-range` | `43,55` mouse/human, `43,54` otherwise | probe length window, nt (60 at most) |
| `--split-tm` | 60 | Tm each arm of the split probe must reach, °C |
| `--target-probes` | 48 | recorded in provenance; the pool cap is `n_probes` in the manifest |
| `--allow-file` | `<codebook>.acceptable.json` | per-gene acceptable off-targets |
| `--list-failed` | — | list targets with no final output, then exit |
| `--list-failed-all` | — | the same, plus each one's top off-target counts |

Match `--workers` to the CPUs you were actually allocated.

The five design settings (`--minimum` through `--split-tm`) live in the
manifest's `design` block, and that is where a panel's values belong: see
[Design settings](#design-settings) below. A flag given on the command line
overrides the manifest for that run, and the run log says which came from
where.

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
- `--tm-range`, `--length-range`, `--split-tm` — the thermodynamic settings,
  as on `run-panel`. Candidates left over from a run under other settings are
  redone rather than reused.

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

- `--minimum` — probes to aim for. Drives the adaptive overlap search: the
  screen tiles at no overlap first, then at 5, 10, ... nt of overlap until
  the count is reached or `--maxoverlap` is hit, writing one
  `_screened_ol<N>_` file per overlap tried.
- `-l, --overlap` — a fixed overlap. `--minimum` overrides it.
- `--maxoverlap` — how far the search may go to reach `--minimum`. Note the
  default differs here (20) from `run-panel` (0).
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

- `--overlap` — which screened file to build from, by its overlap. Left out,
  it takes the one the screen settled on: the largest overlap present, which
  is the first that reached `--minimum`. Before this, construct always read
  the no-overlap file, so `--maxoverlap` produced files nothing used.
- `-N, --target-probes` — recorded in the output's provenance. The number of
  probes that reach the pool is capped by `n_probes` in the manifest at
  assembly, not here. Also accepted as `--target_probes`, the older spelling.
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

## Design settings

Five settings decide how many probes a target can yield before any off-target
check runs. They were calibrated on mammalian transcripts at about 50% GC,
and they are the reason an AT-rich transcriptome (cephalopods sit near 36%
GC) returns thin panels: a window that cannot reach the Tm floor within the
length cap is dropped, and only local GC-rich islands yield candidates.

They live in the manifest, one block per probe set, so a panel is designed
and assembled under one recorded set of values. `mkprobes init` writes the
block with the defaults filled in:

```json
"design": {
  "tm_range": [54, 68],
  "length_range": [43, 54],
  "split_tm": 60,
  "min_probes": 60,
  "max_overlap": 0
}
```

| Setting | What it does | Lever on AT-rich transcripts |
| --- | --- | --- |
| `tm_range` | Tm window (°C at the design formamide) a probe must fall in | lower the floor to admit probes that bind less tightly |
| `length_range` | probe length window, nt | raise the cap (60 at most) so a window can grow long enough to reach the Tm floor |
| `split_tm` | Tm each arm of the split probe must reach | the largest lever, and the one to lower with the most care: both arms must bind for ligation |
| `min_probes` | probes per gene the screen aims for | |
| `max_overlap` | how far neighbouring probes may overlap to reach `min_probes` | 10 or 20 buys probes without touching thermodynamics |

Delete a field to keep its default. `run-panel` reads the block for the
probe set that names its codebook (a `manifest.json` beside the codebook, or
`--manifest`), a flag overrides one setting for one run, and `assemble`
warns when a target's output was designed under anything other than the
manifest's block. Different panels in one manifest can carry different
blocks.

Measured on a 17.9 kb octopus transcript at 34% GC, with off-target
screening on:

| Settings | Screened pairs |
| --- | --- |
| defaults | 58 |
| `split_tm` 55 | 67 |
| `split_tm` 55, `length_range` 43–60, `tm_range` 50–68 | 126 |
| defaults, `max_overlap` 20 | 80 |

Widening the GC window does nothing (GC is not the check that fails), and
lowering the formamide alone makes things worse. Every thermodynamic change
alters how the probes hybridise at the bench, so treat the manifest block as
a decision about chemistry, not a knob to turn until the count looks right.

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
- **Sparse probe counts on many targets** — on an AT-rich transcriptome, the
  [design settings](#design-settings) are the cause; otherwise raise
  `max_overlap`, or accept verified homologs with `--allow`.
- **Changing the design block did nothing** — finished targets are skipped.
  The run warns which were designed under other settings; pass `--overwrite`
  to redesign them.
- **Zero probes for a multi-isoform gene on a custom dataset** — should not
  happen, siblings are auto-allowed. If it does, inspect
  `<target>_offtarget_counts.csv` for a homolog and `--allow` it.

---

Next: {doc}`qc_your_panel`.
