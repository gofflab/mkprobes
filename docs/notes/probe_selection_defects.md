# Probe-tiling selection: two inherited defects

Both defects predate the port and are present verbatim in the upstream
`fishtools` code.

| | status |
| --- | --- |
| Defect 1 — `overlap` sign convention | **fix approved, not yet applied** — awaiting confirmation against a real crawl. Panel-changing. |
| Defect 2 — `selected_global` outside its loop | **done** — dead code removed, `RecursionError` path guarded. Behaviour-preserving. |

Reproduce any number below, and produce the real per-gene figures defect 1 is
waiting on, with:

```bash
python scripts/quantify_overlap_defects.py data/och_test_output --overlap -2
```

The script reports per-gene probe counts under four variants (baseline, sign
fix, accumulate, both) and checks on every run that its instrumented
`handle_overlap` still reproduces the shipped `the_filter` exactly.

---

## Defect 1 — `overlap` has opposite sign conventions in the two selectors

`src/mkprobes/utils/_algorithms.py`:

| selector | compatibility test | used for |
| --- | --- | --- |
| `OverlapWeighted.q` (line 34) | `bisect_left(end[:j-1], start[j-1] + overlap)` → `prev_end < start + overlap` | priority tiers 2‑7 |
| `find_overlap` (line 70) | `start[i] - overlap > end[out[-1]]` → `prev_end < start - overlap` | priority tier 1 only |

The two agree only at `overlap == 0`. For every other value they mean opposite
things.

### The parameter's intended meaning is not ambiguous

`overlap` is meant to be *permitted overlap in nt* — higher is more permissive:

- `--maxoverlap` is documented as "Maximum sequence overlap between probes"
  (`docs/reference/cli.md`, `screen.py`).
- `run_screen(minimum=...)` escalates `overlap` through `(-2, 5, 10, 15, 20)`
  specifically to *increase* the probe count until `--minimum` is met.
- `OverlapWeighted.q`'s own comment: "Guarantees that overlap = 0 means no
  overlap."

Raising `overlap` must therefore yield more probes. Measured on a dense tiling
(50 nt probes every 5 nt):

```
overlap:      -10   -5   -2    0    5   10   20
weighted:       7    8    8    8    9   10   14     <- relaxes, as intended
greedy:        10    9    8    8    8    7    6     <- tightens, inverted
```

`find_overlap` responds to the parameter backwards. At the production default
`overlap=-2`, the weighted selector enforces a 2 nt **gap** while the greedy
selector permits a 2 nt **overlap**:

```
probes A=(0-49) B=(48-97) C=(52-101)
  find_overlap(overlap=-2)          -> [0, 1]   A+B, sharing 2 nt
  find_overlap_weighted(overlap=-2) -> [0, 2]   A+C, 2 nt gap
```

Two adjacent padlock probes cannot both bind overlapping target sequence, so the
weighted convention is the assay-correct one.

### Impact

`find_overlap` is only reached for priority tier 1, and its result only survives
if the tier loop breaks at `i == 1`, which requires tier 1 alone to yield more
than `n` (default 100) probes. So the defect is **live only for long,
high-quality genes**.

On a proxy corpus of long high-quality transcripts (9–16 kb), 7 of 10 genes
break at tier 1:

| gene | baseline pairs | sign-fixed pairs | probes changed | adjacent pairs sharing target seq (baseline) |
| --- | --- | --- | --- | --- |
| Long.700.1 | 154 | 149 | +32 / −37 | 25 (2 nt each) |
| Long.701.1 | 122 | 117 | +14 / −19 | 17 |
| Long.702.1 | 184 | 171 | +28 / −41 | 28 |
| Long.704.1 | 132 | 127 | +29 / −34 | — |
| Long.706.1 | 173 | 163 | +31 / −41 | — |
| Long.707.1 | 165 | 160 | +41 / −46 | — |
| Long.708.1 | 146 | 141 | +14 / −19 | — |

Fixing the sign removes every physical overlap (25→0, 17→0, 28→0) at a cost of
~3–7 % fewer probes. On a corpus of ordinary-length genes (0.7–5 kb) **no gene
breaks at tier 1 and the fix changes nothing at all** — see defect 2.

### The escalation ladder is inverted for exactly the genes it should rescue

`run_screen --minimum N --maxoverlap 20` walks `overlap` upward expecting more
probes. For tier-1-breaking genes it currently gets fewer:

```
gene          variant      ol=-2  ol=5  ol=10  ol=15  ol=20
Long.700.1    baseline       154   147    137    129    122   <- inverted
Long.700.1    sign-fixed     149   157    163    181    196
Long.702.1    baseline       184   169    161    153    146   <- inverted
Long.702.1    sign-fixed     171   190    209    224    243
```

Production defaults (`--minimum 60 --maxoverlap 0`) disable the ladder —
`chain((-2,), range(5, 1, 5))` is just `(-2,)` — so this bites only workflows
that pass a non-zero `--maxoverlap`, such as
`docs/workflows/phase_3_candidate_screen_construct.md`.

### Decision — **fix, after confirming on the real crawl**

Align `find_overlap` with `OverlapWeighted.q`:

```python
if start[i] + overlap > end[out[-1]]:      # was: start[i] - overlap
```

It is the assay-correct convention, it is the one the parameter is documented
and used as, and it makes the two selectors agree. On ordinary-length genes it
is free — zero probes change. It only affects long high-quality genes, where it
removes real overlapping-probe pairs at a cost of a few percent of probes.

**Not yet applied.** The proxy numbers above establish the direction and the
mechanism, but the size of the change on the real panel depends on how many
genes break at tier 1, which is a property of the actual crawl. Run the script
against `data/och_test_output` and read the summary line:

```
break at tier 1 (greedy find_overlap reaches the panel; defect 1 is live): N/M
```

`N == 0` means the fix is a no-op on that panel and can land with no golden-test
regeneration. `N > 0` means those `N` genes' oligos will change and
`tests/test_assembly.py` needs regenerating — which is the intended signal to
confirm the change is wanted before it lands.

`tests/utils/test_algorithms.py` already pins the corrected convention behind
`xfail(strict=True)`. Applying the one-line fix turns those nine tests into
XPASS, which is the prompt to drop the marks and delete
`test_greedy_currently_permits_overlap_at_negative_overlap`.

---

## Defect 2 — `selected_global` is accumulated outside its loop

`src/mkprobes/utils/_filtration.py`, `handle_overlap`: the tier loop computes
`sel_local` per tier, but `selected_global |= sel_local` sits *after* the loop
(line ~108). Consequences as written:

1. Only the final tier's `sel_local` survives.
2. The in-loop `~pl.col("index").is_in(selected_global)` filter is a permanent
   no-op, because `selected_global` stays empty for the whole loop.
3. `sel_local` is read after the loop behind a `# type: ignore`.

### The loop was not meant to accumulate

The obvious reading — move `selected_global |= sel_local` inside the loop so
each tier tops up the panel — was tested and is **wrong**. Each tier re-runs the
selector on a candidate pool with the already-selected probes *removed*, so the
selector has no knowledge of the positions they occupy and happily picks probes
on top of them:

| gene | transcript | max non-overlapping | baseline | "accumulate" | adjacent pairs overlapping | worst overlap |
| --- | --- | --- | --- | --- | --- | --- |
| Och.600.1 | 4 355 nt | ~96 | 39 | 119 | 89 | 53 nt |
| Och.601.1 | 3 397 nt | ~75 | 35 | 119 | 96 | 53 nt |
| Och.603.1 | 1 274 nt | ~28 | 10 | 47 | 39 | 53 nt |
| Och.605.1 | 2 526 nt | ~56 | 19 | 101 | 86 | 53 nt |

Probes are 43–54 nt, so a 53 nt overlap is a near-duplicate probe. Och.603.1
would get 47 pairs onto a transcript that physically holds 28. The baseline has
**zero** overlap violations.

What the loop actually implements is *progressive relaxation with restart*: tier
`i` selects over all probes of priority ≤ `i`, from scratch, so each iteration is
a complete self-consistent tiling drawn from a strictly larger pool, and the loop
stops at the first tier yielding more than `n` probes. Under that reading the
post-loop `selected_global |= sel_local` is exactly right and `selected_global`
is simply vestigial. The weighted selector already handles quality preference via
the `sqrt(len(criteria) + 1 - priority)` weights — it does not need a
carry-forward set.

### Impact

Behaviourally, none — the current output is the intended output. The
consequences are code-health only:

- Dead code: `selected_global` and the `~is_in(...)` filter.
- Sub-claim 3 needs a correction: `sel_local` cannot actually be unbound via
  "every tier hits `continue`". `filter_have_both` guarantees a splint row at
  some tier, so some tier always produces a `run`. The only route to
  `UnboundLocalError` is tier 1 `continue` followed by a `RecursionError` at a
  later tier — and that could not be triggered with realistic input:
  `OverlapWeighted.__init__` raises the limit to 5000 itself, and `dp`'s recursion
  depth is bounded by candidate density within one probe span (~50), not by the
  candidate count. Dense inputs of 5 200 candidates completed fine. It is a
  latent robustness wart, not a reachable crash.

### Decision — **cleaned up, not "fixed"** (done)

Applied in `handle_overlap`:

- `selected_global` and the no-op `~is_in(...)` filter deleted; the tier result
  is now a single `selected` set initialised before the loop, and the `# type:
  ignore` is gone.
- The `RecursionError` branch now logs which tier failed, records it in `stats`,
  and re-raises if no tier ever completed rather than falling through to an
  unbound name. Where a tier had completed, it still falls back to that tier —
  unchanged.
- A comment records why the tiers must not accumulate, so the next reader does
  not re-derive it.

Verified behaviour-preserving: panels are byte-identical across 34 proxy genes
before and after, and the script's fidelity check still reports `OK`.
`stats` keys are unchanged in normal operation, so `.stats.json` output is
unaffected. The accumulation was explicitly **not** moved inside the loop.

---

## Separate observation (neither defect)

The `if len(sel_local) > n: break` early stop can itself reduce the panel: a
tier-1 result of 105 probes stops the search, where falling through to tier 7
would have yielded 118. That is a deliberate quality-over-quantity trade, but it
makes probe count non-monotonic in `overlap` even after defect 1 is fixed
(visible as Long.703.1 in the ladder table). Worth knowing when tuning
`--minimum`.

---

## Caveat on these numbers

`data/och_test_output/` is not reachable from the environment this analysis ran
in, so the tables above come from a proxy corpus: synthetic AT-rich transcripts
(*O. chierchiae* is ~40 % GC) put through the **real** `crawler`, with `oks`,
`hp`, `tm` and `gc_content` computed by the **real** `candidates.py` expressions.
Only `max_tm_offtarget` and `maps_to_pseudo` are simulated, since those need a
bowtie index of the real transcriptome.

The instrumented `handle_overlap` in `scripts/quantify_overlap_defects.py` is
verified to reproduce the shipped `the_filter` exactly on every input it is run
against (`baseline fidelity vs shipped the_filter: OK`), so running the script on
the real crawl will give the real per-gene numbers directly.

The proxy corpus is enough to establish the *mechanism* and the *direction* of
each defect, which is what the recommendations rest on. It is not a substitute
for the real per-gene counts, which is why defect 1 is held pending them.
