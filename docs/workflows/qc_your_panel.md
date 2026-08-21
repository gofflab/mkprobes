# QC your panel

**Step 5 of the {doc}`../getting_started` workflow.**

Every target succeeding individually does not mean the panel is good. A target
that scraped through with 12 probes will be dim and unreliable on the
microscope, and you want to know that before you pay for oligos, not after.

## Count what each target actually got

```bash
mkprobes filter-genes panel_a/output --genes panel_a/genes.converted.tss.txt \
    --min-probes 48 --out panel_a/genes.pass.txt
```

`--genes` is **required** — the command needs to know which targets were
supposed to exist, not just which files happen to be in the directory.

The count comes from each target's `_final_` file: the constructed, encoded
probes, which is exactly what that target would contribute to an oligo order.
(It is not the `_screened_` count, which is larger — screening selects probe
pairs, and construction then caps and encodes them.)

Three kinds of result:

- **At or above `--min-probes`** — fine. The comparison is `>=`, so a target
  with exactly 48 passes at `--min-probes 48`.
- **Below the threshold** — warned about individually, by name and count.
  These are the ones to rework.
- **No constructed probes at all** — reported separately as an error, with the
  count and the first several names. These are not thin targets; they are
  `run-panel` failures that never completed. Go back to
  {doc}`design_probes` and check `codebook.failed.txt` before treating them as
  a QC problem.

`--out` writes the passing targets, one per line, so the next step has an
explicit list rather than an implicit one.

## What threshold to use

48 is the usual floor for a bright, reliably detected target. Below roughly 30
you should expect the target to be unreliable rather than merely dim. The
right number depends on your expression levels and imaging setup; pick one for
the panel and apply it consistently.

## Fixing the thin ones

For each target below threshold, in rough order of what to try:

1. **Check the transcript.** A short isoform simply has less room for probes.
   Try `-m longest` in {doc}`choose_your_targets` and regenerate.
2. **Look at what it was losing probes to.** Run
   `mkprobes run-panel ... --list-failed-all`, or read the target's
   `_crawled.stats.json`. One dominant cross-reactive binder is a different
   problem from diffuse loss.
3. **Accept verified off-targets.** If the dominant binder is a homolog you do
   not mind labelling, `--allow` it — or use the interactive triage in
   {doc}`order_your_oligos`, which records your decisions in
   `codebook.acceptable.json` and lets `run-panel` apply them automatically.
4. **Loosen screening.** Raise `--maxoverlap` so probes may overlap slightly
   to reach the count.
5. **Drop the target.** Sometimes the honest answer. Remove it from the target
   list, regenerate the codebook, and re-run — do not simply delete it from
   the codebook by hand.

Then re-run `filter-genes` and check again.

## Where this leaves you

The per-target `_final_` parquet files are the design deliverable:

```text
panel_a/output/<target>_final_BamHIKpnI_<bits>.parquet
```

They are not orderable yet — turning them into oligos is
{doc}`order_your_oligos`.

To see exactly how one was made:

```bash
mkprobes provenance panel_a/output/Sox2_final_BamHIKpnI_1,2,3.parquet
```

Column meanings: {doc}`../reference/columns`.

## On a cluster

- QC is lightweight; run it as a post-job step, not its own allocation.
- Archive the codebook, target list, `genes.pass.txt` and the final parquet
  files together as one panel release bundle.
- Keep rework jobs scoped to the specific failing targets.

---

Next: {doc}`order_your_oligos`.
