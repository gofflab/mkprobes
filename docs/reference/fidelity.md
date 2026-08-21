# Fidelity to the original implementation

`mkprobes` was extracted from
[chaichontat/fishtools](https://github.com/chaichontat/fishtools), where it
lived as `fishtools/mkprobes/` plus the driver scripts in `probegen/`. This page
records where the extracted code behaves differently from the original, and why.

It exists so that nobody has to re-derive this from a diff. If you are deciding
whether a panel designed today is comparable to one designed before the
extraction, this is the page to read.

The comparison was made against `fishtools@906c15d`, the last commit where the
code was in active use in-tree. The migration commit that followed introduced no
changes of its own.

## What is unchanged

The parts that determine whether a probe binds its target, and nothing else, are
identical:

| Component | Status |
| --- | --- |
| Candidate enumeration (`utils/_crawler.py`) | byte-identical file |
| Thermodynamics — Tm, hairpin, salt and formamide (`utils/seqcalc.py`) | identical |
| Bowtie 2 invocation and every flag (`utils/_alignment.py`) | identical |
| CIGAR / MD:Z parsing (`utils/_pairwise.py`) | byte-identical file |
| Sequence-quality criteria and the seven-tier cascade | identical |
| Readout attachment (`construct_encoding`, `construct_idt`) | identical |
| Codebook generation | identical output across every configuration tested |
| Vendored code matrices, readout and header/footer tables | 53 files, byte-identical |

Every design constant was compared individually: probe length 43–55 nt, GC 0.30–0.70
when tiling and 0.35–0.65 when filtering, Tm 54–68 °C, 390 mM Na⁺, 1 nM oligo,
maximum off-target Tm 20 °C, arm-split target Tm 60 °C, bowtie2 seed length 12
with `-k 200` and `--score-min L,24,0`. None of them changed.

Two empirical checks are worth quoting if you need to justify this in writing:
364 test sequences across 11 thermodynamic functions produced **zero** differences
between the original dependency versions and the current ones, and the refactored
probe-name parser agrees with the original on **351,170 real rows** of production
output.

## What is different

### Probe selection is no longer rounded to whole numbers

**Status: deliberate. The corrected behaviour is canonical.**

The routine that chooses which candidates get tiled onto a transcript kept its
dynamic-programming table in an integer array. Priorities are `sqrt()`-derived
floats, so every cumulative score written into that table was silently truncated
— the optimiser was solving a rounded version of the problem it meant to solve.

This could not be preserved even in principle: on NumPy 2 the original line
raises `OverflowError`, so the code does not run at all. The table is now
`float64`.

On real transcripts this preserves the **number** of probes selected but changes
**which** ones: roughly 10–30 % of positions are shared with the original. Both
sets are drawn from the same candidate pool after identical off-target and
thermodynamic filtering, so every probe in either set is equally specific and
equally active. What changes is which of the acceptable tiling positions get
used.

**Consequence: a panel designed before the extraction cannot be regenerated
exactly.** The lab has accepted this. If you need to reproduce a specific
historical order, use the oligo file from that order rather than re-running the
design.

### Ties in the assembly sort break by input order

**Status: deliberate. The new behaviour is canonical.**

The sort that picks the top *n* probes per gene now passes `maintain_order=True`.
Polars offers no ordering guarantee for an unstable sort, and the sort keys tie
often, so without this the selected probes could vary between runs. This changes
which probes are chosen relative to the original whenever the keys tie.

### Head splint and splint padding are drawn in row order

**Status: fixed. Output now matches the original exactly.**

Each padlock carries a 3 nt head splint drawn from one shared random generator,
and short splints are padded from one shared `ATAAT` cycle. Both were originally
consumed inside a polars UDF — and polars runs those across threads, so the
order of consumption was not reproducible. Three runs of the same panel produced
three different oligo pools.

Both are now drawn sequentially, in row order, outside the UDF. That reproduces
exactly what the original produced single-threaded, and produces it identically
at any thread count. `tests/test_assembly.py` pins the resulting sequences.

### The non-model-species path gained two protections

**Status: deliberate improvements, non-model path only. Mouse and human are
unaffected.**

Sibling isoforms of a target's own gene are now treated as acceptable binders.
Without this every multi-isoform gene fails screening, because each probe
"off-targets" its own siblings with near-perfect matches. The reference path
always allowed this through Ensembl.

The rRNA/tRNA k-mer blocklist is now applied. In the original it was active on
the reference path and commented out on the generic one, so nothing stopped a
probe for a non-model species landing in ribosomal RNA.

### Smaller differences

- `tm` and `hp` are stored as `Float64` rather than `Float32`. A schema change
  for downstream code; the numeric effect is around 1 part in 10⁶.
- Human GENCODE PAR_Y transcript IDs are no longer stripped to match FASTA keys.
  The affected genes fail loudly rather than producing bad probes, and only once
  the annotation cache is rebuilt.
- `filter-genes` crashed in the original on default-overlap filenames, and
  counted screened rather than final probes. Both are fixed; it now reports what
  a target contributes to an order.
- `CodebookPicker(existing=...)` now preserves row order instead of scrambling it
  through a set difference. Nothing in the shipped workflow passes `existing`, so
  this is latent.
- Extending a codebook, validating its inputs, and resolving a dataset path all
  reject malformed input the original silently accepted.

## Known issues inherited from the original

These are present in both versions and have **not** been changed, because fixing
them would alter probe selection:

- `find_overlap` (used for the top priority tier) and `find_overlap_weighted`
  (used for the rest) apply opposite sign conventions to the spacing parameter,
  so the same `--overlap -2` means a 2 nt gap to one and a 2 nt overlap to the
  other.
- In `handle_overlap`, `selected_global` is updated after the tier loop rather
  than inside it, so only the final tier's selection survives and the
  already-selected filter inside the loop never matches anything.

The second largely masks the first. Both are tracked for a decision.
