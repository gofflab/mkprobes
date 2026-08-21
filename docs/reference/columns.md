# Output columns, stage by stage

The probe-design workflow writes one parquet file per target per stage. Those
files carry up to 44 columns, seven of which are DNA sequences, and the column
names are terse. This page says what every one of them means.

Parquet files are not human-readable in a text editor. To look at one, use
Python with `polars` inside your activated environment:

```python
import polars as pl
df = pl.read_parquet("output/MYGENE_crawled.parquet")
print(df.schema)      # column names and types
print(df.head())
```

Terms in **bold** are defined in {doc}`../glossary`; the assay itself is
explained in {doc}`../what_is_solar`.

## Read this first: `seq` does not mean the same thing twice

`seq` is the column people reach for first, and it is the single most
confusing thing about these files. **It changes meaning at every stage.**

| File | What `seq` holds | Typical length |
| --- | --- | --- |
| `*_crawled.parquet` | One **half** of a candidate window — either the padlock's binding arm or the splint's, one per row | 18–29 nt |
| `*_screened_ol*.parquet` | The **whole** binding window for the pair, in target orientation | 86–110 nt — twice the window (see below) |
| `*_final_*.parquet` | The **encoded payload**: the three readout sequences plus the padlock's binding arm | 82–91 nt |

Three consequences:

- In `_crawled`, one probe pair occupies **two rows** — one named
  `..._splint`, one named `..._padlock`. From `_screened` onward, a pair
  occupies **one row**.
- In `_screened`, `seq` is currently the full window written **twice**, once
  contributed by each member of the pair, so its length is exactly double the
  real window length. The window itself is the first half of the string. This
  is a quirk of how the pair is reassembled, not extra sequence. If you want
  the window, take `seq[:len(seq)//2]`, or use `pos_start`/`pos_end`.
- In `_final`, the old `_screened` value is preserved under the name
  `seqori`, and `seq` is the new encoded sequence. If you are comparing
  stages, compare `_final.seqori` against `_screened.seq`.

The sequences you would actually hand to a synthesis vendor are in neither
place: they are produced later, by the assembly step, and written to
`generated/<panel>_final.txt`.

## And the second surprise: paired columns are lists

From `_screened` onward, one row is one probe *pair*. Every per-half value is
therefore collected into a **two-element list**:

```text
priority : [1, 6]
oks      : [6, 3]
tm       : [43.2, 39.8]
```

The two entries are the pair's two halves — but **the order is not
guaranteed** to be splint-then-padlock, and in practice it varies row to row.
Treat these as an unordered pair. In particular, to reduce them to one number
per pair use an explicit aggregate (`.list.min()`, `.list.max()`,
`.list.mean()`), never `.list.get(0)`. The package itself does exactly this,
sorting probes by `priority.list.max()` and `hp.list.min()`.

Only `name`, `gene`, `pad_start`, `splint`, `padlock`, `seq` (and in `_final`,
`code1`/`code2`/`code3` and `seqori`) are plain scalars.

## `<target>_crawled.parquet` — the candidate pool

One row per candidate half that survived the initial alignment screen. Expect
thousands of rows for a normal-length transcript; they overlap heavily,
because the crawler emits a candidate from nearly every starting position.
Nothing here has been selected yet.

| Column | Type | Meaning |
| --- | --- | --- |
| `name` | String | Full probe-half identifier: `{gene}_{transcript}:{start}-{end}_{splint\|padlock}`. The two halves of a pair share everything but the suffix. |
| `seq` | String | This half's binding sequence, in target orientation (not yet reverse-complemented). 18–29 nt. |
| `seq_full` | String | The complete window both halves came from, before splitting. 43–55 nt. |
| `gene` | String | Gene name parsed from `name`. For custom-species datasets this equals the transcript ID. |
| `transcript_ori` | String | The transcript this probe was **designed against**, parsed from `name`. Contrast with `transcript` below. |
| `pos_start` | UInt32 | 0-based start of the full window within the transcript sequence. |
| `pos_end` | UInt32 | 0-based end of the full window, inclusive. |
| `length` | UInt32 | `pos_end - pos_start + 1` — the length of the **full window**, not of `seq`. |
| `pad_start` | Int16 | Offset within `seq_full` at which the second (splint) half begins. Both halves of a pair carry the same value. The construct stage requires it to be greater than 17. |
| `transcript` | String | The transcript this row's **alignment** landed on. In `_crawled` this is always an acceptable target (the intended transcript or one of its sibling isoforms). |
| `id` | UInt32 | Line number of the underlying SAM record. Internal bookkeeping. |
| `flag`, `pos`, `cigar`, `aln_score`, `aln_score_best`, `n_ambiguous`, `n_mismatches`, `n_opens`, `n_extensions`, `edit_distance`, `mismatched_reference` | see below | Standard alignment fields carried through from `bowtie2`. Explained in [The alignment columns](#the-alignment-columns). |
| `match` | Int64 | Number of bases that matched the reference in this alignment (the sum of the run lengths in `mismatched_reference`). |
| `match_consec` | UInt8 | Longest uninterrupted run of matching bases in this alignment. |
| `match_consec_all` | UInt8 | The worst (largest) `match_consec` this half achieved against any **off-target** transcript. `0` means nothing off-target matched substantially. |
| `max_tm_offtarget` | Float32 | The worst-case off-target melting temperature, in °C. See [Off-target columns](#off-target-columns). `0.0` means no meaningful off-target hit. |
| `maps_to_pseudo` | String | For human/mouse datasets, the ID of an allow-listed extra transcript (typically a **pseudogene** you accepted) that this probe also binds. Empty string for custom-species datasets, which do not run pseudogene triage. |
| `gc_content` | Float64 | Fraction of `seq` that is G or C, as a number between 0 and 1. |
| `tm` | Float64 | Predicted melting temperature of `seq` against its RNA target, in °C, RNA:DNA nearest-neighbour model, 390 mM Na, formamide-corrected. Lower than an uncorrected calculator would give. |
| `hp` | Float64 | **Hairpin** score: predicted melting temperature in °C of the most stable self-structure `seq` can form. Higher is worse. Strongly negative values (around −29 in a typical run) mean no hairpin was found at all. This is **not** a homopolymer flag — see `ok_homopolymer`. |
| `ok_quad_c` | Boolean | Sequence-quality flag. **The name is misleading — see below.** |
| `ok_quad_a` | Boolean | Sequence-quality flag. **The name is misleading — see below.** |
| `ok_stack_c` | Boolean | Sequence-quality flag. |
| `ok_comp_a` | Boolean | Sequence-quality flag. **The name is misleading — see below.** |
| `ok_homopolymer` | Boolean | Sequence-quality flag. |
| `ok_gc` | Boolean | Sequence-quality flag. |
| `oks` | UInt32 | How many of the six `ok_*` flags are `true`. Ranges 0–6; higher is better. |

## `<target>_screened_ol*.parquet` — the selected probe pairs

One row per selected probe **pair**. Expect tens of rows, tiled along the
transcript with a small gap between neighbours. This is the file panel QC
counts, and the one the construct stage reads.

It carries every `_crawled` column except `seq_full` (which became `seq`),
with the per-half values collected into lists, plus these:

| Column | Type | Meaning |
| --- | --- | --- |
| `name` | String | The pair's identifier: `{gene}_{transcript}:{start}-{end}`, with no `_splint`/`_padlock` suffix. |
| `splint` | String | The splint's target-binding half, in target orientation. 18–29 nt. |
| `padlock` | String | The padlock's target-binding half, in target orientation. 18–27 nt. |
| `seq` | String | The full binding window, written twice (see the warning above). |
| `priority` | List(UInt8) | Which tier of the filter cascade this half satisfied: `1` is the strictest, `7` the most permissive. See [The filter cascade](#the-filter-cascade). |
| `index` | List(UInt32) | Row number used internally by the tiling algorithm. Not meaningful outside it. |
| `gene`, `pad_start` | String, Int16 | Scalars, as in `_crawled`. |

Everything else — `id`, `flag`, `transcript`, `pos`, `cigar`, `aln_score`,
`aln_score_best`, `n_ambiguous`, `n_mismatches`, `n_opens`, `n_extensions`,
`edit_distance`, `mismatched_reference`, `transcript_ori`, `pos_start`,
`pos_end`, `length`, `match`, `match_consec`, `maps_to_pseudo`,
`max_tm_offtarget`, `match_consec_all`, all six `ok_*` flags, `gc_content`,
`tm`, `hp`, `oks` — is the `_crawled` column wrapped in `List(...)`.

## `<target>_final_*.parquet` — the encoded constructs

One row per constructed probe pair, after the gene's three readout sequences
have been stitched onto the padlock. Slightly fewer rows than `_screened`:
the construct step tries four different spacer choices per probe, and drops
the probe entirely if all four would produce a run of five identical bases
(`AAAAA`, `TTTTT`, `CCCCC`, `GGGGG`) or a BamHI site in the joined sequence.

| Column | Type | Meaning |
| --- | --- | --- |
| `name` | String | Same pair identifier as in `_screened`. |
| `seq` | String | The encoded payload: the three readout sequences (lowercase, reverse-complemented) joined by two-nucleotide spacers, followed by the padlock's binding half (uppercase). 82–91 nt. |
| `code1` | Int64 | First **readout** ID written onto this probe. |
| `code2` | Int64 | Second readout ID. |
| `code3` | Int64 | Third readout ID. |
| `seqori` | String | The `_screened` `seq` value, carried through unchanged. |
| `splint`, `padlock`, `pad_start`, `gene` | String / Int16 | Unchanged from `_screened`. |

All other columns are the `_screened` list columns, unchanged.

About `code1`/`code2`/`code3`: together they are the gene's codebook entry.
Every probe for a gene gets **all three** of its bits — it is not the case that
some probes carry one bit and others carry another. Their *order* varies from
row to row, because the construct step cycles through the permutations of the
triple looking for one that does not create a forbidden motif. So
`{code1, code2, code3}` as a set is the meaningful thing, and it will match the
bits token in the filename. Readout IDs run from 1 to 49.

## The alignment columns

These seven columns come straight from `bowtie2`'s SAM output and are not
specific to this pipeline. If you have never read a SAM file, this is all you
need:

`flag`
: A bitfield describing the alignment. In these files you will see `0`
  (aligned to the forward strand, best alignment for this probe), `16`
  (aligned to the reverse strand), `256` (a secondary alignment — the same
  probe also aligned somewhere else, and this row is one of the other
  places), and `272` (secondary *and* reverse strand).

`pos`
: The 1-based leftmost coordinate of the alignment **on the transcript named
  in `transcript`**. Do not confuse it with `pos_start`, which is a coordinate
  on the transcript the probe was designed against.

`cigar`
: A compact description of the alignment shape. `23M` = 23 aligned bases.
  `1S21M4S` = 1 base soft-clipped at the front, 21 aligned, 4 clipped at the
  end — a partial match, which is what most off-target hits look like.

`aln_score`
: bowtie2's alignment score for *this* alignment (`AS:i` in SAM). Higher is a
  better match.

`aln_score_best`
: **Misleading name.** This is bowtie2's `XS:i` field: the score of the best
  alignment found for this probe *other than* the one on this row. It is
  present only when the probe aligned more than once. When it equals
  `aln_score`, the probe has at least two equally good placements — usually
  sibling isoforms of the same gene.

`n_mismatches`, `n_ambiguous`, `n_opens`, `n_extensions`, `edit_distance`
: Respectively `XM:i` (mismatched bases), `XN:i` (ambiguous reference bases
  overlapped), `XO:i` (gaps opened), `XG:i` (gap extensions), and `NM:i`
  (total edit distance to the reference).

`mismatched_reference`
: The SAM `MD:Z` field, describing the alignment from the reference's point of
  view: numbers are runs of matching bases, letters are the reference base at
  a mismatch. `23` means 23 matches in a row; `10A12` means 10 matches, a
  mismatch where the reference had an A, then 12 more matches. The `match` and
  `match_consec` columns are computed directly from this string — `match` is
  the sum of the numbers, `match_consec` is the largest of them.

## Off-target columns

"Off-target" is not recorded as a yes/no flag. Every candidate is aligned
against the entire transcriptome, alignments to acceptable transcripts (the
target itself plus its sibling isoforms plus anything you explicitly allowed)
are set aside, and the remaining hits are summarised into two numbers.

`match_consec_all`
: The longest uninterrupted match this probe half achieved against any
  unacceptable transcript. `0` means nothing off-target matched at all.

`max_tm_offtarget`
: The important one. For every off-target alignment with a run of more than 16
  consecutive matches, the pipeline finds the longest perfectly matched
  stretch of at least 15 bases and computes its melting temperature; this
  column holds the highest such value across all off-target hits, in °C.
  `0.0` means no off-target hit was substantial enough to be worth scoring.
  This is a thermodynamic measure, not a count: a 25-base off-target match at
  high GC is far more dangerous than a 25-base match at low GC, and this
  column reflects that. The default screening threshold is 20 °C.

`maps_to_pseudo`
: A separate, human/mouse-only channel for hits you decided to tolerate. It
  holds the ID of an allow-listed transcript (typically a pseudogene or close
  paralog) that this probe also binds. The two strictest filter tiers require
  it to be empty, so a probe with a pseudogene hit can still be used, but only
  once the stricter tiers have been exhausted.

## The `ok_*` quality flags, precisely

Six boolean checks are applied to each half's `seq`. `oks` counts how many
passed. **Three of the names do not describe what they test.** Those names are
scheduled to be changed but have not been changed yet, so read the table, not
the name.

| Flag | What it actually tests | Name accurate? |
| --- | --- | --- |
| `ok_quad_c` | `true` when `seq` contains **no `GGGG`**. | **No.** The name says `c`; the test is on **G**. |
| `ok_quad_a` | `true` when `seq` contains **no `TTTT`**. | **No.** The name says `a`; the test is on **T**. |
| `ok_stack_c` | `true` when no 6-nucleotide window within the **last 11 nucleotides** of `seq` contains 4 or more G's. Six overlapping windows are checked, ending at the 3' end. | Partly — again G, not C. |
| `ok_comp_a` | `true` when **fewer than 28% of the bases are T**. | **No.** The name says `a`; the test counts **T**. |
| `ok_homopolymer` | `true` when `seq` contains no run of 4 or more identical bases (no `AAAA`, `TTTT`, `CCCC` or `GGGG`). | Yes. |
| `ok_gc` | `true` when `gc_content` is between 0.35 and 0.65 inclusive. | Yes. |

Note that `ok_homopolymer` implies both `ok_quad_c` and `ok_quad_a`, so the
flags are not independent. `oks` therefore ranges 0–6, and 6 means a clean
sequence on every count.

And once more, because it catches everyone: **`hp` is the hairpin score, not
a homopolymer score.** Homopolymers are `ok_homopolymer`.

## The filter cascade

`priority` records which tier of a seven-step cascade a probe half satisfied.
Tier 1 is strictest; each subsequent tier relaxes the requirements. Screening
first tries to tile the transcript using only tier-1 probes; if that does not
yield enough, it retries allowing tiers 1–2, then 1–3, and so on.

| `priority` | Requirements |
| ---: | --- |
| 1 | `oks` > 5 and `hp` < 32 and `max_tm_offtarget` < 20 and no pseudogene hit |
| 2 | `oks` > 4 and `hp` < 32 and `max_tm_offtarget` < 20 and no pseudogene hit |
| 3 | `oks` > 4 and `hp` < 32 and `max_tm_offtarget` < 20 |
| 4 | `oks` > 4 and `hp` < 37 and `max_tm_offtarget` < 24 |
| 5 | `oks` > 3 and `hp` < 37 and `max_tm_offtarget` < 24 |
| 6 | `oks` > 2 and `hp` < 37 and `max_tm_offtarget` < 24 |
| 7 | `oks` > 1 and `hp` < 37 and `max_tm_offtarget` < 24 |

A candidate that satisfies none of the seven is discarded entirely, so
`priority` is never 0 in a written file. A panel dominated by priority 6–7
probes is a warning sign: it means the strict tiers could not fill the
transcript, and you should look at the gene's `*.stats.json` before trusting
it.

## Decoding the filenames

Filenames encode the parameters that produced them, which is how you can tell
two runs apart in the same directory.

```text
Och.958.1_crawled.parquet
Och.958.1_screened_ol-2_BamHIKpnI.parquet
Och.958.1_final_BamHIKpnI_4,5,6.parquet
```

`<target>`
: The transcript the probes were designed against. For human/mouse reference
  datasets this is the Ensembl transcript *name*; for custom-species datasets
  it is the transcript ID exactly as it appears in your annotation — which is
  why these names can contain dots and underscores.

`_crawled`
: The raw candidate pool.

`_screened_ol<overlap>`
: The selected, tiled probe pairs. `<overlap>` is the spacing parameter, and
  **it is negative by default (`-2`)**.

`_<enzymes>`
: The restriction enzymes candidates were screened against, concatenated with
  no separator between them. `BamHIKpnI` means BamHI **and** KpnI, which is
  the default pair. Absent if screening ran without enzyme filtering.

`_final_<enzymes>_<bits>`
: The encoded constructs. `<bits>` is the gene's three readout IDs, sorted
  ascending and comma-separated — `4,5,6`. If you regenerate the codebook and
  the gene's bits change, you get a *new* file rather than an overwritten one,
  so stale files with old bit triples can accumulate. Delete them, or the
  assembly step may pick up the wrong one.

### The `ol` value is a gap, not an overlap

This one is genuinely counterintuitive, so read it twice.

The parameter is named "overlap", and larger positive values do mean more
overlap between neighbouring probes. But the default is `-2`, and a **negative
overlap means a gap**: adjacent selected probe windows must be separated by at
least that many nucleotides of untouched transcript. `ol-2` therefore means
"leave at least a 2-nucleotide gap between neighbouring probes", not
"let neighbouring probes overlap by 2".

This is confirmed by the output files themselves. In three test panels
produced with `ol-2`, the smallest gap between adjacent selected windows is 2,
3 and 2 nucleotides respectively — no pair of probes overlaps at all.

Positive `ol` values are only reached when a gene cannot supply enough probes
at the default spacing, at which point the screening step walks upward
(`-2`, then 5, 10, 15, 20) allowing progressively more overlap until it hits
the requested probe count. In the default panel configuration that escalation
is switched off, so `_ol-2_` is what you will normally see.

## The other files in the output directory

Not everything the pipeline writes is a probe table.

`<target>_all.parquet`
: Every alignment of every candidate half against the whole transcriptome —
  the raw material the off-target columns are computed from. This is by far
  the largest file per gene (around 2 MB for a 16 kb transcript, versus 140 kB
  for `_crawled`) and is kept mainly so a re-run can skip realignment. Safe to
  delete once a gene is finished.

`<target>_bowtie.parquet`
: A three-column summary of `_all`: `transcript`, `count`, `name` — how many
  candidate halves hit each transcript in the transcriptome, sorted
  descending. The quickest way to see whether a gene has a homologue problem:
  if a transcript you did not expect sits near the top, your probes are not
  specific.

`<target>_crawled.stats.json`
: Transcript length, how many candidates the crawler produced, the list of
  transcripts treated as acceptable, the top off-target counts, and how many
  candidates survived the match filter.

`<target>_screened_ol*.stats.json`
: How many candidates entered screening, how many survived each filter tier
  (`selected_1` through `selected_7`), and how many pairs came out. If a gene
  is short on probes, this file tells you at which tier it ran out.

`<target>_crawled.coverage.txt`
: A plain-text ASCII picture of where along the transcript the candidates
  fall. Useful for spotting a transcript where everything clusters at one end.

`<target>.log`
: The full per-gene log from the parallel design run, with timestamps. The
  first place to look when one gene fails and the rest succeed.

`<target>_offtarget_counts.csv`, `<target>_acceptable_tss.csv`
: Human/mouse reference datasets only. The first lists significant off-target
  binders for interactive triage; the second lists every transcript that was
  treated as an acceptable binder for this target.

See {doc}`file_formats` for the dataset-level files (`dataset.json`,
`solar_intake.yaml`, index and k-mer files) and for which kind of dataset each
command loads.
