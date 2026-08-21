# What SOLAR is and what these probes do

This page is the orientation for someone who has never run a probe-design
pipeline. It explains the assay `mkprobes` designs for, the vocabulary that
shows up in every output file, and what you physically end up ordering. No
computational background is assumed. Terms in **bold** are also collected in
the {doc}`glossary`.

Everything stated here is drawn from what the code in this repository
actually does. Where the software makes a specific choice (three bits per
gene, 20-nucleotide readouts, 148-nucleotide ordered oligos), that choice is
named explicitly.

## The problem: counting transcripts without losing the map

If you grind up a piece of tissue and sequence it, you learn which genes were
on and roughly how much. You lose *where*. If instead you keep the tissue
intact and light up individual RNA molecules under a microscope, you keep the
map — every detected molecule has an *x, y* (and often *z*) coordinate, so you
can see which cell it was in, which layer, which side of a boundary.

That is what spatial transcriptomics buys you: expression measurements with
the geography still attached. Cell types can be defined *in situ*, gradients
can be measured across a structure, and neighbouring cells can be compared
without dissociating them first.

The catch is throughput. A microscope has a handful of colour channels, but a
panel usually has hundreds of genes. You cannot give each gene its own colour.
The solution used here is to spell out each gene's identity as a *pattern*
across several rounds of imaging — a barcode. Most of what `mkprobes` does is
in service of two goals: put enough well-behaved probes on each target
transcript, and attach the right barcode to every one of them.

## One site, two oligos: the splint and the padlock

A **padlock probe** is a single linear DNA oligo designed so that its two ends
end up next to each other on the target, at which point a ligase can seal the
nick and turn the linear oligo into a covalently closed circle. That "clasp
that snaps shut" behaviour is where the name comes from: it goes on like a
padlock and locks around the target site.

SOLAR uses a **splint/padlock pair** — two separately synthesised oligos per
site, working together:

- The **padlock** carries one target-binding arm (roughly 18–28 nucleotides in
  practice, taken from the 5' side of the chosen window on the transcript),
  followed by the **readout** sequences that encode gene identity.
- The **splint** carries its own target-binding arm, which lands on the
  transcript immediately next to the padlock's arm — the software leaves a gap
  of only a couple of nucleotides between the two arms. The splint's other end
  is a 12-nucleotide clamp (6 nucleotides plus 6 nucleotides) that base-pairs
  with the padlock's 5' end and 3' end at the same time, pulling them
  nose-to-tail so ligase can join them.

This is why every output table has both a `splint` column and a `padlock`
column, why the design step always yields probes in *pairs*, and why you order
**two oligos for every probe**. The package enforces the pairing throughout:
a candidate whose partner half fails any filter is dropped along with it.

The pairing also buys specificity. Circularisation needs the padlock's arm and
the splint's arm to find adjacent sites on the *same* RNA molecule, and needs
the splint's clamp to hold the padlock's own ends together. A partial or
off-target match that satisfies only one of those conditions produces no
circle, and therefore no signal.

## Rolling-circle amplification: turning one circle into one bright dot

A single sealed padlock is a tiny DNA circle, physically threaded around the
site it bound. One circle is far too dim to see. **Rolling-circle
amplification (RCA)** fixes that: a strand-displacing polymerase starts around
the circle and simply keeps going, copying the same small circle over and over
without ever letting go. The product is one very long single strand made of
hundreds of tandem copies of the padlock — including hundreds of copies of its
readout sequences.

Because the polymerase never releases the template, that long strand balls up
right where the original transcript was. Under the microscope it appears as a
single bright, sub-micron dot. These dots go by several names in the
literature — **amplicon**, **rolony**, RCA product — and they are the objects
you actually count. One dot, ideally, means one detected molecule at one
location.

## Readouts and bits: how identity is written onto the amplicon

The dot is bright, but nothing about it says *which gene* it came from. That
information is carried by the **readout** sequences stitched onto the padlock.

A readout is a short, orthogonal DNA sequence — in this package, **20
nucleotides**, drawn from a fixed vendored table of **49 available readouts**,
numbered 1 through 49. They are chosen so that a fluorescent detection
oligo complementary to readout 12 sticks to readout 12 and to nothing else.

Imaging then works in rounds. In each round you wash in a set of fluorescent
detection oligos, photograph the sample, strip them off, and repeat. A dot
lights up in a given round if the amplicon contains the readout that round is
looking for. Over the whole series, each dot produces a pattern of ON and OFF
observations. A single readout position in that pattern is a **bit**.

Two details matter for reading the output files:

- `mkprobes` assigns **three bits per gene** by default. The relevant code
  matrices are all "three-on" designs, and the construct step writes three
  readout IDs into the columns `code1`, `code2`, and `code3`.
- All three readouts go onto **every** padlock for that gene — they are
  stitched in a row onto the same molecule, separated by two-nucleotide
  spacers. It is not the case that some probes carry bit A and others carry
  bit B. Every amplicon from that gene should therefore be positive in all
  three of its rounds.

## The codebook: why combinatorial coding wins

A **codebook** is the lookup table that says which three readout bits belong to
which gene. It is a JSON file mapping each target name to a list of exactly
three integers:

```text
{"Och.687.1": [1, 2, 3], "Och.958.1": [4, 5, 6], "Och.576.10": [7, 8, 9]}
```

The reason to bother with a codebook at all is arithmetic. If each gene needed
its own dedicated readout, N readouts would measure N genes. Choosing *three*
readouts per gene instead gives you "N choose 3" distinguishable patterns:

| Readout bits available | Genes distinguishable (3 bits each) |
| ---: | ---: |
| 10 | 120 |
| 16 | 560 |
| 20 | 1140 |
| 24 | 2024 |
| 30 | 4060 |

Those numbers are the exact row counts of the code matrices vendored inside
the package. The codebook generator picks the smallest code whose capacity
exceeds your gene count by at least 5%, so a 300-gene panel lands on a
16-bit code (560 codewords) and a 1000-gene panel on a 20-bit code.

Two more things the generator does, worth knowing because they show up in the
file it writes:

- **Bit ordering is deliberate.** A gene's three bits are spread across the
  readout-ID space rather than being allowed to land together, and a small
  hard-coded set of triplets that the code comments describe as perfectly
  confounding the imaging rounds is swapped out onto unused codewords.
- **Assignment can be expression-informed.** If you supply a table of
  expression values, the generator tries many assignments and keeps the one
  that spreads total expression most evenly across the readout bits, so no
  single bit is dominated by a handful of very highly expressed genes. This is
  entirely optional; without it, assignment is a plain seeded shuffle.

## Error-correcting distance, MHD4, and Blank codes

Imaging is imperfect. A round can be too dim, a dot can be missed, a
neighbouring dot can bleed in. Each of those flips one bit of a barcode. What
happens next depends on how far apart the codewords are.

**Hamming distance** is just the number of bit positions in which two barcodes
differ. The **minimum Hamming distance (MHD)** of a codebook is the smallest
such distance over all pairs of codewords in it, and it determines what you can
recover from:

- **MHD2** — every pair of valid barcodes differs in at least 2 positions. A
  single flipped bit lands you on something that is not a valid barcode, so you
  can *detect* that something went wrong, but you cannot tell which valid
  barcode it came from.
- **MHD4** — every pair differs in at least 4 positions. Now a single flipped
  bit leaves you closer to the true barcode than to any other, so it can be
  *corrected*, and a two-bit error can still be detected. The cost is capacity:
  spacing codewords further apart means fewer of them fit, so you need more
  imaging rounds for the same number of genes.

Be precise about what this package produces. The code matrices that
`make-codebook` actually selects from are the files named
`<N>bit_on3_dist2.csv`: every possible way of choosing 3 bits out of N, with a
minimum Hamming distance of 2. That is not an error-*correcting* code. What it
does give you is the **constant-weight rule**: every valid barcode has exactly
three bits ON, so any dot decoding to two or four ON bits is visibly an error
and can be discarded. MHD4-style matrices are also vendored in the package
directory, but the codebook generator does not pick them up, so a codebook
produced by the normal workflow is a 3-of-N distance-2 code. Plan your
analysis accordingly.

**Blank codes** are the safety net. After every real target has been assigned a
codeword, the remaining codewords in the matrix are handed out as `Blank-1`,
`Blank-2`, and so on, and written into the same JSON file. No probes are ever
synthesised for them — the design step explicitly skips every `Blank-*` entry.
So a blank barcode can only ever appear in your images as a mistake: optical
noise, misassignment, or a genuinely spurious amplicon. Counting how often
blanks show up, relative to how many blank codewords exist, is your empirical
false-positive rate. This is the single most useful QC number a panel gives
you, which is why the generator warns when fewer than 5% of the coding
capacity is left blank.

## How it fits together

```text
  1. Pick a window on the target transcript and split it into two arms

     target RNA  5'---------[  arm 1  ]--[  arm 2  ]---------------3'
                                                (a 0-2 nt gap between arms)

  2. Two oligos bind that window. The padlock covers one arm and trails
     its three readouts; the splint covers the other arm AND clamps the
     padlock's two ends together (6 nt + 6 nt).

     padlock   [ arm 1' ][ readout ][ readout ][ readout ]
                ^ 5' end                              3' end ^
                 \                                          /
     splint       [ 6 nt clamp ][ 6 nt clamp ][   arm 2'   ]

  3. Ligase seals the padlock's ends -> a closed DNA circle, locked in
     place around the transcript. No adjacent pair, no circle.

                        .-------------------.
                       (   circular padlock  )
                        '-------------------'

  4. Rolling-circle amplification copies the circle hundreds of times
     in place -> one bright dot (amplicon / "rolony") per molecule.

              circle ->  [copy][copy][copy][copy][copy][copy]...

  5. Successive imaging rounds detect the readouts. The ON/OFF pattern
     is the barcode; the codebook turns it back into a gene name.

              readout 4 : ON        readout 4,5,6 -> Och.958.1
              readout 5 : ON
              readout 6 : ON
              all others: off
```

## What you physically order at the end

The end product of the workflow is an **oligo pool**: a single tube of many
thousands of synthesised DNA sequences, ordered from a commercial
array-synthesis vendor as a plain list of sequences, one per line.

From a real assembled test panel in this repository (3 genes, 16 probe pairs
each), the output file contains 96 sequences — two per probe pair, splint and
padlock interleaved — and every one of them is **148 nucleotides** long. That
length is not all probe. Each ordered oligo is the working probe wrapped in
constant machinery:

- constant primer sequences at both ends, so the whole pool can be amplified
  from the tiny amounts an array synthesiser delivers;
- restriction sites (KpnI and BamHI in the default configuration) that let you
  cut the amplified product down to the working probe afterwards;
- filler sequence to bring everything to a uniform length.

This is also why probe design filters out any candidate containing a KpnI or
BamHI site: a site inside the probe body would be cut along with the intended
ones and destroy the probe. Restriction-site screening is on by default for
exactly this reason.

Alongside the pool you get a provenance record — the package version, the
codebook and its hash, the number of probe pairs, and the settings used — so a
pool that arrives six months later can still be traced back to the design that
produced it.

## Where to go next

- Practical prerequisites, reference files, and machine requirements:
  {doc}`before_you_start`
- Every term used above, defined once: {doc}`glossary`
- The columns in the parquet files the pipeline writes:
  {doc}`reference/columns`
- The commands, in order: {doc}`getting_started`
- Stage-by-stage internals, written for developers:
  {doc}`workflows/assay_and_under_the_hood`
