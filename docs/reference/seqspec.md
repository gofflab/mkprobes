# The SOLAR seqspec

[seqspec](https://github.com/pachterlab/seqspec) is a machine-readable format
for describing what a library's molecules are made of, region by region.
`mkprobes` ships them for both halves of the assay: the oligo pool you order,
and the readout probes you image with.

They are worth having for three reasons. A pool can be checked against the spec
before it goes to a vendor, which is the last point at which a mistake is free.
The spec is a precise, shareable answer to "what exactly did you order?" — more
useful in a methods section than prose. And it draws, which is the fastest way
to explain the construct to someone new.

| File | What it describes |
| --- | --- |
| `solar_bcidx<N>.seqspec.yaml` | the splint/padlock oligo pool from `mkprobes assemble`, for `bcidx` N |
| `solar_readout_probes.seqspec.yaml` | the 49 readout (detection) probes |
| `solar_readouts.txt` | the 49 readout sequences, referenced by every spec |

They live in the installed package, under `mkprobes/data/seqspec/`. All are
generated from the same tables the pipeline designs against — the header/footer
table and the readout table — so they cannot drift from what the code builds.

There is **one pool spec per `bcidx`**, one for each value the manifest field
accepts. The readout probes do not depend on the index, so there is one of
those. Which index a pool used is worked out from the pool itself, so you
rarely need to name it.

## Validate a pool

```bash
mkprobes validate-pool panel_a/generated/panel_a_final.txt
```

Pass the codebook too, which enables the strongest check of the set:

```bash
mkprobes validate-pool panel_a/generated/panel_a_final.txt -c panel_a/codebook.json
```

A clean pool reports one line. A bad one lists what is wrong and exits
non-zero, naming the probe pair so you can find it.

What it checks:

- Every oligo is exactly **148 nt** and carries the right regions in the right
  order — both primer binding sites, the restriction scars, the homology arm,
  and (on padlocks) three readouts from the vendored table of 49.
- The working probe contains **no KpnI or BamHI site**. Those sites belong at
  the handle boundaries; one inside the probe means the digest that releases it
  also destroys it.
- Splints and padlocks **alternate**, and each splint's 6+6-nt clamp actually
  templates its partner's two ends. A pair that fails this makes no circle,
  so it produces no signal at all.
- With `--codebook`: every codeword in the pool is one the codebook assigns,
  and **no `Blank-*` codeword was synthesised**. Blanks are how you measure
  your false-positive rate, so a blank with real probes against it quietly
  destroys that estimate — this check exists to catch that before you order.

This complements `mkprobes assemble`, which asserts much of the same geometry
as it builds. The difference is that `validate-pool` runs on the file itself,
so it still works on a pool from six months ago, from a colleague, or one that
has been through a spreadsheet.

## Draw the construct

```bash
mkprobes draw-spec --which pool -o solar_pool.png
mkprobes draw-spec --which readout -o solar_readout.png
```

This needs the optional `viz` extra:

```bash
pip install 'mkprobes[viz]'
```

The figure is a **generalized view**: every region is labelled with the length
range it is allowed, not the length it happens to have in one probe. Regions
whose lengths trade off against each other are drawn at a representative split,
so the picture adds up to the 148 nt actually synthesised.

Reading it: grey is a primer binding site, red is a restriction site, amber is
the single base that site leaves behind on the probe, teal is the
target-homology arm, purple is a readout, and pale grey-blue is a linker,
spacer or filler. The splint's long `splint_backfill` and all four primer sites
are there to make synthesis and amplification work — they are cut away before
the probe is used.

Arrows run 5'→3'. Two kinds of region break that:

- The two `*_primer_3p` regions **point the other way**, because the reverse
  primer is the reverse complement of that region and extends back toward the
  5' end.
- The four restriction sites are drawn **square-ended, with no arrow at all**.
  `GGATCC` and `GGTACC` are palindromes — each is its own reverse complement,
  so it reads identically on both strands — and the enzyme cuts the duplex
  rather than one strand. An arrow would assert a direction the feature does
  not have.

seqspec has no field for strand, so this is one thing the drawing says that the
YAML can only say in prose.

To draw a different index:

```bash
mkprobes draw-spec --bcidx 7 -o solar_bcidx7.png
```

## Three things the spec cannot say

seqspec describes one linear molecule as a tree of regions. SOLAR strains that
in three places, so these are checked by `validate-pool` instead and are
written into the comment header of the file itself:

1. **A pool entry is a pair.** `splint_oligo` and `padlock_oligo` are two
   separately synthesised molecules, not one 296-nt strand. seqspec has no way
   to say "these siblings are separate molecules", so they are nested under one
   `dna` region and the pairing rule is enforced by the validator.
2. **Variable regions co-vary.** On the splint, pad + arm is always 33 nt; on
   the padlock, arm + filler is always 27 nt. seqspec sums each child's range
   independently, so it advertises the splint as 137–159 nt when every real one
   is 148.
3. **The readouts must be the gene's codeword.** Any three of the 49 satisfy
   the spec; only the codebook says which three are right, which is what
   `--codebook` is for.

One more detail worth knowing when reading a padlock: the three readouts appear
in a **different order on different probes**. The construct step cycles through
permutations of a gene's three bits, so `padlock_readout_1` means "the readout
in the first slot", not `code1`. Only the set of three is fixed per gene, which
is why validation compares them unordered.

## The restriction sites, and where the scar lands

Each oligo carries one KpnI site at the 5' end and one BamHI site at the 3'
end. Both appear as their own regions in the spec, and both are **positioned so
that the cut falls inside the site**:

```text
     KpnI  GGTAC^C                         BamHI  G^GATCC
   ...primer GGTAC | C spacer ...probe... clamp G | GATCC primer...
             \___/   ^                          ^   \___/
          discarded  the only base kept   kept       discarded
```

Both sites read `GGATCC` / `GGTACC` in the 5'→3' direction of the oligo as
written, exactly once each — and being palindromic, they read the same on the
complementary strand, which is what lets a double-strand cutter work on them
regardless of how the duplex is presented.

The asymmetry of the cut is the design, not an accident. Five of each site's six bases
sit on the handle that the digest throws away; exactly **one** base survives on
the working probe — a `C` at the 5' end, a `G` at the 3' end. So the scar the
enzyme leaves is deliberately spent on the part nobody cares about, and the
probe itself picks up a single base rather than a six-base footprint.

The regions that carry this:

| Region | What it is |
| --- | --- |
| `*_kpni_site` | `GGTAC` — five of the six bases of `GGTAC^C`, discarded with the 5' handle |
| `*_kpni_retained` | `C` — the one base of that site left on the working probe |
| `*_bamhi_site` | `GATCC` — five of the six bases of `G^GATCC`, discarded with the 3' handle |

The retained `G` at the 3' end is not broken out as its own region. On both
oligos it is the last of the six nucleotides that clamp the padlock for
ligation, and those six are the functional unit — splitting one base off would
make the spec tidier and the chemistry harder to read. Its identity is stated
in `splint_clamp_3p` and `padlock_ligation_end` instead.

Because these are the enzyme's sequences rather than the panel's, they are
**identical at every `bcidx`** — the one part of the flanking machinery that
does not move. A test asserts exactly that, and separately reassembles each
site across its cut to confirm the six bases still spell `GGTACC` and `GGATCC`.

## What `bcidx` actually controls

`bcidx` picks the two primer pairs that amplify a probe set, so several panels
can be pooled and each still amplified on its own. It reaches further than the
primers, though, and this is the reason there is a whole spec per index rather
than a footnote. **Ten regions** change with it:

- both primer binding sites on each oligo (four regions),
- the design bases flanking each cut — the spacer after the KpnI site, and the
  six nucleotides at the padlock's 3' end,
- the splint's **clamp**, which has to template the padlock's actual ends.

The restriction sites themselves do not: they are enzyme sequence, and are the
same at every index.

So a pool checked against the wrong index does not merely mismatch at the
primers; it mismatches at the ligation junction, the part that decides whether
a probe circularises at all. Everything else — arms, readouts, spacers,
backfill — is identical across indices.

`validate-pool` detects the index from the pool and says which one it found.
Pass `--bcidx` to assert a particular one instead, which is the way to check
that a panel is the index you meant it to be:

```bash
mkprobes validate-pool panel_a/generated/panel_a_final.txt --bcidx 7
```

## Regenerating

The specs are generated, not hand-written:

```bash
python scripts/generate_seqspec.py
```

A test re-runs the generator and compares every file, so editing the
header/footer table or the readout table without regenerating is caught in CI
rather than at the vendor. Adding rows to the header/footer table raises the
`bcidx` ceiling, and re-running the generator writes the new specs (and deletes
any whose index no longer exists).

## Checking the specs themselves

To verify the files still conform to the seqspec specification (a developer
task — validating a pool does not need it):

```bash
seqspec check src/mkprobes/data/seqspec/solar_bcidx0.seqspec.yaml
```

Expect exactly one error, on every one of these files:

```text
[error 1] Reads must have the same number of files
```

That is not a defect. seqspec compares the file counts of sequencing reads, and
its check cannot pass when there are no reads at all. SOLAR is an imaging
assay: nothing is sequenced, so `sequence_spec` is empty by design. The test
suite allows this one error and no other.
