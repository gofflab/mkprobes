# Glossary

Two kinds of vocabulary appear in `mkprobes` output, and both are collected
here alphabetically.

The first kind is established field vocabulary — *padlock*, *Tm*, *k-mer*,
*GTF*. These words are not going away; they are what the assay and the
surrounding tools are called, and they appear in column names, log lines, and
error messages. The definitions below are the minimum you need to read those
messages.

The second kind is this tool's own shorthand — *crawled*, *screened*,
*construct*, *bcidx*. You will see these in filenames whether or not anyone
ever explains them, so they are explained here.

For the biology behind the first group, start with {doc}`what_is_solar`.

bcidx
: A small integer in a panel's `manifest.json` that selects which pair of
  constant primer/adapter sequences the assembled oligos get. The package
  ships a table of 201 header/footer rows; a probeset with `bcidx: 0` uses
  rows 0 and 1 (splint and padlock respectively), `bcidx: 1` uses rows 2 and
  3, and so on. Give two panels different `bcidx` values if you want to
  amplify them independently from the same synthesis batch.

bit
: One ON/OFF position in a gene's barcode — in practice, one **readout**
  sequence, detected in one imaging round. `mkprobes` assigns three bits per
  gene. A codebook entry such as `[4, 5, 6]` means "bits 4, 5 and 6 are ON for
  this gene, every other bit is OFF".

blocklist
: A set of sequences that probes must not resemble, supplied as FASTA files
  (typically ribosomal and transfer RNA) when a dataset is built. Their
  15-nucleotide **k-mers** are counted into a file named `blocklist15.jf`, and
  any candidate probe sharing even one 15-mer with that set is thrown out.
  Optional but strongly recommended — see {doc}`before_you_start`. When a
  dataset was built without one, the pipeline logs `blocklist filtering is a
  no-op` on every run.

bowtie2 index
: A pre-built, binary search structure that lets the aligner `bowtie2` find
  where a short sequence matches a transcriptome, fast. It is a set of files
  ending in `.bt2` sitting next to the FASTA they were built from, created
  once when the dataset is built. Note the **2**: `bowtie` and `bowtie2` are
  different programs with incompatible indexes, and this pipeline uses
  `bowtie2` exclusively.

canonical transcript
: The single transcript a database designates as a gene's default
  representative, when a gene has several. Reference (human/mouse) datasets
  can resolve gene names through Ensembl's canonical annotation; for custom
  species there is no such designation, so the workflow instead offers
  "longest" or "all" as the selection rule.

codebook
: The JSON file mapping each target to its three readout **bit** IDs, plus
  `Blank-*` entries for every codeword left unused. It is the panel's master
  list: the design driver treats it as the work list of genes to process
  (skipping blanks), and the final assembly step reads it back to know which
  file belongs to which gene. Changing it invalidates everything downstream,
  which is why the package prints a short hash of it.

construct
: The stage that takes screened probe pairs and glues the gene's three readout
  sequences onto the padlock, producing the `*_final_*.parquet` files. Also
  the name of the resulting object — a "construct" is a probe sequence with
  its encoding attached. Constructs that would end up containing a run of five
  identical bases or a BamHI site are silently discarded, which is why a gene
  usually has slightly fewer final probes than screened ones.

crawled
: This tool's word for the raw candidate stage. The "crawler" walks along the
  target transcript one starting position at a time, growing each window until
  it satisfies length, GC, melting-temperature and hairpin limits, and emits
  every window that passes. `*_crawled.parquet` therefore holds thousands of
  overlapping, not-yet-selected candidates — the search space, not the answer.

FASTA
: The plainest possible sequence file format: a `>` line naming a sequence,
  then the sequence itself, repeated. Used here for genomes (one record per
  chromosome or scaffold), transcriptomes (one record per transcript), and
  blocklist inputs. It contains sequence only — no coordinates, no gene
  structure, no annotation.

FPKM / TPM
: Units of RNA abundance from bulk RNA-seq, both normalising for sequencing
  depth and transcript length so that numbers are comparable across genes.
  `mkprobes` never requires them. They are optional inputs in two places: to
  balance expression load across readout bits when generating a codebook, and
  to down-weight highly expressed off-target transcripts during screening. Any
  reasonable per-gene abundance estimate works; the exact unit does not matter
  because only the relative values are used.

gene, transcript, isoform
: A **gene** is a named locus. A **transcript** is one specific RNA molecule
  the locus produces, with a definite sequence and its own identifier. An
  **isoform** is one of several alternative transcripts of the same gene.
  This distinction is load-bearing here: probes are designed against a
  *transcript* sequence, not a gene, so a gene with five isoforms forces a
  choice about which one to target. Probes hitting the target's sibling
  isoforms are treated as acceptable rather than as off-target, since they are
  still the right gene.

GTF / GFF3
: Two closely related plain-text formats describing where genes live on a
  genome: one line per feature (gene, transcript, exon), giving the sequence
  name, start and end coordinates, strand, and a final column of attributes.
  A genome FASTA tells you the letters; the GTF tells you which stretches of
  those letters are exons of which transcript. Both files must describe the
  *same* assembly. GFF3 uses `key=value` attributes and GTF uses
  `key "value";`; `mkprobes` detects which one you gave it and converts GFF3
  to GTF with `gffread` automatically.

hairpin
: A probe folding back and base-pairing with itself instead of with its
  target. Candidates are scored for this, and the score (column `hp`) is the
  predicted melting temperature in °C of the most stable self-structure the
  sequence can form — higher means a more stable, more troublesome hairpin.
  Note the column is `hp`, not `homopolymer`; the two are unrelated.

homopolymer
: A run of the same base repeated — `AAAA`, `GGGG`. Long runs are hard to
  synthesise accurately, prone to slippage during amplification, and can form
  unwanted structures. Candidates containing a run of four or more identical
  bases fail the `ok_homopolymer` check, and the crawler additionally refuses
  windows containing `AAAAAA`, `TTTTTT`, `CCCCC` or `GGGGG` outright.

k-mer
: Every substring of a given length k, taken from a sequence. `GATTACA` has
  five 3-mers. K-mers make "does this probe share any exact 18-letter stretch
  with anything unwanted?" answerable by table lookup instead of by alignment,
  which is why the pipeline counts them (with `jellyfish`) into `.jf` files
  when it builds a dataset. Two k-mer sizes are used: 18-mers for
  repetitive-region screening, and 15-mers for the rRNA/tRNA **blocklist**.

manifest
: Currently used for **three different files**. Read the sentence around it
  before assuming which one is meant.

  1. `manifest.json` — the panel manifest, a list of **probeset** entries, and
     the file you hand to the assembly step.
  2. `dataset.json` — sometimes called the dataset manifest; it records what a
     custom dataset contains (FASTA, GTF, index and k-mer file names,
     blocklist, registered annotation tables).
  3. `solar_intake.yaml` — the provenance manifest written when a genome and
     annotation are ingested, holding input checksums, tool versions, the
     literal command that was run, and QC counts, plus stub fields for the
     operator to complete.

MHD4 / MHD2
: The minimum Hamming distance of a codebook — the smallest number of bit
  positions in which any two valid barcodes differ. MHD2 lets you *detect* a
  single misread bit; MHD4 lets you *correct* one and still detect two. The
  code matrices `mkprobes` actually generates codebooks from are 3-of-N
  distance-2 designs (files named `<N>bit_on3_dist2.csv`), so what you get is
  detection via the constant-weight rule — exactly three bits ON — rather than
  correction. See {doc}`what_is_solar` for the full explanation.

off-target
: A place a probe binds that is not the transcript it was designed for. The
  pipeline does not treat this as a yes/no property. Every candidate is
  aligned against the whole transcriptome, and each unwanted hit is quantified
  three ways: how many bases matched in total (`match`), the longest
  uninterrupted matching stretch (`match_consec`), and the predicted melting
  temperature of the best such stretch (`max_tm_offtarget`). Sibling isoforms
  of the target gene, and any transcripts you explicitly allow, are excluded
  from the off-target accounting.

oligo pool
: The physical deliverable: one tube containing many thousands of distinct
  synthesised DNA sequences, made in parallel on an array and ordered from a
  commercial vendor as a plain list of sequences. Because per-sequence yield is
  tiny, pools are amplified before use, which is why each ordered oligo carries
  constant primer sequences and restriction sites in addition to the probe
  itself.

padlock
: The oligo that gets circularised. It carries one target-binding arm plus the
  gene's readout sequences; when its two ends are brought together by the
  **splint**, a ligase seals it into a closed circle that is then copied by
  **rolling-circle amplification**. Named for the way it clasps shut around
  the site. In the output tables, the `padlock` column holds only the
  target-binding half in target orientation; the readouts are added later, at
  the construct stage.

panel
: The set of genes designed and ordered together — your gene list, its
  codebook, and everything produced from them. Panel-level questions (does
  every gene have enough probes? how much coding capacity is left blank?) are
  asked once across all genes, not per gene.

probeset
: A single entry in a panel `manifest.json`: a name, a species, the codebook it
  uses, its **bcidx**, and how many probes per gene to keep. One manifest can
  hold several probesets, which are assembled together into one pool. Not to
  be confused with "probe set" in the loose sense of "the probes for one
  gene".

pseudogene
: A copy of a gene that has lost its function but retains recognisable
  sequence similarity. Pseudogenes are a classic false-positive source: a probe
  can be perfectly specific to your gene's *sequence* and still bind a
  pseudogene, and you will never see the difference in an image. For human and
  mouse the pipeline explicitly tallies pseudogene hits and lets you accept or
  reject them per gene (column `maps_to_pseudo`); for custom species there is
  no pseudogene annotation to work from, so this screening does not happen.

readout
: A short, orthogonal DNA sequence appended to the padlock, which a
  complementary fluorescent oligo detects during one imaging round. In this
  package readouts are 20 nucleotides long and come from a fixed vendored
  table of 49, numbered 1–49. A gene's three readout IDs are what a codebook
  entry lists, and what the `code1`/`code2`/`code3` columns record.

restriction site
: A short sequence a restriction enzyme recognises and cuts. Assembly uses
  KpnI and BamHI sites deliberately, to trim the constant amplification
  sequences off the ordered oligos. That only works if the probe body contains
  no such site of its own, so candidates containing one are filtered out
  during screening — which is why enzyme names appear in output filenames
  (`_BamHIKpnI`).

rolling-circle amplification (RCA)
: Copying a small circular DNA template continuously with a strand-displacing
  polymerase, which never lets go and so produces one very long strand of
  tandem repeats. Because the product stays attached where it started, it
  collapses into a single bright, localised ball of DNA (an amplicon, often
  called a "rolony") that can be imaged as one dot. This is the amplification
  step that makes a single circularised padlock visible.

SAM / CIGAR
: SAM is the standard text format an aligner emits, one line per alignment,
  reporting where a sequence matched, how well, and with what differences.
  CIGAR is one field of that line, a compact description of the alignment
  shape: `23M` means 23 aligned bases; `1S21M4S` means 1 base clipped off the
  front, 21 aligned, 4 clipped off the end. Several output columns
  (`flag`, `cigar`, `pos`, `aln_score`, `n_mismatches`, `edit_distance`,
  `mismatched_reference`) are simply SAM fields carried through, and are
  described in {doc}`reference/columns`.

screened
: This tool's word for the selection stage. From the thousands of overlapping
  crawled candidates, screening applies a cascade of quality filters and then
  picks a non-overlapping subset tiled along the transcript, keeping only
  pairs where both halves survive. `*_screened_ol*.parquet` holds tens of
  probe pairs, not thousands of candidates.

splint
: The partner oligo to the padlock. It binds the transcript immediately
  adjacent to the padlock's arm, and its other end is a 12-nucleotide clamp
  (6 + 6) that grabs the padlock's 5' and 3' ends simultaneously and holds
  them together for ligation. It is a separate synthesised oligo, so every
  probe costs two sequences in the pool.

Tm (melting temperature)
: The temperature at which half of a duplex has come apart — the standard
  single number for "how tightly does this stick?". Higher Tm means a more
  stable duplex. It depends on length, GC content, salt, and, importantly
  here, formamide, which destabilises duplexes and so lowers Tm; the pipeline
  applies a formamide correction, so the `tm` values you see are lower than an
  uncorrected calculator would report. Tm is used in three distinct places:
  to size candidate windows, to split a window into two balanced arms, and to
  judge how dangerous an off-target hit is (`max_tm_offtarget`).
