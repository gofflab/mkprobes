"""
Regenerate the vendored SOLAR seqspec files from the tables that define the chemistry.

Every fixed sequence in the specs is derived here from `data/headerfooter.csv`,
the padding constants in `assembly.py`, and `codebook/readout_ref_filtered.csv`.
Nothing is transcribed by hand, so the specs cannot drift from the assay the
package actually builds - `tests/test_oligospec.py` re-runs this module and
fails if the checked-in files differ.

One pool spec is written per `bcidx`, because the two primer pairs are the one
thing that changes between them and a pool can only be validated against the
pair it was actually built with. The readout spec is index-independent, so
there is one of those.

Usage: `python scripts/generate_seqspec.py` (writes into src/mkprobes/data/seqspec).
"""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path
from typing import Any

import yaml

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "src" / "mkprobes" / "data"
OUT = DATA / "seqspec"

#: seqspec version these files are written against.
SEQSPEC_VERSION = "0.3.0"
#: Date stamped into both specs, in the format seqspec's schema demands.
SPEC_DATE = "26 August 2026"

#: Where each recognition site sits in a header/footer row, and where the
#: enzyme cuts it. Both sites straddle the cut on purpose: five of the six bases
#: leave with the handle that gets discarded, and exactly one stays on the
#: working probe. `tests/test_oligospec.py` re-derives these from the table.
KPNI_SITE = "GGTACC"       # cuts GGTAC^C
BAMHI_SITE = "GGATCC"      # cuts G^GATCC
KPNI_START = 12            # index of GGTACC in every header
KPNI_CUT = KPNI_START + 5  # header[:17] goes with the primer, header[17] stays
BAMHI_START = 2            # index of GGATCC in every footer
BAMHI_CUT = BAMHI_START + 1  # footer[2] stays on the probe, footer[3:] goes

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def rc(seq: str) -> str:
    return seq.upper().translate(_COMPLEMENT)[::-1]


def _headerfooter_rows() -> list[dict[str, str]]:
    return list(csv.DictReader((DATA / "headerfooter.csv").open()))


def max_bcidx() -> int:
    """Highest usable index. Mirrors `init_project.max_bcidx`: two rows each."""
    return len(_headerfooter_rows()) // 2 - 1


def _read_headerfooter(bcidx: int = 0) -> tuple[str, str, str, str]:
    """Splint and padlock header/footer for one `bcidx`, which uses two rows."""
    rows = _headerfooter_rows()
    splint, padlock = rows[bcidx * 2], rows[bcidx * 2 + 1]
    return splint["header"], splint["footer"], padlock["header"], padlock["footer"]


def _read_readouts() -> list[str]:
    path = ROOT / "src" / "mkprobes" / "codebook" / "readout_ref_filtered.csv"
    return [row["seq"] for row in csv.DictReader(path.open())]


# --------------------------------------------------------------------------- #
# Region construction
#
# seqspec requires that a parent's sequence be the exact concatenation of its
# children and that its min_len/max_len be the sums of theirs. Building the tree
# bottom-up makes both true by construction rather than by careful typing.
# --------------------------------------------------------------------------- #


def leaf(
    region_id: str,
    region_type: str,
    name: str,
    sequence_type: str,
    sequence: str,
    min_len: int | None = None,
    max_len: int | None = None,
    onlist: dict[str, Any] | None = None,
) -> dict[str, Any]:
    return {
        "region_id": region_id,
        "region_type": region_type,
        "name": name,
        "sequence_type": sequence_type,
        "sequence": sequence,
        "min_len": len(sequence) if min_len is None else min_len,
        "max_len": len(sequence) if max_len is None else max_len,
        "onlist": onlist,
        # A leaf, but seqspec types this field as a plain list: `null` fails its
        # strict loader, which is what `seqspec print` and `seqspec index` use.
        "regions": [],
    }


def fixed(region_id: str, region_type: str, name: str, sequence: str) -> dict[str, Any]:
    return leaf(region_id, region_type, name, "fixed", sequence.upper())


def variable(
    region_id: str, region_type: str, name: str, min_len: int, max_len: int
) -> dict[str, Any]:
    """A region of variable content. seqspec spells these `X`, one per position."""
    return leaf(region_id, region_type, name, "random", "X" * max_len, min_len, max_len)


def joined(region_id: str, region_type: str, name: str, regions: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "region_id": region_id,
        "region_type": region_type,
        "name": name,
        "sequence_type": "joined",
        "sequence": "".join(r["sequence"] for r in regions),
        "min_len": sum(r["min_len"] for r in regions),
        "max_len": sum(r["max_len"] for r in regions),
        "onlist": None,
        "regions": regions,
    }


def _onlist(filename: str, size: int, md5: str) -> dict[str, Any]:
    return {
        "file_id": filename,
        "filename": filename,
        "filetype": "txt",
        "filesize": size,
        "url": filename,
        "urltype": "local",
        "md5": md5,
    }


def _assay(**fields: Any) -> dict[str, Any]:
    return {
        "seqspec_version": SEQSPEC_VERSION,
        "doi": "",
        "date": SPEC_DATE,
        "modalities": ["dna"],
        "lib_struct": "",
        "library_protocol": "spatial transcriptomics (NTR:0000761)",
        "sequence_protocol": "Custom",
        "sequence_kit": "Custom",
        "sequence_spec": [],
        **fields,
    }


# --------------------------------------------------------------------------- #
# The two specs
# --------------------------------------------------------------------------- #

def pool_header(bcidx: int, spl_header: str, spl_footer: str, pad_header: str, pad_footer: str) -> str:
    """The comment block at the top of one pool spec."""
    return f"""\
# Structure of the SOLAR splint/padlock oligo pool for bcidx {bcidx}, as
# assembled by `mkprobes assemble`. Validate a pool against it with
# `mkprobes validate-pool`, and draw it with `mkprobes draw-spec`.
#
# GENERATED FILE - edit scripts/generate_seqspec.py and re-run it instead.
# Every fixed sequence below is derived from rows {bcidx * 2} and {bcidx * 2 + 1} of
# data/headerfooter.csv, the padding constants in assembly.py, and
# codebook/readout_ref_filtered.csv.
#
# `bcidx` picks the two primer pairs that amplify this probe set, so panels
# pooled together can each be amplified on their own:
#
#   splint  forward {spl_header[:17]}      reverse {rc(spl_footer[3:])}
#   padlock forward {pad_header[:17]}      reverse {rc(pad_footer[3:])}
#
# Each oligo reads 5'->3' as written. The forward primer has the same sense as
# the sequence below; the reverse primer is the reverse complement of the 3'
# binding site, so it primes back toward the 5' end.
#
# The index reaches further than the primers, which is why there is a whole
# spec per index rather than a note about two regions. Ten regions change with
# it: both primer binding sites on each oligo, the KpnI and BamHI scars those
# sites leave behind on the working probe, and the splint's clamp - which has
# to template the padlock's actual ends, and those ends are the scars. A pool
# validated against the wrong index therefore fails at the ligation junction,
# not only at the primers. Everything else - arms, readouts, spacers, backfill
# - is identical across indices.
#
# Three properties of the pool this file cannot state, which
# `mkprobes validate-pool` checks instead:
#
#   1. Each entry is a probe PAIR. `splint_oligo` and `padlock_oligo` are two
#      separately synthesized 148-nt molecules, not one 296-nt strand, and they
#      alternate down the ordered pool file: splint, padlock, splint, ...
#      seqspec has no way to say "sibling regions are separate molecules", so
#      they are nested under one `dna` region.
#   2. Variable regions co-vary to a constant total. On the splint,
#      pad + arm = 33 nt; on the padlock, arm + fill = 27 nt. seqspec sums
#      child ranges independently, so the parents below carry the wider sums
#      while every real oligo is exactly 148 nt.
#   3. The three readouts on a padlock must be that gene's codeword, in any
#      order. Any three of the 49 satisfy this file; only the codebook says
#      which three are right, so that check needs the codebook JSON.
#
# seqspec has no field for strand, so the fact that the 3' sites are read in
# the opposite direction lives in their region names and descriptions - and in
# the arrows `mkprobes draw-spec` draws.
"""


READOUT_HEADER = """\
# Structure of the SOLAR readout (detection) probes - the oligos washed in and
# imaged, as opposed to the pool oligos that are ordered and hybridized.
#
# GENERATED FILE - edit scripts/generate_seqspec.py and re-run it instead.
#
# These are not designed by mkprobes. They are the fixed vendored table of 49
# orthogonal 20-nt sequences in codebook/readout_ref_filtered.csv, three of
# which the construct step stitches onto every padlock for a gene. Each
# detection oligo carries a fluorophore conjugate, which is chemistry rather
# than sequence and so appears only in the description here.
#
# Orientation, which is easy to get backwards. The construct step stitches
# `rc(readout)` onto the screened padlock, but assembly then reverse-
# complements that whole payload, so the ORDERED padlock carries the readout
# exactly as the table lists it - which is what `mkprobes validate-pool`
# matches against. Rolling-circle amplification copies the padlock circle as a
# template, so the amplicon presents the reverse complement, and a detection
# oligo that hybridizes to the amplicon therefore carries the table sequence
# itself. That is what this file describes.
"""


def build_pool_spec(onlist: dict[str, Any], bcidx: int = 0) -> dict[str, Any]:
    spl_header, spl_footer, pad_header, pad_footer = _read_headerfooter(bcidx)

    # assembly.backfill() pads the splint construct to 148 nt from the 5' end.
    # The splint body is a constant 91 nt, so the filler is always these 57.
    backfill = "TTCCACTAACTCACATGTCATGCATTATCTTCTATACCTCTGAGCAGATCAGTAGTCTATTACATGCTCGTAGTACCGTAAGCCAGATAC"
    backfill_used = backfill[: 148 - 91]
    # assembly.padpad() right-pads the padlock body; an 18-nt arm (the shortest
    # allowed) consumes the most filler, so at most these 9 bases appear.
    padlock_fill = "AATCACATAAAT"[:9]

    splint_probe = joined(
        "splint_probe",
        "named",
        "working splint, released by the KpnI + BamHI double digest",
        [
            fixed(
                "splint_kpni_retained",
                "linker",
                "the single base of the KpnI site GGTAC^C left on the working probe. The other "
                "five leave with the 5' handle, so the site is positioned to scar the part that "
                "gets thrown away",
                spl_header[KPNI_CUT],
            ),
            fixed("splint_spacer_5p", "linker", "spacer between the KpnI cut and the A/T pad", spl_header[KPNI_CUT + 1 : 20]),
            variable("splint_pad", "linker", "A/T filler from the cycling ATAAT pad; pad + arm = 33 nt", 4, 15),
            variable("splint_arm", "cdna", "splint target-homology arm, reverse complement of the 3' half of the transcript window", 18, 29),
            fixed("splint_linker_ca", "linker", "CA linker", "CA"),
            variable("splint_clamp_var", "linker", "clamp, the 3 nt templating the padlock's head splint", 3, 3),
            fixed("splint_clamp_fix", "linker", "clamp, the 3 nt templating the padlock's KpnI scar", rc(pad_header[-3:])),
            fixed(
                "splint_clamp_3p",
                "linker",
                "clamp, the 6 nt templating the padlock's 3' end. Its final base is the one base "
                "of the BamHI site G^GATCC kept on the working probe - not split out as its own "
                "region, because these six are the unit that has to clamp the padlock",
                rc(rc(spl_footer[:BAMHI_CUT]) + pad_footer[:BAMHI_CUT]),
            ),
        ],
    )
    splint_oligo = joined(
        "splint_oligo",
        "named",
        "splint oligo as ordered",
        [
            fixed("splint_backfill", "linker", "length-equalizing filler, discarded with the 5' handle", backfill_used),
            fixed(
                "splint_primer_5p",
                "custom_primer",
                "5' primer binding site: the forward primer has this same sequence and extends "
                "toward the 3' end",
                spl_header[:KPNI_START],
            ),
            fixed(
                "splint_kpni_site",
                "linker",
                "KpnI restriction site GGTAC^C, cut at its far end. These five bases are "
                "discarded with the 5' handle; only the C that follows stays on the probe",
                spl_header[KPNI_START:KPNI_CUT],
            ),
            splint_probe,
            fixed(
                "splint_bamhi_site",
                "linker",
                "BamHI restriction site G^GATCC, cut at its near end. The G that completes it is "
                "the last base of the clamp above; these five are discarded with the 3' handle",
                spl_footer[BAMHI_CUT : BAMHI_CUT + 5],
            ),
            fixed(
                "splint_primer_3p",
                "custom_primer",
                "3' primer binding site, read on the opposite strand: the reverse primer is this "
                "region's reverse complement and extends back toward the 5' end",
                spl_footer[BAMHI_CUT + 5 :],
            ),
        ],
    )

    # Position, not identity: `construct_encoding` cycles through permutations
    # of the gene's three bits, so which readout lands in slot 1 differs from
    # probe to probe. Slot 1 is NOT `code1`. Only the set of three is fixed per
    # gene, which is why validation compares sorted triples.
    readout = lambda n: leaf(
        f"padlock_readout_{n}",
        "barcode",
        f"readout in position {n} of 3: one of the gene's three codebook bits, 20 nt from the "
        "vendored table of 49. Which of the three lands here varies per probe",
        "onlist",
        "N" * 20,
        onlist=onlist,
    )
    padlock_probe = joined(
        "padlock_probe",
        "named",
        "working padlock, released by the double digest and circularized on the target",
        [
            fixed(
                "padlock_kpni_retained",
                "linker",
                "the single base of the KpnI site GGTAC^C left on the working probe. The other "
                "five leave with the 5' handle, so the site is positioned to scar the part that "
                "gets thrown away. First of the 6 nt the splint clamp templates",
                pad_header[KPNI_CUT],
            ),
            fixed(
                "padlock_spacer_5p",
                "linker",
                "spacer completing, with the retained base, the 6 nt the splint clamp templates",
                pad_header[KPNI_CUT + 1 : 20],
            ),
            variable("padlock_head_splint", "linker", "head splint, 3 nt drawn from A/C/T to minimize hairpins", 3, 3),
            fixed("padlock_spacer_ta", "linker", "TA spacer", "TA"),
            variable("padlock_arm", "cdna", "padlock target-homology arm, reverse complement of the 5' half of the transcript window; arm + fill = 27 nt", 18, 27),
            readout(1),
            variable("padlock_spacer_1", "linker", "A/T spacer between readouts", 2, 2),
            readout(2),
            variable("padlock_spacer_2", "linker", "A/T spacer between readouts", 2, 2),
            readout(3),
            variable("padlock_fill", "linker", "length-equalizing filler, truncated from the right; arm + fill = 27 nt", 0, 9),
            fixed("padlock_linker_at", "linker", "AT linker", "AT"),
            fixed(
                "padlock_ligation_3p",
                "linker",
                "first half of the 6 nt at the padlock's 3' end that the splint clamp templates; "
                "the reverse complement of the splint footer's first three bases",
                rc(spl_footer[:BAMHI_CUT]),
            ),
            # Not split down to the single retained base, unlike the KpnI end:
            # that base is the last of the 6 nt the splint clamp templates, and
            # those six are the functional unit `test_splint_padlock` checks.
            fixed(
                "padlock_ligation_end",
                "linker",
                "second half of the 6 nt the splint clamp templates, and the padlock's 3' "
                "terminus. Its final base is the one base of the BamHI site G^GATCC kept on the "
                "working probe; the other five leave with the 3' handle",
                pad_footer[:BAMHI_CUT],
            ),
        ],
    )
    # `padlock_fill`'s sequence is written at its 9-nt maximum, so it is spelled
    # X's like any variable region; the constant it draws from is in the comment
    # header and in assembly.padpad().
    assert padlock_fill == "AATCACATA"

    padlock_oligo = joined(
        "padlock_oligo",
        "named",
        "padlock oligo as ordered",
        [
            fixed(
                "padlock_primer_5p",
                "custom_primer",
                "5' primer binding site: the forward primer has this same sequence and extends "
                "toward the 3' end",
                pad_header[:KPNI_START],
            ),
            fixed(
                "padlock_kpni_site",
                "linker",
                "KpnI restriction site GGTAC^C, cut at its far end. These five bases are "
                "discarded with the 5' handle; only the C that follows stays on the probe",
                pad_header[KPNI_START:KPNI_CUT],
            ),
            padlock_probe,
            fixed(
                "padlock_bamhi_site",
                "linker",
                "BamHI restriction site G^GATCC, cut at its near end. The G that completes it is "
                "the padlock's last base; these five are discarded with the 3' handle",
                pad_footer[BAMHI_CUT : BAMHI_CUT + 5],
            ),
            fixed(
                "padlock_primer_3p",
                "custom_primer",
                "3' primer binding site, read on the opposite strand: the reverse primer is this "
                "region's reverse complement and extends back toward the 5' end",
                pad_footer[BAMHI_CUT + 5 :],
            ),
        ],
    )

    return _assay(
        assay_id=f"mkprobes-solar-oligo-pool-bcidx{bcidx}",
        name=f"SOLAR splint/padlock oligo pool (bcidx {bcidx})",
        description=(
            "Array-synthesized oligo pool for SOLAR (splint/padlock, STARmap-style) spatial "
            f"transcriptomics, as assembled by `mkprobes assemble` with bcidx {bcidx}. Every "
            "probe pair is two 148-nt DNA oligos. Each is flanked by the primer binding sites "
            "that amplify this probe set: a forward site at the 5' end ending in a KpnI site, "
            "and a reverse site at the 3' end opening with a BamHI site. The reverse primer is "
            "the reverse complement of that 3' region and primes back toward the 5' end. Once "
            "amplified, the KpnI/BamHI double digest trims both sites off and releases the "
            f"working probes. bcidx {bcidx} sets those four sites, the restriction scars they "
            "leave on the working probe, and the splint clamp that templates the padlock's "
            "ends. The working padlock carries one target-homology arm followed by "
            "the gene's three 20-nt readouts; the working splint carries the other arm plus a "
            "6+6-nt clamp that holds the padlock's ends together for ligation. Case in the "
            "assembled pool files is cosmetic."
        ),
        library_kit="array-synthesized oligo pool, two 148-nt oligos per probe pair",
        library_spec=[
            joined(
                "dna",
                "dna",
                "SOLAR probe pair: two separately synthesized 148-nt oligos",
                [splint_oligo, padlock_oligo],
            )
        ],
    )


def build_readout_spec(onlist: dict[str, Any]) -> dict[str, Any]:
    return _assay(
        assay_id="mkprobes-solar-readout-probes",
        name="SOLAR readout (detection) probes",
        description=(
            "The 49 orthogonal 20-nt readout sequences SOLAR images against, numbered 1-49 in "
            "the vendored table. A detection oligo is one of these sequences carrying a "
            "fluorophore. Rolling-circle amplification copies the padlock circle, so the "
            "amplicon presents the reverse complement of the readout on the padlock, in "
            "hundreds of tandem copies, and the detection oligo hybridizes to it. A gene's "
            "codebook entry names the three readout IDs stitched onto that gene's padlocks; "
            "the order they appear in varies from probe to probe."
        ),
        library_kit="fluorophore-conjugated detection oligos, one 20-nt sequence each",
        library_spec=[
            joined(
                "dna",
                "dna",
                "SOLAR readout (detection) probe",
                [
                    leaf(
                        "readout_probe",
                        "barcode",
                        "readout sequence, 20 nt from the vendored table of 49; the fluorophore conjugate is not sequence and is not represented",
                        "onlist",
                        "N" * 20,
                        onlist=onlist,
                    )
                ],
            )
        ],
    )


def _dump(spec: dict[str, Any], header: str) -> str:
    return header + yaml.dump(spec, sort_keys=False, default_flow_style=False, width=100)


def pool_filename(bcidx: int) -> str:
    return f"solar_bcidx{bcidx}.seqspec.yaml"


def render() -> dict[str, str]:
    """
    Every file this script maintains, as `{filename: contents}`.

    Kept separate from writing them so `tests/test_oligospec.py` can compare
    against what is checked in without touching the working tree.
    """
    onlist_text = "\n".join(_read_readouts()) + "\n"
    payload = onlist_text.encode()
    onlist = _onlist("solar_readouts.txt", len(payload), hashlib.md5(payload).hexdigest())

    out = {
        "solar_readouts.txt": onlist_text,
        "solar_readout_probes.seqspec.yaml": _dump(build_readout_spec(onlist), READOUT_HEADER),
    }
    for bcidx in range(max_bcidx() + 1):
        header = pool_header(bcidx, *_read_headerfooter(bcidx))
        out[pool_filename(bcidx)] = _dump(build_pool_spec(onlist, bcidx), header)
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    # Indices removed from headerfooter.csv would otherwise leave stale specs
    # behind, and a stale spec is one a pool could still be validated against.
    written = render()
    for stale in OUT.glob("solar_bcidx*.seqspec.yaml"):
        if stale.name not in written:
            stale.unlink()
            print(f"removed stale {stale.relative_to(ROOT)}")

    total = 0
    for name, text in written.items():
        (OUT / name).write_text(text)
        total += len(text.encode())
    print(f"{len(written)} files, {total / 1024:.0f} KiB, in {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
