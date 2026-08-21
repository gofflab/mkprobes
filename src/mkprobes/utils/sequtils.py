import random
import re
from typing import Any, Literal

import colorama
import numpy as np
import numpy.typing as npt
import polars as pl
from Bio.Seq import Seq
from matplotlib.axes import Axes

name_splitter = re.compile(r"(.+)_(.+):(\d+)-(\d+)")


def reject_ambiguous_bases(frame: pl.DataFrame, stage: str, column: str = "seq") -> None:
    """
    Raises if any sequence in `column` contains an ambiguous base.

    An N that reaches synthesis is a defect, not a warning, so this raises
    rather than asserts: `python -O` strips asserts and would ship the oligo.
    """
    offending = frame.filter(pl.col(column).str.contains("N"))
    if not len(offending):
        return
    raise ValueError(
        f"{len(offending)} sequence(s) contain an ambiguous base (N) after {stage}. "
        "The reference usually has masked or unresolved regions in this transcript. "
        f"First offender: {offending[column][0]}"
    )

# Probe names are structured `{gene}_{transcript}:{start}-{end}` with an
# optional `_{splint|padlock}` suffix appended at the melt step. De novo IDs
# may themselves contain underscores (TRINITY_DN123_c0_g1_i1) and dots
# (STRG.1.1), so parsing must not naively split on `_`.
_PROBE_NAME_RE = r"^(.+):(\d+)-(\d+)(?:_[A-Za-z]+)?$"


def probe_coord_exprs() -> list[pl.Expr]:
    """
    Polars expressions extracting gene / transcript_ori / pos_start / pos_end
    from a `name` column, robust to underscores inside IDs.

    The generic pipeline always builds the prefix as `{transcript}_{transcript}`
    (gene == transcript), so an even-duplicate prefix `X_X` unambiguously
    yields gene = transcript = X for ANY X — including underscore-rich de novo
    IDs. Otherwise (reference path, where IDs contain no underscores) the
    prefix splits at its last underscore.
    """
    prefix = pl.col("name").str.extract(_PROBE_NAME_RE, 1)
    plen = prefix.str.len_chars()
    half = ((plen - 1) // 2).cast(pl.Int64)
    first = prefix.str.slice(0, half)
    mid = prefix.str.slice(half, 1)
    second = prefix.str.slice(half + 1)
    is_dup = (plen >= 3) & (mid == pl.lit("_")) & (first == second)
    return [
        pl.when(is_dup).then(first).otherwise(prefix.str.extract(r"^(.+)_([^_]+)$", 1)).alias("gene"),
        pl.when(is_dup)
        .then(first)
        .otherwise(prefix.str.extract(r"^(.+)_([^_]+)$", 2))
        .alias("transcript_ori"),
        pl.col("name").str.extract(_PROBE_NAME_RE, 2).cast(pl.UInt32).alias("pos_start"),
        pl.col("name").str.extract(_PROBE_NAME_RE, 3).cast(pl.UInt32).alias("pos_end"),
    ]


def probe_identity_exprs() -> list[pl.Expr]:
    """
    Polars expressions splitting a suffixed probe name into
    full_name (`{gene}_{transcript}:{start}-{end}`) and probe_type
    (`splint`/`padlock`, the final `_`-segment) — robust to underscores
    inside IDs, unlike a positional split on `_`.
    """
    return [
        pl.col("name").str.extract(r"^(.+)_([^_]+)$", 1).alias("full_name"),
        pl.col("name").str.extract(r"^(.+)_([^_]+)$", 2).alias("probe_type"),
    ]


c = re.compile(r"(\d+)S(\d+)M")
c2 = re.compile(r"(\d+)M")


def parse_cigar(s: str, m_only: bool = False) -> tuple[int, ...]:
    if not m_only:
        try:
            return tuple(map(int, c.findall(s)[0]))
        except IndexError:
            try:
                return 0, max(map(int, (c2.findall(s))))
            except ValueError:
                return 0, 0
            except IndexError:
                return 0, 0
    return tuple(map(int, c2.findall(s)))


def printc(seq: str):
    for c in seq:
        if c == "A" or c == "a":
            print(colorama.Fore.GREEN + c, end="")
        elif c == "T" or c == "t":
            print(colorama.Fore.RED + c, end="")
        elif c == "C" or c == "c":
            print(colorama.Fore.BLUE + c, end="")
        elif c == "G" or c == "g":
            print(colorama.Fore.YELLOW + c, end="")
        else:
            print(colorama.Fore.WHITE + c, end="")
    print(colorama.Fore.RESET)


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    # https://bioinformatics.stackexchange.com/a/3585
    return Seq(seq).reverse_complement().__str__()


def pcr(seq: str, primer: str, primer_rc: str) -> str:
    loc = seq.find(primer)
    if loc == -1:
        raise ValueError(f"Primer {primer} not found in sequence {seq}")
    loc_rc = reverse_complement(seq[loc:]).find(primer_rc)
    if loc_rc == -1:
        raise ValueError(f"Primer {primer_rc} not found in sequence {seq}")
    return seq[loc : None if loc_rc == 0 else -loc_rc]


def is_subsequence(sub_dna: str):
    iupac_dict = {
        "R": "[AG]",
        "Y": "[CT]",
        "S": "[GC]",
        "W": "[AT]",
        "K": "[GT]",
        "M": "[AC]",
        "B": "[CGT]",
        "D": "[AGT]",
        "H": "[ACT]",
        "V": "[ACG]",
        "N": "[ACGT]",
    }

    # Convert sub_dna to regex
    sub_dna_regex = ""
    for base in sub_dna:
        if base in iupac_dict:
            sub_dna_regex += iupac_dict[base]
        else:
            sub_dna_regex += base

    # Check if sub_dna is in main_dna
    sub_dna_regex = re.compile(sub_dna_regex)

    def inner(main_dna: str):
        if match := sub_dna_regex.search(main_dna):
            return match.span()
        return None

    return inner


def gen_random_base(n: int) -> str:
    """Generate a random DNA sequence of length n."""
    return "".join(random.choices("ACGT", k=n))


def gc_content(seq: str) -> float:
    return (seq.count("G") + seq.count("C")) / len(seq)


def equal_distance(total: int, choose: int) -> npt.NDArray[np.int_]:
    return np.linspace(0, total - 1, choose).astype(np.int_)


def plot_gc_content(ax: Axes, seq: str, window_size: int = 50, **kwargs: Any):
    """Plot windowed GC content on a designated Matplotlib ax."""
    if len(seq) < window_size:
        raise ValueError("Sequence shorter than window size")
    out = np.zeros(len(seq) - window_size + 1, dtype=float)
    curr = gc_content(seq[:window_size]) * window_size
    out[0] = curr
    for i in range(1, len(seq) - window_size):
        curr += (seq[i + window_size - 1] in "GC") - (seq[i - 1] in "GC")
        out[i] = curr
    out /= window_size

    ax.fill_between(np.arange(len(seq) - window_size + 1), out, **(dict(alpha=0.3) | kwargs))  # type: ignore
    ax.set_ylim(bottom=0, top=1)
    ax.set_ylabel("GC (%)")


def stripplot(**kwargs: Any) -> Axes:
    import pandas as pd
    import seaborn as sns

    sns.set()

    df = pd.concat(pd.DataFrame({"x": v, "y": k}) for k, v in kwargs.items())
    return sns.stripplot(data=df, x="x", y="y", **(dict(orient="h", alpha=0.6)))  # type: ignore


def gen_idt(name: str, seq: str, scale: str = "25nm", purification: str = "STD") -> str:
    return f"{name}\t{seq}\t\t{scale}\t{purification}"


def gen_plate(name: str | list[str], seqs: list[str], order: Literal["C", "F"] = "C") -> pl.DataFrame:
    if order == "C":
        wells = [f"{row}{col:02d}" for row in "ABCDEFGH" for col in range(1, 13)]
    else:
        wells = [f"{row}{col:02d}" for col in range(1, 13) for row in "ABCDEFGH"]

    return pl.DataFrame(
        {
            "Well Position": wells[: len(seqs)],
            "Name": name,
            "Sequence": seqs,
        }
    )


def kmers(seq: str, n: int) -> list[str]:
    return [seq[i : i + n] for i in range(len(seq) - n + 1)]
