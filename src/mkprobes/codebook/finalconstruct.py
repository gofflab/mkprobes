# %%
import json
import re
from itertools import chain, cycle, permutations
from pathlib import Path
from typing import Annotated, Callable, Collection, Final, Iterable, Sequence, cast

import click
import polars as pl
from Bio import Restriction, Seq
from loguru import logger

from ..candidates import _run_bowtie
from ..ext.dataset import Dataset, ReferenceDataset, load_dataset
from ..starmap.starmap import rotate, test_splint_padlock
from ..utils.provenance import provenance_metadata
from ..utils.sequtils import reject_ambiguous_bases
from ..utils.sequtils import reverse_complement as rc

READOUTS: Final[dict[int, str]] = {
    x["id"]: x["seq"] for x in pl.read_csv(Path(__file__).parent / "readout_ref_filtered.csv").to_dicts()
}


def _check_pad_start(name: str, pad_start: int) -> None:
    """
    The padlock arm has to begin past position 17, or there is no room for the
    readouts. Raised rather than asserted so `python -O` cannot strip it.
    """
    if pad_start <= 17:
        raise ValueError(
            f"Probe {name} has its padlock arm starting at {pad_start}, but readout "
            "attachment needs it past 17. Re-run screening for this target."
        )


def assign_overlap(
    output: Path | str,
    gene: str,
    *,
    target_probes: int = 16,
    max_overlap: int = 5,
    restriction: str = "BsaI",
) -> int:
    if max_overlap % 5 != 0 or max_overlap < 0:
        raise ValueError("max_overlap must be a multiple of 5")
    output = Path(output)
    for ol in chain((-2,), range(5, max_overlap + 1, 5)):
        df = pl.read_parquet(output / f"{gene}_screened_ol{ol}{restriction}.parquet")
        if len(df) >= target_probes:
            return ol

    # if len(df) >= min_probes:  # type: ignore
    return ol  # type: ignore

    # raise ValueError(f"Gene {gene} cannot be fixed")


def stitch(seq: str, codes: Sequence[int], sep: str = "TT") -> str:
    return sep.join(rc(READOUTS[c]).lower() for c in codes) + seq


def construct_idt(seq_encoding: pl.DataFrame, idxs: Sequence[int]):
    def gen_starpad(bit: int, s: tuple[str, str]):
        for x in "atcg":
            out_pad = rotate(
                READOUTS[bit].lower() + x + s[1] + "tattcaat"[: max(0, 46 - len(s[1]) - 20)].lower(), 12 + 1
            )
            out_pad_upper = out_pad.upper()
            if all(seq not in out_pad_upper for seq in ["AAAAA", "TTTTT", "CCCCC", "GGGGG"]):
                break
        else:
            raise ValueError("Homopolymers")

        splint = s[0] + (footer := "ta" + rc(out_pad[:6]) + rc(out_pad[-6:]))
        # Raised, not asserted: `python -O` would strip these and emit a padlock
        # that cannot circularise.
        if not test_splint_padlock(footer[2:], out_pad):
            raise ValueError(f"Splint does not template the padlock's ends: {out_pad}")
        if not 46 <= len(out_pad) <= 48:
            raise ValueError(f"Padlock payload must be 46-48 nt, got {len(out_pad)}: {out_pad}")
        return splint, out_pad

    assert len(idxs) == 1
    out = dict(name=[], code=[], cons_pad=[], cons_splint=[], seq=[])

    for name, splint, pad, pad_start in seq_encoding[["name", "splint", "padlock", "pad_start"]].iter_rows():
        _check_pad_start(name, pad_start)

        out["name"].append(name)  # f"{name};;{sep}{','.join(map(str,codes))}")
        out["code"].append(idxs[0])
        cons_splint, cons_pad = gen_starpad(idxs[0], (rc(splint), rc(pad)))
        out["cons_pad"].append(cons_pad)
        out["cons_splint"].append(cons_splint)
        out["seq"].append(cons_pad)

    return pl.DataFrame(out)


def construct_encoding(
    seq_encoding: Annotated[pl.DataFrame, ["name", "padlock", "pad_start"]], idxs: Sequence[int], n: int = -1
):
    if len(idxs) == 1:
        return construct_idt(seq_encoding, idxs)
    if n == -1:
        n = len(idxs)
    if len(idxs) != n:
        raise ValueError(f"Invalid number of readouts: {idxs}")
    if any(idx not in READOUTS for idx in idxs):
        raise ValueError(f"Invalid readout indices: {idxs}")

    perms = cast(Iterable[tuple[int, ...]], cycle(permutations(idxs, n)))

    out = dict(name=[], seq=[])
    for i in range(n):
        out[f"code{i + 1}"] = []

    for name, pad, pad_start in seq_encoding[["name", "padlock", "pad_start"]].iter_rows():
        # for codes, _ in zip(perms, range(4)):
        _check_pad_start(name, pad_start)
        for sep, codes in zip(["AA", "TA", "AT", "TT"], perms):
            stitched = stitch(pad, codes, sep=sep)
            stitched_upper = stitched.upper()
            if any(hp in stitched_upper for hp in ["AAAAA", "TTTTT", "CCCCC", "GGGGG"]):
                continue
            if Restriction.BamHI.search(Seq.Seq(stitched)):
                continue
            out["name"].append(name)  # f"{name};;{sep}{','.join(map(str,codes))}")
            out["seq"].append(stitched)
            for i, code in enumerate(codes):
                out[f"code{i + 1}"].append(code)
            break  # only one separator is needed. This is a split design.

    return pl.DataFrame(out)


def check_offtargets(dataset: Dataset, constructed: pl.DataFrame, acceptable_tss: list[str]):
    sam = _run_bowtie(dataset, constructed, ignore_revcomp=True)[0]
    acc = sam.agg_tm_offtarget(acceptable_tss)
    logger.debug(acc)
    df = (
        acc  # .filter(
        # pl.col("max_tm_offtarget").lt(40)
        #  ~pl.col("seq").map_elements(lambda x: dataset.check_kmers(cast(str, x)))
        #
        .drop("seq")
        .join(constructed, on="name", suffix="padlock", how="inner")
        .with_columns(name=pl.col("name").str.split(";;").list.first())
        .sort("match", descending=True)
        .unique("name", keep="first", maintain_order=True)
    )
    logger.info(f"Padlocks w/o offtarget: {len(df)}")
    return df


# %%


# fmt: off
@click.command("construct")
@click.argument("path", type=click.Path(exists=True, dir_okay=True, file_okay=False, path_type=Path))
@click.argument("output_path", type=click.Path(dir_okay=True, file_okay=False, path_type=Path))
@click.option("--gene", "-g", type=str, required=True, help="Target transcript to build probes for.")
@click.option("--codebook", "-c", required=True, type=click.Path(exists=True, dir_okay=False, file_okay=True, path_type=Path), help="Codebook JSON assigning readout bits to each target.")
@click.option("--target-probes", "--target_probes", "-N", type=int, help="Maximum number of probes per gene", default=72, show_default=True)
@click.option("--restriction", multiple=True, type=str, help="Restriction enzymes to exclude sites for. Repeatable.")
# fmt: on
def click_construct(
    path: Path,
    output_path: Path,
    gene: str,
    codebook: Path,
    target_probes: int = 72,
    restriction: list[str] | str | None = None,
):
    """Attach readout sequences to screened probes for one target.

    Reads that target's screened probes from OUTPUT_PATH and writes
    `<target>_final_<enzymes>_<bits>.parquet` beside them.
    """
    from ..constants import validate_restriction

    try:
        validate_restriction(restriction)
    except ValueError as e:
        raise click.BadParameter(str(e), param_hint="--restriction") from e

    construct(
        load_dataset(path),
        output_path,
        transcript=gene,
        codebook=json.loads(codebook.read_text()),
        target_probes=target_probes,
        restriction=restriction,
    )


def construct(
    dataset: Dataset,
    output_path: Path | str,
    *,
    transcript: str,
    codebook: dict[str, Collection[int]],
    target_probes: int = 64,
    restriction: list[str] | str | None = None,
    construction_function: Callable[[pl.DataFrame, Collection[int]], pl.DataFrame] = construct_encoding,
    overwrite: bool = False,
):
    output_path = Path(output_path)
    if isinstance(restriction, (list, tuple)) and restriction:
        restriction = "_" + "".join(restriction)
    restriction = restriction or ""

    if (final_path := Path(output_path / f"{transcript}_final{restriction}_{','.join(map(str, sorted(codebook[transcript])))}.parquet")).exists():
        if overwrite:
            final_path.unlink()
        else:
            return
    final_path.parent.mkdir(exist_ok=True, parents=True)
    #     pl.read_parquet(final_path)[["code1", "code2"]].to_numpy().flatten()
    # ) != set(codebook[transcript]):
    #     logger.critical(f"Codebook for {transcript} has changed.")
    # exit(1)

    # assign_overlap(output_path, transcript, target_probes=target_probes, restriction=restriction)
    overlap = -2
    screened = pl.read_parquet(
        scr_path := output_path / f"{transcript}_screened_ol{overlap}{restriction}.parquet"
    )
    logger.debug(f"Using {scr_path} for {transcript}.")
    logger.debug(f"Screened probes: {len(screened)}")

    reject_ambiguous_bases(screened, "screening")
    # acceptable_tss = pl.read_csv(next(output_path.glob(f"{transcript}_acceptable_tss.csv")))[
    #     "transcript_id"
    # ].to_list()
    # mrna = dataset.ensembl.get_seq(transcript)

    # screened = (
    #     screened.with_columns(splitted=pl.col("seq").map_elements(lambda pos: split_probe(pos, 58), return_dtype=pl.List(pl.Utf8)))
    #     .with_columns(
    #         splint=pl.col("splitted").list.get(1),  # need to be swapped because split_probe is not rced.
    #         padlock=pl.col("splitted").list.get(0),
    #         padstart=pl.col("splitted").list.get(2).cast(pl.Int16),
    #     ).drop("splitted")
    # )

    # screened = screened.filter(
    #     (pl.col("splint").str.len_chars() > 0)
    #     # & (pl.col("splint").map_elements(lambda x: hp(x, "dna")) < 50)
    #     # & (pl.col("padlock").map_elements(lambda x: hp(x, "dna")) < 50)
    # )

    logger.debug(f"With proper pad_start: {len(screened)}")

    res = construction_function(screened, codebook[transcript]).join(
        screened.rename(dict(seq="seqori")), on="name", how="left"
    )

    logger.info(f"Constructed {len(res)} probes for {transcript}.")
    if res["seq"].is_null().any():
        raise ValueError(
            f"{res['seq'].is_null().sum()} constructed probe(s) for {transcript} have no "
            "sequence. The screened input and the codebook disagree about this target."
        )

    final_path.parent.mkdir(exist_ok=True, parents=True)
    res.write_parquet(
        final_path,
        metadata=provenance_metadata(
            dataset.path,
            stage="construct",
            transcript=transcript,
            bits=sorted(codebook[transcript]),
            # `restriction` has already been folded into its filename form here.
            restriction=restriction.lstrip("_") or None,
            target_probes=target_probes,
        ),
    )
    # res.write_csv(final_path.with_suffix(".tsv"), separator="\t")  # deal with nested data
    logger.info(f"Written to {final_path}")
    return res


def count_final_probes(output_path: Path, gene: str) -> int | None:
    """
    Probes a target actually contributes to the pool, or `None` if it never
    reached the construct stage.

    Counts `_final_` rather than `_screened_` output: screening produces the
    candidate pool, construct decides how many of those survive readout
    attachment, and it is the latter that gets ordered.
    """
    finals = list(output_path.glob(f"{gene}_final_*.parquet"))
    if not finals:
        return None
    # One target can have several, if it was rebuilt under different enzymes or
    # bit assignments. The newest is the one the current panel refers to.
    newest = max(finals, key=lambda path: path.stat().st_mtime)
    return len(pl.read_parquet(newest))


@click.command("filter-genes")
@click.argument("output_path", type=click.Path(exists=True, dir_okay=True, file_okay=False, path_type=Path))
@click.option(
    "--genes",
    "-g",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, path_type=Path),
    help="Target list to check, one per line.",
)
@click.option(
    "--min-probes", type=int, help="Minimum number of probes per target", default=48, show_default=True
)
@click.option(
    "--out",
    "-o",
    type=click.Path(dir_okay=False, file_okay=True, path_type=Path),
    help="Write the targets that pass to this file, one per line.",
)
def filter_genes(output_path: Path, genes: Path, min_probes: int, out: Path | None = None):
    """Report how many probes each target ended up with, and flag the thin ones.

    Counts probes in each target's `_final_` output - what the target would
    contribute to an oligo order. Targets below --min-probes are warned about
    individually; rework or drop them before assembly.
    """
    gene_list = [line.strip() for line in genes.read_text().splitlines() if line.strip()]
    counts: dict[str, int] = {}
    missing: list[str] = []

    for gene in gene_list:
        count = count_final_probes(output_path, gene)
        if count is None:
            missing.append(gene)
            continue
        if count < min_probes:
            logger.warning(f"{gene} has {count} probes, less than {min_probes}.")
        counts[gene] = count

    if missing:
        logger.error(
            f"{len(missing)} target(s) have no constructed probes in {output_path}: "
            f"{', '.join(missing[:10])}{' ...' if len(missing) > 10 else ''}. "
            "Run `mkprobes run-panel` first, or check its `.failed.txt`."
        )

    passed = [gene for gene, n in counts.items() if n >= min_probes]
    logger.info(f"{len(passed)}/{len(gene_list)} target(s) have at least {min_probes} probes.")
    if out:
        out.write_text("\n".join(passed) + "\n" if passed else "")


# readouts = pl.read_csv("data/readout_ref_filtered.csv")
# # genes = Path("zach_26.txt").read_text().splitlines()
# # genes = Path(f"{name}_25.txt").read_text().splitlines()
# smgenes = ["Cdc42", "Neurog2", "Ccnd2"]
# genes, _, _ = gtf_all.check_gene_names(smgenes)
# acceptable_tss = {g: set(pl.read_csv(f"output/{g}_acceptable_tss.csv")["transcript"]) for g in genes}
# n, short_threshold = 67, 65
# # %%
# dfx, overlapped = {}, {}
# for gene in genes:
#     dfx[gene] = GeneFrame.read_parquet(f"output/{gene}_final.parquet")
# dfs = GeneFrame.concat(dfx.values())
# short = dfs.count("gene").filter(pl.col("count") < short_threshold)


# # %%
# fixed_n = {}
# short_fixed = {}


# def run_overlap(genes: Iterable[str], overlap: int):
#     def runpls(gene: str):
#         subprocess.run(
#             ["python", "scripts/new_postprocess.py", gene, "-O", str(overlap)],
#             check=True,
#             capture_output=True,
#         )

#     with ThreadPoolExecutor(32) as executor:
#         for x in as_completed(
#             [
#                 executor.submit(runpls, gene)
#                 for gene in genes
#                 if not Path(f"output/{gene}_final_overlap_{overlap}.parquet").exists()
#             ]
#         ):
#             print("ok")
#             x.result()


#     needs_fixing = set(short["gene"])

#     for ol in [5, 10, 15, 20]:
#         print(ol, needs_fixing)
#         run_overlap(needs_fixing, ol)
#         for gene in needs_fixing.copy():
#             df = GeneFrame.read_parquet(f"output/{gene}_final_overlap_{ol}.parquet")
#             if len(df) >= short_threshold or ol == 20:
#                 needs_fixing.remove(gene)
#                 fixed_n[gene] = ol
#                 short_fixed[gene] = df
#     # else:
#     #     raise ValueError(f"Gene {gene} cannot be fixed")

#     short_fixed = GeneFrame.concat(short_fixed.values())
#     # %%
#     cutted = GeneFrame.concat([dfs.filter(~pl.col("gene").is_in(short["gene"])), short_fixed[dfs.columns]])
#     cutted = GeneFrame(
#         cutted.sort(["gene", "priority"]).groupby("gene").agg(pl.all().head(n)).explode(pl.all().exclude("gene"))
#     )

# %%
