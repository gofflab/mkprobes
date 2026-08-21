"""
SOLAR probeset assembly: manifest of ProbeSets -> orderable oligo pool.

Package integration of scripts/probegen/2_assemble_manifest.py (verbatim
logic): `short` triages under-provisioned genes (with interactive off-target
acceptance into <codebook>.acceptable.json), `gen` assembles final
splint/padlock oligos (header/footer stitching, optional RepeatMasker,
BamHI/KpnI-free enforcement, geometry and length assertions) and writes the
orderable pool plus a provenance sidecar.

CLI: ``mkprobes assemble <manifest.json> {short|gen}``.
"""

import json
import subprocess
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from importlib.resources import files
from itertools import cycle, islice
from pathlib import Path

import numpy as np
import polars as pl
import pyfastx
import questionary
import rich
import rich.rule
import rich_click as click
from Bio import Seq
from Bio.Restriction import BamHI, KpnI  # type: ignore
from loguru import logger

from .codebook.codebook import ProbeSet
from .constants import RESTRICTION_TOKEN
from .starmap.starmap import generate_head_splint, test_splint_padlock
from .utils._alignment import gen_fasta
from .utils.provenance import encode, provenance_record, read_provenance
from .utils.sequtils import reverse_complement as rc

pl.Config.set_fmt_str_lengths(100)
# Splint/padlock header/footer table vendored inside the package;
# override with --headerfooter on the CLI group.
DEFAULT_HEADERFOOTER = Path(str(files("mkprobes") / "data" / "headerfooter.csv"))
hfs = pl.read_csv(DEFAULT_HEADERFOOTER)
console = rich.get_console()
species_mapping = {"mouse": "mus musculus", "human": "homo sapiens"}
# %%


# def until_first_g(seq: str, target: str = "G"):
#     r, target = rc(seq.upper()), target.upper()
#     res, truncated = r[:6], r[6:]
#     res += "" if res[-1] == target else truncated[: truncated.index(target) + 1]
#     if len(res) > 12:
#         raise ValueError("No G found")
#     assert res[-1] == target
#     return res


def backfill(seq: str, target: int = 148):
    return (
        "TTCCACTAACTCACATGTCATGCATTATCTTCTATACCTCTGAGCAGATCAGTAGTCTATTACATGCTCGTAGTACCGTAAGCCAGATAC"[
            : max(0, target - len(seq))
        ]
        + seq
    )


def run(
    path: Path,
    probeset: ProbeSet,
    n: int = 16,
    toolow: int = 4,
    low: int = 12,
    rm_species: str | None = None,
    skip_repeatmasker: bool = False,
):
    rand = np.random.default_rng(0)
    idx = probeset.bcidx
    codebook = probeset.load_codebook(path)
    logger.info(f"Loaded {probeset.codebook} with {len(codebook)} genes.")

    tss = list(codebook)
    dfs_ = []
    cols = []
    bads = []
    lows = []

    design_source: Path | None = None

    for ts in tss:
        try:
            final = path / f"output/{ts}_final_{RESTRICTION_TOKEN}_{','.join(map(str, sorted(codebook[ts])))}.parquet"
            df = pl.read_parquet(final).sort(
                [
                    pl.col("priority").list.max(),
                    pl.col("hp").list.min(),
                ],
                # Stable sort: ties previously landed in arbitrary order, making
                # the selected top-n (and thus the oligo pool) nondeterministic.
                maintain_order=True,
            )[:n]
        except FileNotFoundError as e:
            logger.critical(e)
            continue

        if len(df) < toolow:
            bads.append({"name": ts, "count": len(df)})
            logger.error(f"Too few probes ({len(df)}) for {ts}.")
            continue

        if len(df) < toolow * 2:
            # Resample to prevent probe dropout. Capping at 3x coverage.
            lows.append({"name": ts, "count": len(df)})
            df = pl.concat([df, df[: min(low - len(df), len(df) * 2)]])

        if not cols:
            cols = df.columns
            design_source = final

        if df["gene"].dtype == pl.List:
            df = df.with_columns(gene=pl.col("gene").list.get(0))

        dfs_.append(df[:, cols])

    # if len(dfs_) != len(tss):
    #     raise Exception("Gene count mismatch")

    dfs = pl.concat(dfs_)

    outpath = Path(path / "generated" / probeset.name)
    outpath.mkdir(exist_ok=True, parents=True)
    # RepeatMasker: explicit taxon (--rm-species) > built-in mouse/human
    # mapping > skip. --skip-repeatmasker silences the skip warning.
    rm_taxon = rm_species or species_mapping.get(probeset.species or "")
    if skip_repeatmasker:
        rm_taxon = None
        logger.info("RepeatMasker skipped (--skip-repeatmasker).")
    if rm_taxon:
        with ThreadPoolExecutor() as exc:
            for col_name in ["splint", "padlock"]:
                (outpath / f"{col_name}.fasta").write_text(gen_fasta(dfs[col_name]).getvalue())
                exc.submit(
                    subprocess.run,
                    f'RepeatMasker -pa 16 -norna -s -no_is -species "{rm_taxon}" {outpath / f"{col_name}.fasta"}',
                    shell=True,
                    check=True,
                )

        dfs = dfs.with_columns({
            col_name: [seq for name, seq in pyfastx.Fastx((outpath / f"{col_name}.fasta.masked").as_posix())]
            for col_name in ["splint", "padlock"]
            if (outpath / f"{col_name}.fasta.masked").exists()
        }).filter(~pl.col("splint").str.contains("N") & ~pl.col("padlock").str.contains("N"))
    elif not skip_repeatmasker:
        logger.warning(
            f"No RepeatMasker taxon for species {probeset.species!r}; skipping. "
            "Pass --rm-species '<taxon>' to run it (any taxon RepeatMasker's library supports), "
            "or --skip-repeatmasker to silence this warning."
        )

    if len(bads):
        msg = f"Too few probe pairs ({toolow=}) for {len(bads)} genes.\n{pl.DataFrame(bads).sort('count')}"
        logger.error(msg)

    if len(lows):
        msg = f"Low count genes.\n{pl.DataFrame(lows).sort('count')}"
        logger.warning(msg)

    # Before
    spl_idx = idx * 2
    pad_idx = idx * 2 + 1

    def padpad(s: str, target: int = 99):
        if len(s) > target + 2:
            raise ValueError("Too long")
        if len(s) > target:
            return s
        return s + "AATCACATAAAT"[: target - len(s)]

    # This is padlock.
    logger.info("Generating splint header.")
    # `rand` and the `it` cycle below are shared across rows, so the order in
    # which they are consumed decides the sequence each row gets. polars runs
    # `map_elements` across threads, so drawing from them inside the UDF is not
    # reproducible - repeated runs of the same panel emitted different oligos.
    # Drawing here, in row order, reproduces exactly what the original single-
    # threaded run produced, at any thread count.
    res = dfs.with_columns(
        _head_splint=pl.Series(
            [generate_head_splint(padlock, rand) for padlock in dfs["padlock"]], dtype=pl.Utf8
        )
    ).with_columns(
        pad_cut=(
            # head
            hfs[pad_idx, "header"][-3:]
            + pl.col("_head_splint").str.to_lowercase()
            + "ta"  # what the paper uses
            + pl.col("seq").map_elements(rc, return_dtype=pl.Utf8)
        ).map_elements(padpad, return_dtype=pl.Utf8)
        + "at"
        + rc(hfs[spl_idx, "footer"][:3])
        + hfs[pad_idx, "footer"][:3]
    )

    it = cycle("ATAAT")

    def splint_pad(seq: str, target: int = 47):
        if len(seq) > target:
            return seq
        return "".join(islice(it, target - len(seq))) + seq

    # Splint
    res = res.with_columns(
        _spl_unpadded=(
            # "TGTTGATGAGGTGTTGATGAATA"
            pl.col("splint").map_elements(rc, return_dtype=pl.Utf8)
            + "ca"
            + pl.col("pad_cut").str.slice(0, 6).map_elements(rc, return_dtype=pl.Utf8)
            + pl.col("pad_cut").str.slice(-6, 6).map_elements(rc, return_dtype=pl.Utf8)
        )
    )
    res = (
        res.with_columns(
            # Padded in row order, for the reason given above the head splint.
            spl_cut=pl.Series([splint_pad(s) for s in res["_spl_unpadded"]], dtype=pl.Utf8)
        ).drop("_head_splint", "_spl_unpadded")
    ).filter(
        (
            pl.col("spl_cut")
            .map_elements(lambda x: BamHI.search(Seq.Seq(x)), return_dtype=pl.List(pl.UInt32))
            .list.len()
            == 0
        )
        & (
            pl.col("spl_cut")
            .map_elements(lambda x: KpnI.search(Seq.Seq(x)), return_dtype=pl.List(pl.UInt32))
            .list.len()
            == 0
        )
        & (
            pl.col("pad_cut")
            .map_elements(lambda x: BamHI.search(Seq.Seq(x)), return_dtype=pl.List(pl.UInt32))
            .list.len()
            == 0
        )
        & (
            pl.col("pad_cut")
            .map_elements(lambda x: KpnI.search(Seq.Seq(x)), return_dtype=pl.List(pl.UInt32))
            .list.len()
            == 0
        )
    )

    def double_digest(s: str) -> str:
        return BamHI.catalyze(KpnI.catalyze(Seq.Seq(s))[1])[0].__str__()

    # These guard the molecule that gets synthesized, so they are raised rather
    # than asserted: `python -O` would strip an assert and ship the bad oligo.
    for s, r in zip(res["spl_cut"], res["pad_cut"]):
        if not test_splint_padlock(s, r, lengths=(6, 6)):
            raise ValueError(
                "Splint does not template the padlock's ends (6 nt each side), so the "
                f"padlock could not be circularised.\n  splint:  {s}\n  padlock: {r}"
            )

    out: pl.DataFrame = res.with_columns(
        # restriction scar already accounted for
        splintcons=hfs[spl_idx, "header"] + pl.col("spl_cut") + hfs[spl_idx, "footer"][3:],
        padlockcons=hfs[pad_idx, "header"][:-3].lower() + pl.col("pad_cut") + hfs[pad_idx, "footer"][3:],
    ).with_columns(splintcons=pl.col("splintcons").map_elements(backfill, return_dtype=pl.Utf8))

    for s, r in zip(out["splintcons"], out["padlockcons"]):
        if not test_splint_padlock(*map(double_digest, (s, r)), lengths=(6, 6)):
            raise ValueError(
                "After the KpnI/BamHI double digest the splint no longer templates the "
                f"padlock's ends.\n  splint:  {s}\n  padlock: {r}"
            )

    lengths = out["padlockcons"].str.len_chars()
    if not lengths.is_between(139, 150).all():
        raise ValueError(
            f"Padlock oligos must be 139-150 nt to synthesize; got "
            f"{lengths.min()}-{lengths.max()} nt. The header/footer table and the "
            "probe length have to agree - check --headerfooter."
        )

    from .codebook.codebook import hash_codebook

    # Carry the design parameters forward from a per-gene construct output, so the
    # ordered pool records the thresholds it was built under and not just its own.
    design = read_provenance(design_source) if design_source else None
    record = provenance_record(
        stage="assemble",
        probeset=probeset.model_dump(),
        codebook_hash=hash_codebook(codebook),
        n_probe_pairs=len(out),
        repeatmasker=rm_taxon or "skipped",
        headerfooter=str(DEFAULT_HEADERFOOTER),
        design=design,
    )

    (gen_path := path / "generated").mkdir(exist_ok=True, parents=True)
    out.write_parquet(gen_path / (probeset.name + ".parquet"), metadata=encode(record))
    logger.info(f"{len(out)} probe pairs written to {gen_path / (probeset.name + '.parquet')}")

    (gen_path / (probeset.name + "_pad.fasta")).write_text(
        gen_fasta(out["padlockcons"], names=range(len(out))).getvalue()
    )
    (gen_path / (probeset.name + "_splint.fasta")).write_text(
        gen_fasta(out["splintcons"], names=range(len(out))).getvalue()
    )

    (gen_path / (probeset.name + ".provenance.json")).write_text(
        json.dumps(record, indent=2, default=str) + "\n"
    )
    return out


# %%
@click.group("assemble")
@click.argument("manifest", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--headerfooter",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=f"Override the splint/padlock header/footer table (default: vendored {DEFAULT_HEADERFOOTER.name}).",
)
@click.pass_context
def cli(ctx: click.Context, manifest: Path, headerfooter: Path | None):
    """Turn designed probes into an orderable oligo pool.

    MANIFEST is a JSON list of probe sets, each naming a codebook and the
    header/footer row to build against. Use `gen` to assemble the pool, or
    `short` first to triage targets that came up thin.
    """
    from .init_project import check_manifest

    ctx.ensure_object(dict)
    if headerfooter is not None:
        global hfs
        hfs = pl.read_csv(headerfooter)
        logger.info(f"Using header/footer table {headerfooter}.")
    # Validated here rather than deep in `gen`: an out-of-range bcidx or a
    # missing codebook used to surface only after the pool had been built.
    try:
        mfs = check_manifest(manifest)
    except ValueError as e:
        raise click.ClickException(str(e)) from e
    ctx.obj["manifest"] = mfs
    ctx.obj["path"] = manifest.parent


def handle_checks(ts: str, offtargets: pl.DataFrame):
    offtargets = offtargets.with_columns(
        label=pl.when(pl.col("transcript_name").is_not_null())
        .then(pl.col("transcript_name"))
        .otherwise(pl.col("transcript"))
        + " "
        + pl.col("count").cast(pl.Utf8)
        + pl.when(pl.col("acceptable")).then(pl.lit(" (already ok)")).otherwise(pl.lit(""))
    )
    if offtargets.is_empty():
        return []
    selected_files = questionary.checkbox(
        f"Select acceptable genes for {ts}", choices=offtargets["label"]
    ).ask()
    return sorted({x.split(" ")[0].rsplit("-", 1)[0] for x in selected_files})


def manual_accept(path: Path, probeset: ProbeSet, *, ts: str):
    try:
        offtargets = pl.read_csv(path / f"output/{ts}_offtarget_counts.csv")
    except FileNotFoundError:
        ...
    else:
        path_accept = path / (probeset.codebook.rsplit(".", 1)[0] + ".acceptable.json")
        curr = json.loads(path_accept.read_text()) if path_accept.exists() else {}
        if ts in curr:
            logger.info(f"{ts} already in acceptable genes. Skipping manual check.")
        elif questionary.confirm(f"{ts} manual check?", default=True).ask():
            user_ok = handle_checks(ts, offtargets)
            if user_ok:
                path_accept.write_text(json.dumps({**curr, **{ts: user_ok}}, indent=2, sort_keys=True))
            logger.info("Outputted acceptable genes to " + str(path_accept))
            return user_ok
        return None


@cli.command()
@click.argument("short", type=int)
@click.option("--verbose", "-v", is_flag=True)
# @click.option("--output", type=click.Path(dir_okay=False, file_okay=True, path_type=Path))
@click.option(
    "--delete",
    is_flag=True,
    help="Delete the probes from a .tss.txt file that are too short.",
)
@click.pass_context
def short(
    ctx: click.Context, short: int, verbose: bool = False, delete: bool = False, permanent: bool = False
):
    """
    Identifies and optionally removes transcripts with fewer probes than a specified threshold.

    This function iterates through probe sets defined in the manifest. For each probe set,
    it loads the corresponding codebook and associated parquet files containing probe data.
    It then counts the number of probes per gene.

    If the '--delete' flag is set, genes with probe counts below the 'short' threshold
    are removed from a copy of the .tss.txt file.
    If '--permanent' is also set, the original .tss.txt file is overwritten.
    Otherwise, a new file with the suffix '.tss.ok.txt' is created.

    Args:
        ctx: The Click context, containing the manifest and path.
        short: The minimum number of probes a gene must have.
        verbose: If True, prints detailed information about genes with too few probes.
        delete: If True, removes genes with too few probes from the .tss.txt file.

    Raises:
        ValueError: If genes marked for deletion are not found in the .tss.txt file.
    """

    manifest: list[ProbeSet] = ctx.obj["manifest"]
    path_main: Path = ctx.obj["path"]

    for probeset in manifest:
        print("\n")
        console.print(rich.rule.Rule(title=probeset.name, align="left"))
        baddies = []

        codebook = probeset.load_codebook(path_main)

        path = (path_main / probeset.codebook).parent

        tss = list(codebook)
        dfs_ = []

        for ts in tss:
            try:
                _df = pl.read_parquet(
                    path / f"output/{ts}_final_{RESTRICTION_TOKEN}_{','.join(map(str, sorted(codebook[ts])))}.parquet"
                )
                _df = _df.sort([pl.col("priority").list.min(), pl.col("hp").list.max()])
                dfs_.append(
                    _df.select([
                        "name",
                        "seq",
                        "code1",
                        "code2",
                        "code3",
                        "index",
                        "id",
                        "flag",
                        "transcript",
                        "pos",
                        "cigar",
                        "aln_score",
                        "aln_score_best",
                        "n_ambiguous",
                        "n_mismatches",
                        "n_opens",
                        "n_extensions",
                        "edit_distance",
                        "mismatched_reference",
                        "gene",
                        "transcript_ori",
                        "pos_start",
                        "pos_end",
                        "length",
                        "match",
                        "match_consec",
                        "pad_start",
                        "maps_to_pseudo",
                        "max_tm_offtarget",
                        "match_consec_all",
                        "ok_quad_c",
                        "ok_quad_a",
                        "ok_stack_c",
                        "ok_comp_a",
                        "gc_content",
                        "ok_gc",
                        "tm",
                        "hp",
                        "oks",
                        "priority",
                        "splint",
                        "padlock",
                        "seqori",
                    ])
                )
                # .sample(shuffle=True, seed=4, fraction=1)
                if len(_df) < short:
                    baddies.append(ts)
                    manual_accept(path, probeset, ts=ts)
            except FileNotFoundError:
                baddies.append(ts)
                logger.warning(
                    "File "
                    + f"output/{ts}_final_{RESTRICTION_TOKEN}_{','.join(map(str, sorted(codebook[ts])))}.parquet"
                    + " not found. This usually means that there are no probes for this gene."
                )
                manual_accept(path, probeset, ts=ts)

        # if not len(dfs_):
        #     logger.warning(f"No data for {probeset.name} at {probeset.codebook}")
        #     continue

        # dfs: pl.DataFrame = pl.concat(dfs_)
        # print(dfs)
        # counts = dfs.group_by(COL_NAME).len(name="count")

        if baddies:
            print(probeset.name)
            print("\n".join(sorted(baddies)))
            logger.warning(f"Found {len(baddies)} genes with fewer than {short} probes.")

            # path_accept = path / (probeset.name + ".acceptable.json")
            # curr = json.loads(path_accept.read_text()) if path_accept.exists() else {}
            # path_accept.write_text(json.dumps({**curr, **out_dict}, indent=2, sort_keys=True))
            # logger.info("Outputted acceptable genes to " + str(path_accept))

        else:
            logger.info(f"All genes have at least {short} probes.")

        if delete and baddies:
            tss = probeset.load_codebook(path_main)
            # Check
            genes = set(tss)
            baddies = set(baddies)
            if not baddies.issubset(genes):
                raise ValueError(f"{baddies - genes} not found in codebook file. Wrong file?")

            goodies = [ts for ts in tss if ts not in baddies]

            out_path = path.joinpath(probeset.name).with_suffix(".good.txt")
            out_path.write_text("\n".join(sorted(goodies)))
            logger.info(
                f"Outputted {len(goodies)} genes with >={short} probes to {out_path}. Please remake the codebook."
            )


@cli.command()
@click.option(
    "--rm-species",
    type=str,
    default=None,
    help="RepeatMasker taxon passed verbatim to -species (e.g. 'octopus', 'mollusca'). "
    "Overrides the built-in mouse/human mapping.",
)
@click.option(
    "--skip-repeatmasker",
    is_flag=True,
    help="Skip RepeatMasker explicitly (silences the no-taxon warning for non-model species).",
)
@click.pass_context
def gen(ctx: click.Context, rm_species: str | None, skip_repeatmasker: bool):
    mfs: list[ProbeSet] = ctx.obj["manifest"]
    path_main: Path = ctx.obj["path"]
    logger.info(f"{path_main=}")

    total_probes = 0
    for x in mfs:
        path = (path_main / x.codebook).parent
        if isinstance(x.n_probes, int):
            logger.info(f"Using {x.n_probes} probes for {x.name}.")
            n = x.n_probes
            low = min(n, 15)
        elif x.n_probes is not None:
            n = 34 if x.n_probes == "high" else 16
            low = 24 if x.n_probes == "high" else 12
        else:
            n = 34 if x.species == "human" else 24
            low = 24 if x.species == "human" else 16

        probes = run(path, x, n=n, toolow=4, low=low, rm_species=rm_species, skip_repeatmasker=skip_repeatmasker)

        total_probes += len(probes) * 2
        logger.info(f"Cumulative probes: {total_probes}")

    superout = []
    Path(path_main / "generated").mkdir(exist_ok=True, parents=True)
    for m in mfs:
        out = []
        path = (path_main / m.codebook).parent
        paths = [
            (path / "generated" / f"{m.name}_splint.fasta"),
            (path / "generated" / f"{m.name}_pad.fasta"),
        ]

        for i in range(2):
            if not paths[i].exists():
                raise FileNotFoundError(f"File {paths[i]} not found.")
                # paths[i] = paths[i].with_name(paths[i].name[:-7])

        for s, p in zip(*[pyfastx.Fastx(p.as_posix()) for p in paths]):
            if "N" not in s[1] and "N" not in p[1]:
                out.append(s[1])
                out.append(p[1])
        Path(path_main / "generated" / f"{m.name}_final.txt").write_text("\n".join(out))
        superout.extend(out)

    now_str = datetime.now().replace(microsecond=0).isoformat()
    Path(path_main / "generated" / f"_allout{now_str}.txt").write_text("\n".join(superout))


if __name__ == "__main__":
    cli()


# t7 = "TAATACGACTCACTATAGGG"
# assert out["padlockcons"].str.contains(rc(t7)[:5]).all()

# %%


# for name in ["genestarpad.fasta", "genestarsplint.fasta"]:

# cons = dfs.with_columns(constructed=header + pl.col("seq") + footer)
# cons = cons.with_columns(constructed=pl.col("constructed").map_elements(backfill))
# # %%
# import pyfastx


# # %%

# Path("starwork/genestar_out.txt").write_text("\n".join(out))


# %%

# %%
