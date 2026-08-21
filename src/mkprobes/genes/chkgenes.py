import difflib
import json
import re
from collections.abc import Sequence
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Annotated, Literal

import polars as pl
import requests
import rich_click as click
from loguru import logger

from ..ext.dataset import Dataset, ReferenceDataset, load_dataset
from ..utils.targets import read_target_list
from ..ext.external_data import MockGTF, get_ensembl
from ..utils.printing import jprint

GENERIC_MODES = ("longest", "all")
REFERENCE_MODES = ("gencode", "ensembl", "canonical", "appris", "apprisalt")


def find_outdated_ts(ts: str) -> tuple[Annotated[str, "gene_name"], Annotated[str, "gene_id"]]:
    if not ts.startswith("ENST"):
        raise ValueError(f"{ts} is not a human Ensembl transcript ID")
    out = set()
    for x in requests.get(
        "https://dev-tark.ensembl.org/api/transcript/",
        params={
            "stable_id": ts,
            "source_name": "ensembl",
            "assembly_name": "GRCh38",
            "expand": "genes",
        },
    ).json()["results"]:
        for y in x["genes"]:
            if y["name"]:
                out.add((y["name"], y["stable_id"]))

    if len(out) != 1:
        logger.error(f"Found {len(out)} genes for {ts}: {out}")
    return list(out)[0]


@click.command()
@click.argument("path", type=click.Path(exists=True, dir_okay=True, file_okay=False, path_type=Path))
@click.argument("genes", type=click.Path(exists=True, dir_okay=False, file_okay=True, path_type=Path))
def chkgenes(path: Path, genes: Path):
    """Validate/check gene names (Ensembl for reference datasets; offline vs the GTF otherwise)"""
    from ..ext.fix_gene_name import check_gene_names

    ds = load_dataset(path)
    del path
    gs = read_target_list(genes)
    for gene in gs:
        if not gene.isascii():
            raise ValueError(f"{gene} not ASCII.")

    if len(gs) != len(s := set(gs)):
        [gs.remove(x) for x in s]
        logger.critical(f"Non-unique genes found: {', '.join(gs)}.\n")
        genes.with_suffix(".unique.txt").write_text("\n".join(sorted(list(s))))
        logger.error(f"Unique genes written to {genes.with_suffix('.unique.txt')}.\n")
        return

    if not isinstance(ds, ReferenceDataset):
        # Offline check against the dataset's own annotation (plus registered
        # alias/ortholog tables), no Ensembl/mygene network dependency.
        try:
            res = get_transcripts_generic(ds, gs, mode="all")
        except ValueError as e:
            raise SystemExit(str(e)) from None
        logger.info(f"{len(gs)} targets resolved to {res['transcript_id'].n_unique()} transcripts.")
        genes.with_suffix(".converted.txt").write_text("\n".join(sorted(set(gs))))
        logger.info(f"Written to {genes.with_suffix('.converted.txt')}.")
        return

    if not ds.ensembl:
        raise ValueError("Not a ReferenceDataset. Cannot check genes.")
    converted, mapping, no_fix_needed = check_gene_names(ds.ensembl, gs, species=ds.species)
    print(converted)
    if mapping:
        logger.info("Mappings:")
        jprint(mapping)
        logger.info(f"Mapping written to {genes.with_suffix('.mapping.json')}.")
        genes.with_suffix(".mapping.json").write_text(json.dumps(mapping))
        logger.info(f"Converted genes written to {genes.with_suffix('.converted.txt')}.")
        genes.with_suffix(".converted.txt").write_text("\n".join(sorted(converted)))
    elif not no_fix_needed:
        logger.warning("Some genes cannot be found.")
    else:
        logger.info(f"{len(s)} genes checked out. No changes needed")
        genes.with_suffix(".converted.txt").write_text("\n".join(sorted(converted)))


def _resolve_via_annotations(dataset: Dataset, token: str) -> list[str]:
    """
    Resolves a token (e.g. a human ortholog symbol) through the dataset's
    registered annotation tables.

    Searches every non-join column of each table for a case-insensitive exact
    match; hits map back through the table's `transcript_id`/`gene_id` join
    column. Returns matching values from the join columns (transcript IDs
    and/or gene IDs), empty if nothing matched.
    """
    hits: list[str] = []
    for name in sorted(dataset.annotation_paths):
        table = dataset.annotation(name)
        join_cols = [c for c in ("transcript_id", "gene_id") if c in table.columns]
        for col in table.columns:
            if col in join_cols or table[col].dtype != pl.Utf8:
                continue
            matched = table.filter(pl.col(col).str.to_lowercase() == token.lower())
            if not matched.is_empty():
                for jc in join_cols:
                    hits.extend(matched[jc].drop_nulls().to_list())
                logger.info(f"Resolved {token!r} via annotation table {name!r} column {col!r}.")
    return list(dict.fromkeys(hits))


def get_transcripts_generic(
    dataset: Dataset,
    tokens: Sequence[str],
    mode: Literal["longest", "all"] = "longest",
) -> pl.DataFrame:
    """
    Annotation-driven transcript selection for custom (non-reference) datasets.

    Each token is resolved in order: exact transcript ID (passthrough), exact
    gene name or gene ID (expands to that gene's isoforms), then a
    case-insensitive search of registered annotation tables (e.g. ortholog or
    alias mappings). `mode="longest"` keeps, per gene, the isoform with the
    longest sequence (from the FASTA, so introns don't distort the choice);
    `mode="all"` keeps every isoform.

    Returns:
        pl.DataFrame[[gene_name, gene_id, transcript_id, transcript_name]]
    """
    if not len(tokens):
        raise ValueError("No genes provided")
    fasta_keys = set(dataset.data.fa.keys())

    if isinstance(dataset.data.gtf, MockGTF):
        missing = [t for t in tokens if t not in fasta_keys]
        if missing:
            raise ValueError(
                f"Dataset has no GTF; targets must be FASTA record IDs. Not found: {missing[:10]}. "
                "Re-create the dataset with --gtf (or mkprobes ingest) to enable gene-name lookups."
            )
        return pl.DataFrame({
            "gene_name": tokens,
            "gene_id": tokens,
            "transcript_id": tokens,
            "transcript_name": tokens,
        })

    gtf = dataset.data.gtf
    cols = ["gene_name", "gene_id", "transcript_id", "transcript_name"]
    tx = gtf.filter(pl.col("feature") == "transcript")[cols] if "feature" in gtf.columns else gtf[cols]

    frames: list[pl.DataFrame] = []
    missing: list[str] = []
    for token in tokens:
        hit = tx.filter(pl.col("transcript_id") == token)
        if hit.is_empty():
            hit = tx.filter((pl.col("gene_name") == token) | (pl.col("gene_id") == token))
        if hit.is_empty():
            via = _resolve_via_annotations(dataset, token)
            if via:
                hit = tx.filter(
                    pl.col("transcript_id").is_in(via)
                    | pl.col("gene_id").is_in(via)
                    | pl.col("gene_name").is_in(via)
                )
        if hit.is_empty():
            missing.append(token)
            continue
        if mode == "longest" and hit["transcript_id"].n_unique() > 1:
            lengths = {
                tid: len(dataset.data.fa[tid]) if tid in fasta_keys else 0
                for tid in hit["transcript_id"].to_list()
            }
            per_gene = []
            for _, group in hit.group_by("gene_id", maintain_order=True):
                best = max(group["transcript_id"].to_list(), key=lambda t: lengths.get(t, 0))
                per_gene.append(group.filter(pl.col("transcript_id") == best))
            hit = pl.concat(per_gene)
        frames.append(hit.unique(subset=["transcript_id"], maintain_order=True))

    if missing:
        vocabulary = set(tx["gene_name"].drop_nulls()) | set(tx["gene_id"].drop_nulls()) | set(
            tx["transcript_id"].drop_nulls()
        )
        for token in missing:
            suggestions = difflib.get_close_matches(token, vocabulary, n=3, cutoff=0.7)
            logger.warning(
                f"Could not resolve {token!r}"
                + (f"; close matches: {suggestions}" if suggestions else " (no close matches)")
            )
        raise ValueError(
            f"{len(missing)}/{len(tokens)} targets could not be resolved: {missing[:10]}. "
            "Use transcript IDs / gene IDs from the annotation, or register an ortholog/alias "
            "annotation table (create-dataset/ingest --annotation-table)."
        )
    return pl.concat(frames)


def get_transcripts(
    dataset: Dataset,
    genes: Sequence[str],
    mode: Literal["gencode", "ensembl", "canonical", "appris", "apprisalt", "longest", "all"] = "canonical",
    output_path: Path | None = None,
    overwrite: bool = False,
) -> pl.DataFrame:
    """Get transcript ID from gene name or gene ID
    Returns:
        pl.DataFrame[[transcript_id, transcript_name, tag]]
        pl.DataFrame[[transcript_id, transcript_name, annotation, tag]] if appris
    """
    if not len(genes):
        raise ValueError("No genes provided")

    if not isinstance(dataset, ReferenceDataset):
        if mode == "canonical":
            logger.info("Custom dataset: using mode='longest' (canonical requires Ensembl).")
            mode = "longest"
        if mode not in GENERIC_MODES:
            raise ValueError(
                f"Mode {mode!r} requires a human/mouse reference dataset. "
                f"Custom datasets support: {GENERIC_MODES}."
            )
        return get_transcripts_generic(dataset, genes, mode)

    if mode in GENERIC_MODES:
        return get_transcripts_generic(dataset, genes, mode)  # type: ignore[arg-type]

    if not dataset.ensembl:
        raise ValueError("Not a ReferenceDataset. Cannot get transcripts.")

    df_genes = dataset.ensembl.filter(pl.col("gene_name").is_in(genes))[
        ["gene_name", "gene_id", "transcript_name", "transcript_id"]
    ]

    if dataset.appris is not None:
        df_genes = df_genes.join(
            dataset.appris.filter(pl.col("gene_id").is_in(df_genes["gene_id"]))[
                ["transcript_id", "annotation"]
            ],
            on="transcript_id",
            how="left",
        ).sort("transcript_name")

    to_return = ["gene_name", "gene_id", "transcript_id", "transcript_name", "tag"]

    match mode:
        case "canonical":
            with ThreadPoolExecutor(3) as exc:
                from functools import partial

                res = exc.map(
                    partial(get_ensembl, output_path or "output/", overwrite=overwrite), df_genes["gene_id"]
                )
                canonical = [r["canonical_transcript"].split(".")[0] for r in res]

            res = dataset.ensembl.filter(pl.col("transcript_id").is_in(canonical))[to_return]
        case "gencode":
            res = dataset.data.filter(pl.col("gene_id").is_in(df_genes["gene_id"]))[to_return]
        case "ensembl":
            res = dataset.ensembl.filter(pl.col("gene_id").is_in(df_genes["gene_id"]))[to_return]
        case "appris":
            if dataset.appris is None:
                raise ValueError("No APPRIS data found.")
            appris = df_genes.filter(pl.col("annotation").is_not_null())
            # if len(principal := appris.filter(pl.col("annotation").str.contains("PRINCIPAL"))):
            #     logger.info("Principal transcripts: " + "\n".join(principal["transcript_id"]))
            res = appris.join(dataset.ensembl[["transcript_id", "tag"]], on="transcript_id", how="left")

            def handle_transcripts(group: pl.DataFrame):
                if len(group) > 1:
                    canonical = group.filter(pl.col("tag") == "Ensembl_canonical")
                    if len(canonical) == 1:
                        return canonical
                    principal = group.filter(pl.col("annotation").str.contains("PRINCIPAL"))
                    if len(principal) == 1:
                        return principal
                return group

            res = res.group_by("gene_name").map_groups(handle_transcripts)

        case "apprisalt":
            if dataset.appris is None:
                raise ValueError("No APPRIS data found.")
            appris = df_genes.filter(pl.col("annotation").is_not_null())
            res = appris.join(dataset.ensembl[["transcript_id", "tag"]], on="transcript_id", how="left")
        case _:  # type: ignore
            raise ValueError(f"Unknown mode: {mode}")

    # for gene, tss in sorted(res.group_by("gene_name"), key=lambda x: x[0]):
    #     if len(tss) > 1:
    #         print(
    #             f"Multiple transcripts found for {gene[0]}. See https://useast.ensembl.org/Mouse/Search/Results?q={gene[0]};site=ensembl;facet_species=Mouse"
    #         )
    #         print(f"Please pick one: {tss.with_row_index()}.")
    #         picked = input("Enter the index of the correct transcript: ")
    #         out.append(tss[int(picked)])
    #     else:
    #         out.append(tss)
    # try:
    #     out = pl.concat(out)
    # except ValueError:
    #     raise ValueError(f"No transcripts found for {genes}")
    # return out
    if len(res) != len(genes):
        logger.warning(f"Found {len(res)} transcripts for {len(genes)} genes.")
        logger.warning(f"Missing genes: {', '.join(set(genes) - set(res['gene_name']))}")
    return res


@click.command()
@click.argument("path", type=click.Path(exists=True, dir_okay=True, file_okay=False, path_type=Path))
@click.argument("genes", type=click.Path(exists=True, dir_okay=False, file_okay=True, path_type=Path))
@click.option(
    "-m",
    "--mode",
    type=click.Choice([*REFERENCE_MODES, *GENERIC_MODES]),
    default="canonical",
)
def convert_to_transcripts(
    path: Path,
    genes: Path,
    mode: Literal["gencode", "ensembl", "canonical", "appris", "apprisalt", "longest", "all"] = "canonical",
):
    """Convert gene names to transcript IDs (canonical for reference datasets; longest/all otherwise)"""
    ds = load_dataset(path)
    del path
    gene_names = read_target_list(genes)
    res = get_transcripts(ds, gene_names, mode=mode)

    if mode != "all":
        res = res.group_by("gene_name", maintain_order=True).agg(pl.all().first())

    genes.with_suffix(".tss.txt").write_text("\n".join(sorted(res["transcript_name"])))
    if isinstance(ds, ReferenceDataset):
        # GENCODE-style names (Sox2-201): one transcript per base gene name.
        assert len(res["transcript_name"]) == len({
            x["transcript_name"].split("-")[0] for x in res.iter_rows(named=True)
        })
    logger.info(f"Written to {genes.with_suffix('.tss.txt')}.")


# fmt: off
@click.command()
@click.argument("path", type=click.Path(exists=True, dir_okay=True, file_okay=False, path_type=Path))
@click.option("--gene", type=str)
@click.option("--genefile", type=click.Path(exists=True, dir_okay=False, file_okay=True, path_type=Path))
@click.option("--canonical", "mode", flag_value="canonical", default=True, help="Outputs canonical transcript only (reference datasets; falls back to --longest for custom datasets)")
@click.option("--gencode"  , "mode", flag_value="gencode", help="Outputs all transcripts from GENCODE basic")
@click.option("--ensembl"  , "mode", flag_value="ensembl", help="Outputs all transcripts from Ensembl")
@click.option("--appris"   , "mode", flag_value="appris" , help="Outputs all principal transcripts from APPRIS (dominant coding transcripts)")
@click.option("--longest"  , "mode", flag_value="longest", help="Per gene, the isoform with the longest sequence (custom datasets; no network)")
@click.option("--all"      , "mode", flag_value="all"    , help="Every isoform of each gene (custom datasets; no network)")
@click.option("--verbose", "-v", is_flag=True, help="Verbose output")
# fmt: on
def transcripts(
    path: Path,
    gene: str | None = None,
    genefile: Path | None = None,
    mode: Literal["gencode", "ensembl", "canonical", "appris", "apprisalt", "longest", "all"] = "canonical",
    verbose: bool = False,
):
    """Get transcript ID from gene name or gene ID"""

    if genefile:
        genes = read_target_list(genefile)
    elif gene:
        genes = [gene]
    else:
        raise ValueError("No gene provided")
    res = get_transcripts(load_dataset(path), genes, mode)
    if verbose:
        click.echo(res)
    else:
        click.echo("\n".join(sorted(res["transcript_name"].to_list())))
