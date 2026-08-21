"""
Optional target selection: ``mkprobes suggest-targets``.

A panel has room for a fixed number of genes, and picking them by hand tends to
produce a list of markers that all report the same thing. Given expression data
and whatever genes you have already committed to, this suggests the genes that
add the most information those choices do not already carry, then scores how
much of the data's structure the resulting panel captures.

This is a *suggestion* step, not a design step. Nothing downstream requires it:
its output is an ordinary target list, and you are expected to edit it. It sits
before `chkgenes` in the workflow, and is skipped entirely if you already know
which genes you want.

The numerical work lives in `gene_selection`; this module is the command around
it, kept free of heavy imports so the CLI stays fast to start.
"""

from pathlib import Path

import rich_click as click
from loguru import logger

from .utils.targets import read_target_list

#: Below this share of top-PC variance, the suggestions are almost certainly
#: noise genes rather than genes reporting structure. Calibrated on synthetic
#: data with known programmes: properly filtered input scores ~67%, input
#: carrying unstructured genes scores ~1%.
LOW_SIGNAL_PERCENT = 5.0


def _load_expression(path: Path, layer: str | None):
    """Reads an AnnData file, with the heavy imports deferred to call time."""
    try:
        import anndata as ad
    except ImportError as e:  # pragma: no cover - anndata is a hard dependency
        raise click.ClickException(
            "Reading expression data needs `anndata`, which is missing from this "
            "environment. Reinstall mkprobes, or `pip install anndata`."
        ) from e

    try:
        adata = ad.read_h5ad(path)
    except Exception as e:
        raise click.ClickException(
            f"Could not read {path} as an AnnData (.h5ad) file: {e}\n"
            "Expression data must be samples (cells) as rows and genes as columns."
        ) from e

    if layer is not None and layer not in adata.layers:
        available = ", ".join(map(str, adata.layers.keys())) or "none"
        raise click.ClickException(
            f"{path} has no layer {layer!r}. Available layers: {available}. "
            "Omit --layer to use the main expression matrix."
        )
    return adata


def _check_known(adata, known: list[str], source: Path) -> None:
    """Every gene you already committed to has to be in the expression data."""
    present = set(adata.var_names)
    missing = [gene for gene in known if gene not in present]
    if missing:
        shown = ", ".join(missing[:10])
        raise click.ClickException(
            f"{len(missing)} gene(s) from {source} are not in the expression data: {shown}"
            f"{' ...' if len(missing) > 10 else ''}. Gene names have to match "
            "`adata.var_names` exactly - check for a different naming convention "
            "(symbols vs Ensembl IDs) or a species mismatch."
        )


@click.command("suggest-targets")  # fmt: off
@click.argument("expression", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option("--add", "-n", "n_add", type=int, required=True,
              help="How many genes to suggest, on top of any you already have.")
@click.option("--have", type=click.Path(exists=True, dir_okay=False, path_type=Path), default=None,
              help="Targets you have already committed to, one per line. Suggestions are chosen to "
                   "add information these do not already carry.")
@click.option("--layer", type=str, default=None,
              help="Expression layer to use. Defaults to the main matrix (adata.X).")
@click.option("--out", "-o", type=click.Path(dir_okay=False, path_type=Path), default=None,
              help="Write the combined target list here, ready for `mkprobes chkgenes`.")
@click.option("--score-pcs", type=int, default=30, show_default=True,
              help="How many principal components to score the panel against. 0 to skip scoring.")
@click.option("--n-components", type=int, default=None,
              help="Components of the residual expression to select against. Defaults to "
                   "max(--add + 20, 50). Lower it if your data has only a few distinct "
                   "programmes; the choice materially changes which genes come out.")
@click.option("--seed", type=int, default=42, show_default=True,
              help="Seed for the randomized SVD, so a run is reproducible.")
def suggest_targets(
    expression: Path,
    n_add: int,
    have: Path | None,
    layer: str | None,
    out: Path | None,
    score_pcs: int,
    n_components: int | None,
    seed: int,
):
    """Suggest panel genes that add information your current targets do not.

    EXPRESSION is an AnnData (`.h5ad`) file with samples as rows and genes as
    columns. Each suggested gene is chosen to be as independent as possible of
    the genes already selected, so the panel spans more of the biology rather
    than measuring the same axis repeatedly.

    **Filter to informative genes first.** A gene that correlates with nothing
    looks maximally independent to this method, so raw expression data yields
    suggestions that are merely noisy. Feed it highly variable genes - scanpy's
    `highly_variable_genes` is the usual route. The command warns when the
    output looks like this happened.

    Optional. Skip it if you already know which genes you want; the output is an
    ordinary target list you are expected to review and edit.
    """
    if n_add < 1:
        raise click.BadParameter("must be at least 1.", param_hint="--add")

    from .gene_selection import calculate_variance_capture_in_global_pcs, select_orthogonal_genes

    adata = _load_expression(expression, layer)
    logger.info(f"{expression}: {adata.n_obs} samples x {adata.n_vars} genes.")

    known = read_target_list(have) if have else []
    if known:
        _check_known(adata, known, have)
        logger.info(f"Holding {len(known)} target(s) from {have}.")

    available = adata.n_vars - len(known)
    if n_add > available:
        raise click.ClickException(
            f"Asked for {n_add} suggestions but only {available} gene(s) remain after the "
            f"{len(known)} you already have. Lower --add."
        )

    suggested = select_orthogonal_genes(
        adata,
        pre_selected_gene_names=known,
        k=n_add,
        layer=layer,
        n_components_rs_svd=n_components,
        random_state_rs_svd=seed,
    )
    if len(suggested) < n_add:
        logger.warning(
            f"Returned {len(suggested)} of the {n_add} requested: the residual expression has "
            "fewer independent dimensions left than that."
        )

    if score_pcs > 0:
        combined = known + suggested
        suggestions_only, _, _ = calculate_variance_capture_in_global_pcs(
            adata, suggested, n_global_pcs_target=score_pcs, layer=layer
        )
        after, _, _ = calculate_variance_capture_in_global_pcs(
            adata, combined, n_global_pcs_target=score_pcs, layer=layer
        )
        message = f"Panel of {len(combined)} captures {after:.1f}% of the variance in the top {score_pcs} PCs"
        if known:
            before, _, _ = calculate_variance_capture_in_global_pcs(
                adata, known, n_global_pcs_target=score_pcs, layer=layer
            )
            message += f" (up from {before:.1f}% with your {len(known)} alone)"
        logger.info(message + ".")

        if suggestions_only < LOW_SIGNAL_PERCENT:
            logger.warning(
                f"The suggested genes account for only {suggestions_only:.1f}% of the variance in "
                f"the top {score_pcs} PCs, which almost always means the expression data still "
                "contains many uninformative genes. A gene that correlates with nothing looks "
                "maximally independent to this method, so unfiltered input yields suggestions "
                "that are merely noisy rather than genes reporting real biology.\n"
                "Filter to highly variable genes first - scanpy's `highly_variable_genes` is the "
                "usual route - and run this again. Treat the list above as unreliable until you do."
            )

    click.echo("\n".join(suggested))
    if out:
        out.write_text("\n".join(known + suggested) + "\n")
        logger.info(f"{len(known) + len(suggested)} target(s) written to {out}.")
        click.echo(f"\nReview {out}, then: mkprobes chkgenes <dataset> {out}")
