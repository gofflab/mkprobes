"""Quantify the panel-level impact of the two probe-selection defects.

Both defects are inherited verbatim from the upstream ``fishtools`` code and both
change which probes end up in a panel, so this script measures the change rather
than applying it.

Defect 1 -- opposite sign conventions for ``overlap``
    ``OverlapWeighted.q``  uses ``start[j - 1] + overlap``  (a negative value
    forces a gap), while ``find_overlap`` uses ``start[i] - overlap`` (a negative
    value permits an overlap). ``find_overlap`` is used for priority tier 1 only.

Defect 2 -- ``selected_global`` accumulated outside its loop
    ``handle_overlap`` computes ``sel_local`` per priority tier but applies
    ``selected_global |= sel_local`` after the loop, so only the last tier's
    selection survives and the in-loop ``~is_in(selected_global)`` filter is a
    no-op.

Usage
-----
    python scripts/quantify_overlap_defects.py DATA_DIR [-o report.json] [--overlap -2]

``DATA_DIR`` is a directory of ``<gene>_crawled.parquet`` files, e.g.
``data/och_test_output``. The script reports, per gene and per variant, how many
probe pairs are selected and how many differ from the current behaviour.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any, cast

import numpy as np
import polars as pl

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from mkprobes.utils._algorithms import find_overlap, find_overlap_weighted
from mkprobes.utils._filtration import filter_have_both, pair_name
from mkprobes.utils._filtration import the_filter as the_filter_real

# --------------------------------------------------------------------------- #
# Selector variants
# --------------------------------------------------------------------------- #


def find_overlap_fixed(start: Sequence[int], end: Sequence[int], overlap: int = 0) -> list[int]:
    """``find_overlap`` with the sign convention aligned to ``OverlapWeighted.q``.

    ``OverlapWeighted`` treats probe ``i`` as compatible with ``j`` when
    ``end[i] < start[j] + overlap``. This is the same predicate, applied greedily.
    """
    out = [0]
    curr_end = end[0]
    for i in range(1, len(start)):
        if end[i] < curr_end:
            raise ValueError("Ends not sorted")
        if start[i] + overlap > end[out[-1]]:
            out.append(i)
    return out


# --------------------------------------------------------------------------- #
# Instrumented handle_overlap
# --------------------------------------------------------------------------- #


def handle_overlap_variant(
    df: pl.DataFrame,
    criteria: list[pl.Expr],
    overlap: int = -2,
    n: int = 100,
    *,
    fix_sign: bool = False,
    loop_mode: str = "replace",
) -> tuple[pl.DataFrame, dict[str, Any]]:
    """Faithful re-implementation of ``handle_overlap`` with the defects switchable.

    fix_sign
        ``False`` keeps ``find_overlap``'s ``start[i] - overlap`` convention.
        ``True`` uses ``start[i] + overlap``, matching the weighted selector.
    loop_mode
        ``"replace"``   -- current behaviour: the last tier's ``sel_local`` wins.
        ``"accumulate"`` -- ``selected_global |= sel_local`` moved inside the loop,
                            i.e. the reading in which each tier tops up the panel.
    """
    if df.select(pl.col("gene").n_unique()).item() > 1:
        raise ValueError("More than one gene in filtered")

    df = df.sort(by=["pos_end", "tm"], descending=[False, True])
    criteria = criteria or [pl.col("*")]

    df_ori = df.with_row_index("index")

    ddf = df_ori.lazy().with_columns(priority=pl.lit(0, dtype=pl.UInt8))
    for priority, criterion in reversed(list(enumerate(criteria, 1))):
        ddf = ddf.update(
            ddf.filter(criterion).with_columns(priority=pl.lit(priority, dtype=pl.UInt8)),
            on="index",
        )
    ddf = ddf.collect()
    df = ddf.filter(pl.col("priority") > 0)
    df = filter_have_both(df)

    stats: dict[str, Any] = {"start": len(df), "match_any": len(df)}
    selected_global: set[int] = set()
    if not len(df):
        return df, stats | {"selected_pair": 0, "empty": True}

    greedy = find_overlap_fixed if fix_sign else find_overlap
    sel_local: set[int] = set()
    tiers_used: list[int] = []

    for i in range(1, len(criteria) + 1):
        run = (
            df.filter((pl.col("priority") <= i) & ~pl.col("index").is_in(selected_global))
            .filter(pl.col("name").str.ends_with("splint"))
            .select(["index", "pos_start", "pos_end", "priority"])
            .sort(["pos_end", "pos_start"])
        )
        if not len(run):
            continue

        priorities = np.sqrt((len(criteria) + 1 - run["priority"]).cast(pl.Float64).to_numpy())
        try:
            if i == 1:
                ols = greedy(
                    cast(Sequence[int], run["pos_start"]),
                    cast(Sequence[int], run["pos_end"]),
                    overlap=overlap,
                )
            else:
                ols = find_overlap_weighted(
                    cast(Sequence[int], run["pos_start"]),
                    cast(Sequence[int], run["pos_end"]),
                    cast(Sequence[int], priorities),
                    overlap=overlap,
                )
            sel_local = set(run[ols]["index"].to_list())
            tiers_used.append(i)
            stats[f"selected_{i}"] = len(sel_local)
            if loop_mode == "accumulate":
                selected_global |= sel_local
                if len(selected_global) > n:
                    break
            elif len(sel_local) > n:
                break
        except RecursionError:
            stats["recursion_error"] = i
            break

    if loop_mode != "accumulate":
        selected_global |= sel_local

    stats["tiers_used"] = tiers_used
    stats["break_tier"] = tiers_used[-1] if tiers_used else None
    selected_names = df.filter(pl.col("index").is_in(selected_global)).select(name=pair_name)["name"]
    out = df.filter(pair_name.is_in(selected_names.to_list()))
    stats["selected_pair"] = len(out) // 2
    return out, stats


def build_criteria(max_tm_offtarget: float = 20.0, max_hp: float = 32.0) -> list[pl.Expr]:
    """The exact criteria ladder from ``the_filter``."""
    clean = pl.col("maps_to_pseudo").is_null() | pl.col("maps_to_pseudo").eq("")
    return [
        (pl.col("oks") > 5) & (pl.col("hp") < max_hp) & pl.col("max_tm_offtarget").lt(max_tm_offtarget) & clean,
        (pl.col("oks") > 4) & (pl.col("hp") < max_hp) & pl.col("max_tm_offtarget").lt(max_tm_offtarget) & clean,
        (pl.col("oks") > 4) & (pl.col("hp") < max_hp) & pl.col("max_tm_offtarget").lt(max_tm_offtarget),
        (pl.col("oks") > 4) & (pl.col("hp") < max_hp + 5) & pl.col("max_tm_offtarget").lt(max_tm_offtarget + 4),
        (pl.col("oks") > 3) & (pl.col("hp") < max_hp + 5) & pl.col("max_tm_offtarget").lt(max_tm_offtarget + 4),
        (pl.col("oks") > 2) & (pl.col("hp") < max_hp + 5) & pl.col("max_tm_offtarget").lt(max_tm_offtarget + 4),
        (pl.col("oks") > 1) & (pl.col("hp") < max_hp + 5) & pl.col("max_tm_offtarget").lt(max_tm_offtarget + 4),
    ]


VARIANTS: dict[str, dict[str, Any]] = {
    "baseline": dict(fix_sign=False, loop_mode="replace"),
    "fix1_sign": dict(fix_sign=True, loop_mode="replace"),
    "fix2_accumulate": dict(fix_sign=False, loop_mode="accumulate"),
    "fix1+fix2": dict(fix_sign=True, loop_mode="accumulate"),
}


def selected_pair_names(df: pl.DataFrame) -> set[str]:
    if not len(df):
        return set()
    return set(df.select(name=pair_name)["name"].unique().to_list())


def run_gene(df: pl.DataFrame, overlap: int, n: int) -> dict[str, Any]:
    criteria = build_criteria()
    res: dict[str, Any] = {}
    names: dict[str, set[str]] = {}
    for label, kw in VARIANTS.items():
        out, stats = handle_overlap_variant(df, criteria, overlap=overlap, n=n, **kw)
        names[label] = selected_pair_names(out)
        res[label] = {
            "n_pairs": stats["selected_pair"],
            "break_tier": stats.get("break_tier"),
            "tiers_used": stats.get("tiers_used"),
        }

    base = names["baseline"]
    for label in VARIANTS:
        if label == "baseline":
            continue
        res[label]["added"] = len(names[label] - base)
        res[label]["dropped"] = len(base - names[label])
        res[label]["jaccard"] = (
            round(len(base & names[label]) / len(base | names[label]), 4) if (base | names[label]) else 1.0
        )
    return res


def verify_baseline_matches_real(df: pl.DataFrame, overlap: int) -> bool | None:
    """Confirm the instrumented baseline reproduces the shipped ``the_filter``.

    Returns ``None`` if the shipped code bailed out (it calls ``exit(0)`` when no
    probe passes any filter), so the check can be retried on the next gene.
    """
    try:
        real, _ = the_filter_real(df, overlap=overlap)
    except SystemExit:
        return None
    mine, _ = handle_overlap_variant(df, build_criteria(), overlap=overlap, **VARIANTS["baseline"])
    return selected_pair_names(real) == selected_pair_names(mine)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_dir", type=Path, help="directory containing <gene>_crawled.parquet files")
    ap.add_argument("--overlap", type=int, default=-2, help="overlap parameter (production default: -2)")
    ap.add_argument("--n", type=int, default=100, help="handle_overlap early-stop threshold")
    ap.add_argument("-o", "--out", type=Path, default=None, help="write the full report as JSON")
    ap.add_argument("--no-verify", action="store_true", help="skip the baseline-fidelity check")
    args = ap.parse_args()

    files = sorted(args.data_dir.glob("*_crawled.parquet"))
    if not files:
        print(f"No *_crawled.parquet under {args.data_dir}", file=sys.stderr)
        return 1

    report: dict[str, Any] = {}
    verified = None
    header = f"{'gene':<24}{'base':>7}{'fix1':>7}{'fix2':>7}{'both':>7}   {'tier':>4}  {'fix1 +/-':>11}  {'fix2 +/-':>11}"
    print(header)
    print("-" * len(header))
    for f in files:
        gene = f.name.removesuffix("_crawled.parquet")
        df = pl.read_parquet(f)
        if not args.no_verify and verified is None:
            verified = verify_baseline_matches_real(df, args.overlap)
            if verified is False:
                print("WARNING: instrumented baseline does not match the shipped the_filter", file=sys.stderr)
        try:
            r = run_gene(df, overlap=args.overlap, n=args.n)
        except SystemExit:
            print(f"{gene:<24}  no probes passed filters")
            continue
        report[gene] = r
        print(
            f"{gene:<24}{r['baseline']['n_pairs']:>7}{r['fix1_sign']['n_pairs']:>7}"
            f"{r['fix2_accumulate']['n_pairs']:>7}{r['fix1+fix2']['n_pairs']:>7}"
            f"   {r['baseline']['break_tier']!s:>4}"
            f"  {'+' + str(r['fix1_sign']['added']) + '/-' + str(r['fix1_sign']['dropped']):>11}"
            f"  {'+' + str(r['fix2_accumulate']['added']) + '/-' + str(r['fix2_accumulate']['dropped']):>11}"
        )

    changed1 = [g for g, r in report.items() if r["fix1_sign"]["added"] or r["fix1_sign"]["dropped"]]
    changed2 = [g for g, r in report.items() if r["fix2_accumulate"]["added"] or r["fix2_accumulate"]["dropped"]]
    tier1 = [g for g, r in report.items() if r["baseline"]["break_tier"] == 1]
    print("-" * len(header))
    print(f"genes: {len(report)}")
    print(
        f"  break at tier 1 (greedy find_overlap reaches the panel; defect 1 is live): "
        f"{len(tier1)}/{len(report)}"
    )
    print(f"  panel changed by fix1 (sign convention):     {len(changed1)}/{len(report)}")
    print(f"  panel changed by fix2 (accumulate reading):  {len(changed2)}/{len(report)}")
    if changed1:
        net = sum(report[g]["fix1_sign"]["n_pairs"] - report[g]["baseline"]["n_pairs"] for g in changed1)
        print(f"  fix1 net probe-pair change across those genes: {net:+d}")
    if verified is not None:
        print(f"baseline fidelity vs shipped the_filter: {'OK' if verified else 'MISMATCH'}")

    if args.out:
        args.out.write_text(json.dumps({"overlap": args.overlap, "n": args.n, "genes": report}, indent=2))
        print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
