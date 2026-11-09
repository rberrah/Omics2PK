#!/usr/bin/env python
"""Name the feature a selector falls back on when the target signal is dead.

Why this exists
---------------
The destroyed-signal diagnostic of protocol §5.5 permutes the donor order of
the one column standing for the target gene, then reports how well each arm
predicts. The number alone cannot be interpreted, because three very different
things produce a *lower* error for an arm that reads the endpoint, and this
cohort contains all three:

* a **non-human duplicate** of the destroyed transcript, which a qualified pool
  refuses -- on `CYP2D6_activity` the feature annotated `CYP2D17`, chosen in ten
  folds out of ten;
* an **admissible paralogue** whose probe plausibly cross-hybridises the same
  transcript, which no annotation can expose as a duplicate because the
  paralogue is a real human gene with an identity of its own -- on
  `CYP2A6_activity` the feature annotated `CYP2A7`;
* a **different gene** whose association with the endpoint may be real biology
  -- on `CYP3A4M_activity` a feature whose accessions resolve to `IL12B-AS1`.

Reading the diagnostic as "the method resisted a destroyed mechanism" is
therefore only admissible once the replacement is named. This script names it,
per fold, with the verdict of the probe-identity qualification, and changes no
arm, loss, margin or contrast: it reports a choice the reported runs already
made internally.

    python scripts/name_destroyed_signal_replacement.py \\
        --hlc "<archive dir>" --knowledge "<knowledge dir>" \\
        --pools registry/QUALIFIED_PROBE_POOLS.json \\
        --out reports/data/h2a_destroyed_signal_replacements.json
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from datetime import date
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

from omics2pk.mechanistic_graph_selection import (  # noqa: E402
    seeds_from_endpoint_label,
)

REFUSED_STRATA = (
    ("undetermined", "undetermined_features"),
    ("conflicting", "conflicting_features"),
    ("symbol_names_another_gene", "symbol_names_another_gene_features"),
    ("platform_self_contradictory", "platform_self_contradictory_features"),
)


class ReplacementError(RuntimeError):
    pass


def load_h2a():
    path = REPO_ROOT / "scripts" / "run_h2a_hlc_selection_comparison.py"
    spec = importlib.util.spec_from_file_location("h2a_runner", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def verdict_index(pools):
    """feature id -> which stratum of the qualification it belongs to."""
    index = {}
    for feature in pools["admissible_features"]["A_qualified_human"]:
        index[feature] = "admissible_pool_a"
    for name, key in REFUSED_STRATA:
        for feature in pools.get(key, ()):
            index.setdefault(feature, name)
    return index


def name_replacements(H2A, data, universe, endpoint, graph, *, verdicts,
                      pool_name):
    """For each fold, the feature the univariate arm picks once the target is dead.

    Everything is taken from the H2a runner: the folds, the seeds, the matrix
    assembly, the permutation and the tie-aware rank statistic.
    """
    import pandas as pd

    target = next(iter(seeds_from_endpoint_label(endpoint, graph)))
    y_all = pd.to_numeric(data["phenotype"].loc[endpoint],
                          errors="coerce").to_numpy(float)
    keep = np.nonzero(np.isfinite(y_all))[0]
    y = y_all[keep]
    folds = H2A.repeated_folds(data["values"].shape[1],
                               n_splits=H2A.OUTER_FOLDS,
                               n_repeats=H2A.OUTER_REPEATS, seed=H2A.SEED)
    donor_folds = [(np.searchsorted(keep, np.intersect1d(keep, tr)),
                    np.searchsorted(keep, np.intersect1d(keep, te)))
                   for tr, te in folds]
    donor_folds = [(tr, te) for tr, te in donor_folds if tr.size and te.size]

    per_fold, counts = [], {}
    for fold_number, (train, _test) in enumerate(donor_folds):
        chosen_feature = {}
        matrix, names, is_decoy = H2A.fold_matrix(data, universe, train, keep,
                                                  chosen=chosen_feature)
        position_of = {str(name): index for index, name in enumerate(names)}
        if target not in position_of:
            return {"endpoint": endpoint, "target_gene": target,
                    "pool": pool_name,
                    "not_assessable": (
                        f"{target} has no feature in this universe, so its "
                        "donor order cannot be permuted")}
        wrecked = matrix.copy()
        shuffle = np.random.default_rng(
            H2A.SEED + 31337 + fold_number).permutation(y.size)
        wrecked[:, position_of[target]] = matrix[shuffle, position_of[target]]
        strength = H2A._spearman_strength(wrecked[train], y[train])
        index = int(np.argmax(strength))
        symbol = str(names[index])
        feature = chosen_feature.get(symbol, "")
        stratum = ("decoy" if bool(is_decoy[index])
                   else verdicts.get(feature, "unclassified"))
        per_fold.append({"fold": fold_number + 1, "symbol": symbol,
                         "feature_id": feature,
                         "abs_spearman": round(float(strength[index]), 4),
                         "identity_stratum": stratum,
                         "is_the_target_gene": symbol == target})
        key = (symbol, feature, stratum)
        counts[key] = counts.get(key, 0) + 1

    return {
        "endpoint": endpoint,
        "target_gene": target,
        "pool": pool_name,
        "folds": len(per_fold),
        "per_fold": per_fold,
        "by_replacement": [
            {"symbol": symbol, "feature_id": feature,
             "identity_stratum": stratum, "folds": n}
            for (symbol, feature, stratum), n in
            sorted(counts.items(), key=lambda item: -item[1])],
        "picked_the_dead_gene_in_n_folds": sum(
            row["is_the_target_gene"] for row in per_fold),
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hlc", required=True)
    parser.add_argument("--knowledge", required=True)
    parser.add_argument("--pools", default="registry/QUALIFIED_PROBE_POOLS.json")
    parser.add_argument("--probe-pool", default=None,
                        help="Restrict every candidate to this pool, as a "
                             "reported run would.")
    parser.add_argument("--pool-name", default="A_qualified_human")
    parser.add_argument("--endpoints", nargs="*", default=None)
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)

    H2A = load_h2a()
    knowledge = Path(args.knowledge)
    graph = H2A.load_string_symbol_graph(
        knowledge / "9606.protein.info.v12.0.txt.gz",
        knowledge / "9606.protein.links.v12.0.txt.gz")
    data = H2A.load_hlc(Path(args.hlc))
    universe = H2A.build_universe(data, graph, np.random.default_rng(H2A.SEED))
    pool_name = "unrestricted"
    if args.probe_pool:
        universe = H2A.restrict_universe_to_qualified_pool(
            universe, data, args.probe_pool, args.pool_name)
        pool_name = args.pool_name
    verdicts = verdict_index(json.loads(
        Path(args.pools).read_text(encoding="utf-8")))

    endpoints = args.endpoints or list(H2A.ENDPOINTS)
    report = {
        "report": "What each selector falls back on when the target signal is dead",
        "run_on": date.today().isoformat(),
        "protocol_version": "1.4.7",
        "universe": pool_name,
        "candidate_genes": len(universe["by_symbol"]),
        "why": ("The destroyed-signal number cannot be read without the name of "
                "the replacement: a non-human duplicate, an admissible "
                "paralogue of the same transcript and a genuinely different "
                "gene all produce a lower error for an arm that reads the "
                "endpoint, and this cohort contains all three."),
        "what_this_does_not_establish": (
            "A named replacement is not a demonstrated mechanism, and a high "
            "rank correlation is not a causal relation. Cross-hybridisation is "
            "never checked here, because no source distributes the "
            "oligonucleotide sequence; a paralogue is called a plausible "
            "duplicate, never a confirmed one."),
        "endpoints": {},
    }
    for endpoint in endpoints:
        print(f"  {endpoint}", flush=True)
        report["endpoints"][endpoint] = name_replacements(
            H2A, data, universe, endpoint, graph, verdicts=verdicts,
            pool_name=pool_name)
        for row in report["endpoints"][endpoint].get("by_replacement", []):
            print(f"    {row['folds']:2d} folds  {row['symbol']:16s} "
                  f"{row['feature_id']:14s} {row['identity_stratum']}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n",
                   encoding="utf-8", newline="\n")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
