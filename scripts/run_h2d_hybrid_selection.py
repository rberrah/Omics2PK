#!/usr/bin/env python
"""H2d: hybrid selection, wired to the same cross-validation as H2a.

Frozen by `OMICS2PK_H2D_HYBRID_SELECTION_PRESPECIFICATION_20261108.md` v1.0.0
with erratum `OMICS2PK_H2D_ERRATUM_20261108.md` v1.0.1.

Knowledge proposes a shortlist from the outcome-blind mechanistic ordering; the
training donors then keep, demote or abstain. Everything else -- folds, donors,
predictor, panel-size grid, loss hierarchy, bootstrap, the one-feature-per-gene
universe, the missing-value rule, the baseline covariates -- is imported from
the H2a runner unchanged, so the only difference between the hybrid arm and the
blind mechanistic arm is the retention step.

Three numbers are reported together for every panel size, because the hybrid may
retain fewer genes than its budget and comparing it against arms that always
retain `k` would credit selection quality for dimension reduction:

  1. performance at budget `k`;
  2. the number actually retained, with the abstention count;
  3. performance at matched dimension, where the comparators are cut to the
     number the hybrid actually used in that fold.

Agreement with the univariate arm is reported as containment, Jaccard and exact
set equality, which answer different questions; the declared collapse verdict
uses exact equality at `k >= 5`.

The generating endpoint is run separately and its result is capped at
exploratory, because the rule was designed after seeing it.

Run from the repository root:

    python scripts/run_h2d_hybrid_selection.py \
        --hlc "<archive directory>" --knowledge "<knowledge directory>" \
        --out reports/data/h2d_confirmatory.json --set confirmatory
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from datetime import date
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

from omics2pk.hybrid_mechanistic_selection import (  # noqa: E402
    FROZEN as HYBRID_FROZEN,
    agreement_with,
    select_hybrid,
)


def load_h2a():
    """Import the H2a runner so every shared choice is literally the same code."""
    path = REPO_ROOT / "scripts" / "run_h2a_hlc_selection_comparison.py"
    spec = importlib.util.spec_from_file_location("h2a_runner", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


H2A = load_h2a()

# The endpoint whose result generated the hypothesis. It cannot confirm it.
GENERATING_ENDPOINT = "CYP2D6_activity"
CONFIRMATORY_ENDPOINTS = tuple(e for e in H2A.ENDPOINTS if e != GENERATING_ENDPOINT)

# Two substrates of one enzyme, correlated at 0.841 on these donors. They count
# as one endpoint in any tally of how many endpoints agreed.
SAME_ENZYME_PAIR = ("CYP3A4M_activity", "CYP3A4T_activity")

COLLAPSE_VERDICT_MIN_K = 5
COLLAPSE_THRESHOLD = 0.90
STATISTICAL = ("univariate_spearman", "lasso_entry_order", "xgboost_gain")

# Measured on the donors carrying expression data, before any hybrid run. Kept
# as a named constant so the caveat travels with every report and can be tested
# as a value rather than as a fragment of source text.
VALIDATION_SET_INDEPENDENCE = {
    "claim": ("validation on endpoints not used to design the rule, "
              "NOT independent replication"),
    "same_cohort": True,
    "pairwise_spearman_median": 0.446,
    "pairwise_spearman_range": [0.144, 0.841],
    "first_component_share_of_variance": 0.496,
    "effective_independent_endpoints": 3.55,
    "of_endpoints": 10,
    "same_enzyme_pair": list(SAME_ENZYME_PAIR),
    "same_enzyme_pair_correlation": 0.841,
    "counting_rule": (
        "The two substrates of one enzyme count as one endpoint. No text may "
        "say 'confirmed on nine endpoints'; the admissible form is nine "
        "endpoint measurements from a single cohort of effective independent "
        "dimensionality about three and a half, none of which was used to "
        "design the rule."),
    "why_not_replication": (
        "Independent replication would need another cohort. Section 4.1 of the "
        "parent protocol records that no admissible public cohort supplies one "
        "for this axis."),
}


class H2dError(RuntimeError):
    pass


def run_endpoint(endpoint, data, universe, orders, *, folds, fold_label,
                 exploratory, log):
    import pandas as pd

    y_all = pd.to_numeric(data["phenotype"].loc[endpoint],
                          errors="coerce").to_numpy(float)
    keep = np.nonzero(np.isfinite(y_all))[0]
    y = y_all[keep]
    if y.size < 50:
        raise H2dError(f"{endpoint}: only {y.size} donors with a measured value")

    donor_folds = [(np.searchsorted(keep, np.intersect1d(keep, train)),
                    np.searchsorted(keep, np.intersect1d(keep, test)))
                   for train, test in folds]
    donor_folds = [(tr, te) for tr, te in donor_folds if te.size and tr.size]
    repeats, repeat_of = H2A.repeat_assignment(
        len(donor_folds), n_splits=H2A.OUTER_FOLDS,
        repeated=fold_label.endswith("repeats"))

    sizes = H2A.PANEL_SIZES
    blank = lambda: {k: np.full((repeats, y.size), np.nan) for k in sizes}  # noqa: E731
    errors = {name: blank() for name in
              ("hybrid_at_budget", "blind_mechanistic") + STATISTICAL}
    errors["statistical_at_matched_dimension"] = {
        arm: blank() for arm in STATISTICAL}
    predictions = {"hybrid_at_budget": blank()}
    retained_counts = {k: [] for k in sizes}
    abstentions = {k: 0 for k in sizes}
    agreement = {k: [] for k in sizes}
    demoted_examples = {k: [] for k in sizes}
    destroyed = {"blind": np.full((repeats, y.size), np.nan),
                 "hybrid": np.full((repeats, y.size), np.nan),
                 "hybrid_kept_the_dead_gene": 0,
                 "folds_examined": 0}

    blind_order = orders[endpoint]["netpk_target_seeded"]
    target_gene = orders[endpoint]["target_gene_only"][0]

    started = time.time()
    for fold_number, (train, test) in enumerate(donor_folds):
        repeat = repeat_of[fold_number]
        matrix, names, _ = H2A.fold_matrix(data, universe, train, keep)
        position_of = {str(name): index for index, name in enumerate(names)}
        clinical = H2A.baseline_matrix(data, keep, train)
        statistical = H2A.statistical_orders(matrix[train], names, y[train],
                                             seed=H2A.SEED + fold_number)

        def predict(column_names, block=matrix):
            picked = [position_of[name] for name in column_names
                      if name in position_of]
            design = (np.column_stack([clinical, block[:, picked]])
                      if picked else clinical)
            return H2A.ridge_prediction(design[train], y[train], design[test],
                                        inner_folds=H2A.INNER_FOLDS,
                                        seed=H2A.SEED + fold_number)[0]

        for k in sizes:
            selection = select_hybrid(blind_order, position_of, matrix[train],
                                      y[train], clinical[train], budget=k)
            kept = list(selection.retained)
            retained_counts[k].append(len(kept))
            abstentions[k] += int(selection.abstained)
            if fold_number == 0:
                demoted_examples[k] = list(selection.demoted)[:10]

            predicted = predict(kept)
            predictions["hybrid_at_budget"][k][repeat, test] = predicted
            errors["hybrid_at_budget"][k][repeat, test] = np.abs(predicted - y[test])

            blind = [g for g in blind_order if g in position_of][:k]
            errors["blind_mechanistic"][k][repeat, test] = np.abs(
                predict(blind) - y[test])

            for arm in STATISTICAL:
                at_budget = statistical[arm][:k]
                errors[arm][k][repeat, test] = np.abs(
                    predict(at_budget) - y[test])
                # Matched dimension: the comparator is cut to what the hybrid
                # actually used, so neither side is credited for being smaller.
                matched = statistical[arm][: len(kept)]
                errors["statistical_at_matched_dimension"][arm][k][repeat, test] = (
                    np.abs(predict(matched) - y[test]))

            agreement[k].append(agreement_with(
                selection, statistical["univariate_spearman"][:k]))

        # Destroyed-signal check, the observation that produced this rule. The
        # hybrid must demote the dead gene where the blind arm cannot.
        if target_gene in position_of:
            wrecked = matrix.copy()
            shuffle = np.random.default_rng(
                H2A.SEED + 31337 + fold_number).permutation(y.size)
            wrecked[:, position_of[target_gene]] = matrix[
                shuffle, position_of[target_gene]]
            destroyed["blind"][repeat, test] = np.abs(
                predict([target_gene], block=wrecked) - y[test])
            wrecked_choice = select_hybrid(
                blind_order, position_of, wrecked[train], y[train],
                clinical[train], budget=1)
            destroyed["hybrid"][repeat, test] = np.abs(
                predict(list(wrecked_choice.retained), block=wrecked) - y[test])
            destroyed["hybrid_kept_the_dead_gene"] += int(
                target_gene in wrecked_choice.retained)
            destroyed["folds_examined"] += 1

        log(f"    fold {fold_number + 1}/{len(donor_folds)} "
            f"({time.time() - started:.0f}s)")

    per_donor = lambda block: np.nanmean(np.asarray(block, float), axis=0)  # noqa: E731
    summary = {
        "endpoint": endpoint,
        "status_cap": "exploratory" if exploratory else "confirmatory_on_unused_endpoint",
        "is_generating_endpoint": bool(exploratory),
        "fold_design": fold_label,
        "repeats": repeats,
        "n_donors": int(y.size),
        "folds": len(donor_folds),
        "target_gene": target_gene,
        "losses_at_budget": {},
        "losses_at_matched_dimension": {},
        "retained": {},
        "agreement_with_univariate": {},
        "demoted_on_first_fold": demoted_examples,
    }
    if exploratory:
        summary["why_exploratory"] = (
            "This endpoint's destroyed-signal result generated the hybrid "
            "hypothesis, so the rule cannot be confirmed here whatever it "
            "produces. Reported for completeness and capped at exploratory.")

    hybrid_errors = {k: per_donor(errors["hybrid_at_budget"][k]) for k in sizes}
    for name in ("hybrid_at_budget", "blind_mechanistic") + STATISTICAL:
        summary["losses_at_budget"][name] = {
            str(k): H2A.loss_summary(
                per_donor(errors[name][k]), y,
                per_donor(predictions["hybrid_at_budget"][k])
                if name == "hybrid_at_budget" else per_donor(errors[name][k]))
            for k in sizes}
    for arm in STATISTICAL:
        summary["losses_at_matched_dimension"][arm] = {
            str(k): {"mae": float(np.nanmean(per_donor(
                errors["statistical_at_matched_dimension"][arm][k])))}
            for k in sizes}

    for k in sizes:
        counts = np.asarray(retained_counts[k], dtype=float)
        identical = [a["identical"] for a in agreement[k] if a["identical"] is not None]
        summary["retained"][str(k)] = {
            "budget": k,
            "median_retained": float(np.median(counts)),
            "mean_retained": float(counts.mean()),
            "min_retained": int(counts.min()),
            "max_retained": int(counts.max()),
            "abstentions": abstentions[k],
            "folds": len(counts),
        }
        summary["agreement_with_univariate"][str(k)] = {
            "mean_containment": float(np.nanmean(
                [a["containment"] for a in agreement[k]])),
            "mean_jaccard": float(np.nanmean([a["jaccard"] for a in agreement[k]])),
            "identical_set_fraction": (float(np.mean(identical))
                                       if identical else None),
            "folds_with_a_panel_to_compare": len(identical),
        }

    verdict = {k: summary["agreement_with_univariate"][str(k)]["identical_set_fraction"]
               for k in sizes if k >= COLLAPSE_VERDICT_MIN_K}
    assessable = [v for v in verdict.values() if v is not None]
    summary["collapse_verdict"] = {
        "assessed_at_panel_sizes": [k for k in sizes if k >= COLLAPSE_VERDICT_MIN_K],
        "threshold": COLLAPSE_THRESHOLD,
        "measure": "exact set equality with the univariate top-k",
        "collapsed": bool(assessable and max(assessable) >= COLLAPSE_THRESHOLD),
        "worst_identical_fraction": (max(assessable) if assessable else None),
        "why_small_k_is_excluded": (
            "At k = 1 every arm already chose the same single transcript in the "
            "H2a run, so exact equality there is expected and says nothing "
            "about collapse. All three measures are still reported at every k."),
    }

    contrasts = []
    for k in sizes:
        contrasts.append(H2A.paired_contrast(
            hybrid_errors[k], per_donor(errors["blind_mechanistic"][k]),
            seed=H2A.SEED + 9000 + k,
            label=f"hybrid_minus_blind_mechanistic_at_budget_k{k}"))
        for offset, arm in enumerate(STATISTICAL):
            contrasts.append(H2A.paired_contrast(
                hybrid_errors[k], per_donor(errors[arm][k]),
                seed=H2A.SEED + 10000 + 1000 * offset + k,
                label=f"hybrid_minus_{arm}_at_budget_k{k}"))
            contrasts.append(H2A.paired_contrast(
                hybrid_errors[k],
                per_donor(errors["statistical_at_matched_dimension"][arm][k]),
                seed=H2A.SEED + 20000 + 1000 * offset + k,
                label=f"hybrid_minus_{arm}_at_matched_dimension_k{k}"))
    summary["contrasts"] = contrasts
    summary["non_inferiority_against_every_statistical_arm"] = {
        str(k): {
            "at_budget": all(
                next(c for c in contrasts
                     if c["label"] == f"hybrid_minus_{arm}_at_budget_k{k}"
                     )["non_inferior_within_margin"] for arm in STATISTICAL),
            "at_matched_dimension": all(
                next(c for c in contrasts
                     if c["label"] ==
                     f"hybrid_minus_{arm}_at_matched_dimension_k{k}"
                     )["non_inferior_within_margin"] for arm in STATISTICAL),
        } for k in sizes}

    if destroyed["folds_examined"]:
        blind_mae = float(np.nanmean(per_donor(destroyed["blind"])))
        hybrid_mae = float(np.nanmean(per_donor(destroyed["hybrid"])))
        summary["destroyed_signal_repair"] = {
            "gene": target_gene,
            "folds_examined": destroyed["folds_examined"],
            "blind_arm_mae": blind_mae,
            "hybrid_arm_mae": hybrid_mae,
            "hybrid_kept_the_dead_gene_in_n_folds":
                destroyed["hybrid_kept_the_dead_gene"],
            "contrast": H2A.paired_contrast(
                per_donor(destroyed["hybrid"]), per_donor(destroyed["blind"]),
                seed=H2A.SEED + 30000,
                label="hybrid_minus_blind_on_destroyed_signal"),
            "how_to_read": (
                "The target gene's donor order is permuted, so it keeps its "
                "place in the mechanistic ordering while carrying nothing. The "
                "blind arm must keep choosing it; the hybrid must demote it. A "
                "negative contrast excluding zero is the repair this rule "
                "exists to deliver."),
        }
    return summary


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hlc", required=True)
    parser.add_argument("--knowledge", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--set", dest="which", default="confirmatory",
                        choices=("confirmatory", "generating", "custom"))
    parser.add_argument("--endpoints", nargs="*", default=None)
    parser.add_argument("--design", choices=("repeated", "centre"),
                        default="repeated")
    args = parser.parse_args(argv)

    def log(message):
        print(message, flush=True)

    if args.which == "confirmatory":
        endpoints = list(CONFIRMATORY_ENDPOINTS)
    elif args.which == "generating":
        endpoints = [GENERATING_ENDPOINT]
    else:
        if not args.endpoints:
            raise H2dError("--set custom needs --endpoints")
        endpoints = list(args.endpoints)

    knowledge = Path(args.knowledge)
    graph = H2A.load_string_symbol_graph(
        knowledge / "9606.protein.info.v12.0.txt.gz",
        knowledge / "9606.protein.links.v12.0.txt.gz")
    data = H2A.load_hlc(Path(args.hlc))
    universe = H2A.build_universe(data, graph, np.random.default_rng(H2A.SEED))
    orders = H2A.mechanistic_orders(graph, knowledge)
    log(f"network {graph.n_nodes} nodes; candidates "
        f"{len(universe['by_symbol'])}; endpoints {len(endpoints)}")

    if args.design == "centre":
        folds, fold_label = H2A.centre_folds(data["centre"]), "leave_one_centre_out"
    else:
        folds = H2A.repeated_folds(data["values"].shape[1],
                                   n_splits=H2A.OUTER_FOLDS,
                                   n_repeats=H2A.OUTER_REPEATS, seed=H2A.SEED)
        fold_label = f"{H2A.OUTER_FOLDS}_fold_x{H2A.OUTER_REPEATS}_repeats"

    report = {
        "analysis": "H2d: hybrid mechanistic selection",
        "prespecification": "OMICS2PK_H2D_HYBRID_SELECTION_PRESPECIFICATION_20261108.md v1.0.0",
        "erratum": "OMICS2PK_H2D_ERRATUM_20261108.md v1.0.1",
        "run_on": date.today().isoformat(),
        "endpoint_set": args.which,
        "endpoint_class": "functional_in_vitro",
        "endpoint_class_note": (
            "Measured microsomal activity, not clinical pharmacokinetics. The "
            "functional analogue of H2; never pooled with a clinical PK result "
            "and cannot promote the clinical claim above PARTIALLY_SUPPORTED."),
        "generating_endpoint": GENERATING_ENDPOINT,
        "validation_set_independence": VALIDATION_SET_INDEPENDENCE,
        "frozen": {
            "hybrid": dict(HYBRID_FROZEN),
            "inherited_from_h2a": {
                "seed": H2A.SEED, "outer_folds": H2A.OUTER_FOLDS,
                "outer_repeats": H2A.OUTER_REPEATS,
                "inner_folds": H2A.INNER_FOLDS,
                "panel_sizes": list(H2A.PANEL_SIZES),
                "bootstrap": H2A.N_BOOTSTRAP,
                "non_inferiority_margin": H2A.NON_INFERIORITY_MARGIN,
                "baseline_covariates": ["age", "sex"],
                "fold_design": fold_label,
            },
            "collapse_verdict": {"threshold": COLLAPSE_THRESHOLD,
                                 "min_panel_size": COLLAPSE_VERDICT_MIN_K,
                                 "measure": "exact set equality"},
        },
        "endpoints": {},
    }

    for endpoint in endpoints:
        log(f"  {endpoint}")
        report["endpoints"][endpoint] = run_endpoint(
            endpoint, data, universe, orders, folds=folds,
            fold_label=fold_label,
            exploratory=(endpoint == GENERATING_ENDPOINT), log=log)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n",
                   encoding="utf-8", newline="\n")
    log(f"wrote {out}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
