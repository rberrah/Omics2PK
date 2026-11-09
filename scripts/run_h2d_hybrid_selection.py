#!/usr/bin/env python
"""H2d: hybrid selection, wired to the same cross-validation as H2a.

Frozen by `OMICS2PK_H2D_HYBRID_SELECTION_PRESPECIFICATION_20261108.md` v1.0.0
with erratum `OMICS2PK_H2D_ERRATUM_20261108.md` v1.0.2.

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
    "dependence_index_not_a_replication_count": 3.55,
    "dependence_index_meaning": (
        "A spectral index of dependence between endpoints, computed from the "
        "eigenvalue spread of their correlation matrix. It is NOT a count of "
        "independent replications and must never be described as one. It is "
        "reported so a tally of agreeing endpoints is read against the "
        "dependence among them."),
    "of_endpoints": 10,
    "same_enzyme_pair": list(SAME_ENZYME_PAIR),
    "same_enzyme_pair_correlation": 0.841,
    "counting_rule": (
        "The two substrates of one enzyme count as one endpoint. No text may "
        "say 'confirmed on nine endpoints'; the admissible form is nine "
        "endpoint measurements from a single cohort carrying a dependence "
        "index of about three and a half, none of which was used to design "
        "the rule."),
    "why_not_replication": (
        "Independent replication would need another cohort. Section 4.1 of the "
        "parent protocol records that no admissible public cohort supplies one "
        "for this axis."),
}


class H2dError(RuntimeError):
    pass


def run_endpoint(endpoint, data, universe, orders, *, folds, fold_label,
                 exploratory, log, pool_features=None):
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
    # Out-of-fold predictions are kept for EVERY arm, not only the hybrid. The
    # log loss is computed from predictions and observations; feeding it
    # absolute errors in place of predictions leaves the plain MAE correct and
    # makes the log loss meaningless, and the log loss is a mandatory
    # robustness co-criterion rather than a decoration.
    predictions = {name: blank() for name in
                   ("hybrid_at_budget", "blind_mechanistic") + STATISTICAL}
    predictions["statistical_at_matched_dimension"] = {
        arm: blank() for arm in STATISTICAL}
    retained_counts = {k: [] for k in sizes}
    abstentions = {k: 0 for k in sizes}
    agreement = {k: [] for k in sizes}
    demoted_examples = {k: [] for k in sizes}
    destroyed = {"blind": np.full((repeats, y.size), np.nan),
                 "hybrid": np.full((repeats, y.size), np.nan),
                 "univariate": np.full((repeats, y.size), np.nan),
                 "hybrid_kept_the_dead_gene": 0,
                 "univariate_kept_the_dead_gene": 0,
                 "folds_examined": 0}

    blind_order = orders[endpoint]["netpk_target_seeded"]
    target_gene = orders[endpoint]["target_gene_only"][0]

    # Probe-level traceability. The probe standing for a gene is chosen on
    # training donors only, so the choice can differ between folds and a
    # selection is not traceable from the gene symbol alone. Every feature any
    # arm actually used is checked against the qualified pool, and the features
    # behind the hybrid's own panels are recorded fold by fold.
    trace = {"features_by_symbol": {}, "selections_examined": 0,
             "outside_the_pool": {}, "decoy_selections": 0}

    def note(symbols, *, chosen, record):
        for symbol in symbols:
            feature = chosen.get(symbol)
            trace["selections_examined"] += 1
            if feature is None:
                continue
            if feature == H2A.DECOY_FEATURE_ID:
                trace["decoy_selections"] += 1
            elif pool_features is not None and feature not in pool_features:
                trace["outside_the_pool"][symbol] = feature
            if record:
                counts = trace["features_by_symbol"].setdefault(symbol, {})
                counts[feature] = counts.get(feature, 0) + 1

    started = time.time()
    for fold_number, (train, test) in enumerate(donor_folds):
        repeat = repeat_of[fold_number]
        chosen_feature = {}
        matrix, names, _ = H2A.fold_matrix(data, universe, train, keep,
                                           chosen=chosen_feature)
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
            note(kept, chosen=chosen_feature, record=True)
            retained_counts[k].append(len(kept))
            abstentions[k] += int(selection.abstained)
            if fold_number == 0:
                demoted_examples[k] = list(selection.demoted)[:10]

            predicted = predict(kept)
            predictions["hybrid_at_budget"][k][repeat, test] = predicted
            errors["hybrid_at_budget"][k][repeat, test] = np.abs(predicted - y[test])

            blind = [g for g in blind_order if g in position_of][:k]
            note(blind, chosen=chosen_feature, record=False)
            blind_predicted = predict(blind)
            predictions["blind_mechanistic"][k][repeat, test] = blind_predicted
            errors["blind_mechanistic"][k][repeat, test] = np.abs(
                blind_predicted - y[test])

            for arm in STATISTICAL:
                note(statistical[arm][:max(k, len(kept))],
                     chosen=chosen_feature, record=False)
                at_budget = predict(statistical[arm][:k])
                predictions[arm][k][repeat, test] = at_budget
                errors[arm][k][repeat, test] = np.abs(at_budget - y[test])
                # Matched dimension: the comparator is cut to what the hybrid
                # actually used, so neither side is credited for being smaller.
                matched = predict(statistical[arm][: len(kept)])
                predictions["statistical_at_matched_dimension"][arm][k][
                    repeat, test] = matched
                errors["statistical_at_matched_dimension"][arm][k][repeat, test] = (
                    np.abs(matched - y[test]))

            agreement[k].append(agreement_with(
                selection, statistical["univariate_spearman"][:k]))

        # Destroyed-signal check, the observation that produced this rule. The
        # hybrid must demote the dead gene where the blind arm cannot.
        if target_gene in position_of:
            note([target_gene], chosen=chosen_feature, record=True)
            wrecked = matrix.copy()
            column = position_of[target_gene]
            # Permute WITHIN each block. A single permutation over all donors
            # would move training values into the test block and back; it
            # cannot leak the endpoint, but it breaks the separation this
            # diagnostic is supposed to respect.
            shuffler = np.random.default_rng(H2A.SEED + 31337 + fold_number)
            for block_rows in (train, test):
                wrecked[block_rows, column] = matrix[
                    block_rows[shuffler.permutation(block_rows.size)], column]
            destroyed["blind"][repeat, test] = np.abs(
                predict([target_gene], block=wrecked) - y[test])
            wrecked_choice = select_hybrid(
                blind_order, position_of, wrecked[train], y[train],
                clinical[train], budget=1)
            destroyed["hybrid"][repeat, test] = np.abs(
                predict(list(wrecked_choice.retained), block=wrecked) - y[test])
            destroyed["hybrid_kept_the_dead_gene"] += int(
                target_gene in wrecked_choice.retained)
            # The prespecification's first success criterion compares the
            # hybrid against the univariate arm on this same wrecked matrix, so
            # the univariate arm has to be computed here too.
            strength = H2A._spearman_strength(wrecked[train], y[train])
            univariate_choice = str(names[int(np.argmax(strength))])
            destroyed["univariate"][repeat, test] = np.abs(
                predict([univariate_choice], block=wrecked) - y[test])
            destroyed["univariate_kept_the_dead_gene"] += int(
                univariate_choice == target_gene)
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
        "target_gene_available_in_folds": destroyed["folds_examined"],
        "losses_at_budget": {},
        "losses_at_matched_dimension": {},
        "retained": {},
        "agreement_with_univariate": {},
        "demoted_on_first_fold": demoted_examples,
        "probe_traceability": {
            "probe_choice_rule": (
                "The probe standing for a gene is chosen on training donors "
                "only, by the highest mean among the probes passing the "
                "missingness rule, so it can differ between folds. A selection "
                "is therefore not traceable from the gene symbol alone."),
            "qualified_pool_applied": pool_features is not None,
            "selections_checked": trace["selections_examined"],
            "selected_features_outside_the_pool": dict(
                sorted(trace["outside_the_pool"].items())),
            "decoy_selections": trace["decoy_selections"],
            "features_used_by_the_hybrid_and_the_target": {
                symbol: dict(sorted(counts.items()))
                for symbol, counts in sorted(trace["features_by_symbol"].items())},
        },
    }
    if pool_features is not None and trace["outside_the_pool"]:
        raise H2dError(
            f"{endpoint}: {len(trace['outside_the_pool'])} selected features "
            "are outside the qualified pool, so the restriction did not hold: "
            f"{sorted(trace['outside_the_pool'].items())[:3]}")
    if exploratory:
        summary["why_exploratory"] = (
            "This endpoint's destroyed-signal result generated the hybrid "
            "hypothesis, so the rule cannot be confirmed here whatever it "
            "produces. Reported for completeness and capped at exploratory.")

    hybrid_errors = {k: per_donor(errors["hybrid_at_budget"][k]) for k in sizes}
    for name in ("hybrid_at_budget", "blind_mechanistic") + STATISTICAL:
        summary["losses_at_budget"][name] = {
            str(k): H2A.loss_summary(per_donor(errors[name][k]), y,
                                     per_donor(predictions[name][k]))
            for k in sizes}
    for arm in STATISTICAL:
        summary["losses_at_matched_dimension"][arm] = {
            str(k): H2A.loss_summary(
                per_donor(errors["statistical_at_matched_dimension"][arm][k]), y,
                per_donor(predictions["statistical_at_matched_dimension"][arm][k]))
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

    if not destroyed["folds_examined"]:
        # Under a qualified pool an endpoint can lose its own transcript. The
        # diagnostic that generated this hypothesis then cannot be run, and
        # saying so is not the same as reporting that it passed.
        summary["destroyed_signal_repair"] = {
            "gene": target_gene,
            "not_assessable": (
                f"{target_gene} has no admissible feature in this universe, so "
                "its donor order cannot be permuted and the destroyed-signal "
                "diagnostic cannot be run for this endpoint. The endpoint is "
                "kept as a secondary analysis of other mechanisms; its own "
                "transcript is not a measurable predictor here."),
        }
    if destroyed["folds_examined"]:
        blind_errors = per_donor(destroyed["blind"])
        hybrid_errors_wrecked = per_donor(destroyed["hybrid"])
        univariate_errors = per_donor(destroyed["univariate"])
        against_blind = H2A.paired_contrast(
            hybrid_errors_wrecked, blind_errors, seed=H2A.SEED + 30000,
            label="hybrid_minus_blind_on_destroyed_signal")
        against_univariate = H2A.paired_contrast(
            hybrid_errors_wrecked, univariate_errors, seed=H2A.SEED + 30001,
            label="hybrid_minus_univariate_on_destroyed_signal")
        summary["destroyed_signal_repair"] = {
            "gene": target_gene,
            "folds_examined": destroyed["folds_examined"],
            "permutation": "within each block separately, per-fold seed",
            "blind_arm_mae": float(np.nanmean(blind_errors)),
            "hybrid_arm_mae": float(np.nanmean(hybrid_errors_wrecked)),
            "univariate_arm_mae": float(np.nanmean(univariate_errors)),
            "hybrid_kept_the_dead_gene_in_n_folds":
                destroyed["hybrid_kept_the_dead_gene"],
            "univariate_kept_the_dead_gene_in_n_folds":
                destroyed["univariate_kept_the_dead_gene"],
            "contrast_against_blind": against_blind,
            "contrast_against_univariate": against_univariate,
            # The prespecification's first success criterion needs both halves.
            "repair_criterion_met": bool(
                against_univariate["non_inferior_within_margin"]
                and against_blind["interval_excludes_zero"]
                and against_blind["delta_mae"] < 0.0),
            "how_to_read": (
                "The target gene's donor order is permuted within each block, "
                "so it keeps its place in the mechanistic ordering while "
                "carrying nothing. The blind arm must keep choosing it; the "
                "hybrid must demote it. The repair criterion needs BOTH halves: "
                "the hybrid within the frozen margin of the univariate arm, and "
                "its interval against the blind arm excluding zero on the "
                "favourable side. Either one alone is not the repair."),
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
    parser.add_argument("--probe-pool", default=None,
                        help="A qualified-probe-pool file. Restricts every "
                             "arm to the same analytically admissible probes, "
                             "exactly as the H2a sensitivity runs do, so the "
                             "two analyses share one universe.")
    parser.add_argument("--pool-name", default="A_qualified_human",
                        help="Which admissible set inside the pool file "
                             "(default: %(default)s).")
    parser.add_argument("--restrict-to-network", action="store_true",
                        help="Sensitivity analysis: let every arm select only "
                             "from genes the public network contains.")
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
    pool_features = None
    if args.probe_pool:
        pool_features = H2A.qualified_pool_features(args.probe_pool,
                                                    args.pool_name)
        universe = H2A.restrict_universe_to_qualified_pool(
            universe, data, args.probe_pool, args.pool_name)
        log(f"QUALIFIED POOL {args.pool_name}: "
            f"{universe['qualified_pool']['symbols_after']} symbols "
            f"(was {universe['qualified_pool']['symbols_before']}), "
            f"{universe['qualified_pool']['probes_dropped']} probes dropped, "
            "no technical control present")
    if args.restrict_to_network:
        universe = H2A.restrict_universe_to_network(universe, graph)
        log(f"SENSITIVITY: restricted to the network -> "
            f"{universe['candidates_after_restriction']} candidates "
            f"(was {universe['candidates_before_restriction']})")
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
        "erratum": "OMICS2PK_H2D_ERRATUM_20261108.md v1.0.2",
        "run_on": date.today().isoformat(),
        "endpoint_set": args.which,
        "endpoint_class": "functional_in_vitro",
        "endpoint_class_note": (
            "Measured microsomal activity, not clinical pharmacokinetics. The "
            "functional analogue of H2; never pooled with a clinical PK result "
            "and cannot promote the clinical claim above PARTIALLY_SUPPORTED."),
        "generating_endpoint": GENERATING_ENDPOINT,
        "validation_set_independence": VALIDATION_SET_INDEPENDENCE,
        "universe": {
            "candidate_genes": len(universe["by_symbol"]),
            "candidate_genes_on_network": len(universe["genes_on_network"]),
            "restricted_to_network": bool(args.restrict_to_network),
            "qualified_pool": universe.get("qualified_pool"),
            "same_pool_as_h2a": (
                "H2a and H2d are only comparable when both run on the same "
                "universe. The pool file and pool name recorded here must "
                "match the H2a run the comparison cites; an unrestricted H2d "
                "result is not comparable with a pool-restricted H2a result."),
        },
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
            exploratory=(endpoint == GENERATING_ENDPOINT), log=log,
            pool_features=pool_features)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n",
                   encoding="utf-8", newline="\n")
    log(f"wrote {out}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
