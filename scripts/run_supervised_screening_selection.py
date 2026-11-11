#!/usr/bin/env python
"""Does the methodology choose better biological information? Two cohorts, one harness.

    python -I scripts/run_supervised_screening_selection.py \
        --mavedb-dir "<dir>" --pdb-dir "<dir>" \
        --hlc-expression "<dir>/curatedExpressionLiver.zip" \
        --hlc-phenotype "<dir>/curatedPhenotypesynapse.zip"

Either cohort may be omitted; what is absent is reported `NOT_EVALUABLE` with
the reason, never silently skipped.

The two questions the commissioning note separates
--------------------------------------------------
**After the measurements** — which omic variables genuinely add information to
what is already in hand? Answered by giving every selector the same pool, the
same folds and the same estimator, with the information already in hand pinned
into every panel as the baseline, and settling value by out-of-sample ablation.

**Before the measurements** — can knowledge say which layer and which markers
would be worth measuring? Answered by `knowledge_agreement`: does an
outcome-blind ordering anticipate what the data ends up keeping, against a
random-ordering null? That is a weaker question than predictive value and is
reported as such.

What is the baseline, and why it decides the whole answer
--------------------------------------------------------
"Which variables are useful" has no answer until "useful on top of what" is
fixed. On MaveDB the baseline is **protein abundance**, which was already by
far the best single predictor, so a selector has to beat the measurement a lab
would make anyway. On HLC it is **age and sex**, the covariates available
without any assay. Pinning the baseline into every panel is what makes the
contrast about *additional* information.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from omics2pk import structural_annotation as SA  # noqa: E402
from omics2pk import supervised_screening_selection as SSS  # noqa: E402
from omics2pk import variant_mechanism_omic_selection as VMS  # noqa: E402

#: Declared before the run. Family size covers the contrasts quoted per cohort.
FAMILY_SIZE = 8

#: The contrasts that decide the commissioned question, per cohort.
PRIMARY = ("S6_knowledge_then_shap", "S1_univariate_outcome")

#: Every selector contrasted against using everything, which is the arm a
#: practitioner gets for free by not selecting at all.
AGAINST_EVERYTHING = ("S1_univariate_outcome", "S2_residual_on_baseline",
                      "S3_elastic_net", "S4_shap_xgboost",
                      "S5_knowledge_only", "S6_knowledge_then_shap",
                      "S7_shap_deduplicated")

#: Full-array protocol, pre-registered in
#: OMICS2PK_HIGH_DIMENSION_PREREGISTRATION_20261111.md.
#: Selectors computable at p = 30 000 without a variance filter.
CHEAP_AT_FULL_WIDTH = ("S0_everything", "S1_univariate_outcome",
                       "S2_residual_on_baseline", "S8_pool_then_univariate")
#: Selectors that need the declared in-fold variance filter to run at all.
NEEDS_PREFILTER = ("S0_everything", "S1_univariate_outcome",
                   "S2_residual_on_baseline", "S3_elastic_net",
                   "S4_shap_xgboost", "S7_shap_deduplicated",
                   "S8_pool_then_univariate")
#: Declared before execution: the number of most-variable candidates kept
#: inside each training fold for the arms that cannot run on the whole array.
PREFILTER_TOP_VARIANCE = 2000

#: Matched random pools, declared before execution. The control that separates
#: biological knowledge from mere dimension reduction: a pool of the same size,
#: matched on level and spread, chosen with no biology at all.
N_MATCHED_POOLS = 40

#: The permuted-outcome control only needs the arms it reports, which is the
#: strongest data-driven selector against no selection at all. Running all
#: seven under permutation would double the cost of the whole run to produce
#: six numbers nobody reads. Declared, not quietly trimmed.
CONTROL_SELECTORS = ("S0_everything", "S4_shap_xgboost")

#: HLC endpoints: three prespecified, one negative and one specificity control.
HLC_PRESPECIFIED = ("CYP3A4T_activity", "CYP2D6_activity", "CYP2C9_activity")
HLC_NEGATIVE_CONTROL = "CYP2E1_activity"
HLC_SPECIFICITY_CONTROL = "CYP3A4M_activity"


class RunError(RuntimeError):
    pass


def _contrasts(result, groups):
    out = {}
    for better in AGAINST_EVERYTHING:
        if better not in result["losses"]:
            continue
        out[f"{better}_vs_S0_everything"] = SSS.contrast(
            result, better, "S0_everything", groups=groups,
            family_size=FAMILY_SIZE)
    if all(name in result["losses"] for name in PRIMARY):
        out["PRIMARY_hybrid_vs_univariate"] = SSS.contrast(
            result, PRIMARY[0], PRIMARY[1], groups=groups,
            family_size=FAMILY_SIZE)
    for better, worse in (("S4_shap_xgboost", "S1_univariate_outcome"),
                          ("S2_residual_on_baseline", "S1_univariate_outcome"),
                          ("S6_knowledge_then_shap", "S4_shap_xgboost"),
                          ("S6_knowledge_then_shap", "S5_knowledge_only"),
                          # isolates deduplication: same ordering, duplicates
                          # removed, nothing else different
                          ("S7_shap_deduplicated", "S4_shap_xgboost"),
                          # the pre-registered primary of the high-dimension
                          # protocol: a mechanistic pool ranked by Spearman
                          # against exact TreeSHAP over the filtered array
                          ("S8_pool_then_univariate", "S4_shap_xgboost"),
                          ("S8_pool_then_univariate", "S1_univariate_outcome")):
        if better in result["losses"] and worse in result["losses"]:
            out[f"{better}_vs_{worse}"] = SSS.contrast(
                result, better, worse, groups=groups, family_size=FAMILY_SIZE)
    return out


def _report(result, groups, *, knowledge_order=None, baseline_k=None):
    block = {
        "n_scored": result["n_scored"],
        "n_candidates": result["n_candidates"],
        "estimator": result["estimator"],
        "k_grid": result["k_grid"],
        "mean_loss": result["mean_loss"],
        "k_chosen_per_fold": result["k_chosen_per_fold"],
        "n_folds_that_abstained": result["n_folds_that_abstained"],
        "abstention_was_available": result["abstention_was_available"],
        "pool_size": result.get("pool_size"),
        "prefilter_top_variance": result.get("prefilter_top_variance"),
        "prefiltered_candidates_per_fold": result.get(
            "prefiltered_candidates_per_fold"),
        "stability": result["stability"],
        "curve": result["curve"],
        "parsimony": {
            name: SSS.parsimony_point(curve, baseline_k=baseline_k)
            for name, curve in result["curve"].items() if len(curve) >= 2},
        "the_curve_is_descriptive": result["the_curve_is_descriptive"],
        "candidate_names": result["candidate_names"],
        "panels_chosen_per_fold": {
            name: result["panels_per_fold"][name]
            for name in ("S4_shap_xgboost", "S6_knowledge_then_shap",
                         "S1_univariate_outcome")
            if name in result["panels_per_fold"]},
        "contrasts": _contrasts(result, groups),
        "baseline_columns": result["baseline_columns"],
    }
    if knowledge_order is not None and "S4_shap_xgboost" in result["losses"]:
        block["can_knowledge_anticipate_the_data_choice"] = {
            str(k): SSS.knowledge_agreement(result, knowledge_order, k=k)
            for k in (3, 8) if k <= result["n_candidates"]}
    return block


# --------------------------------------------------------------------------
# cohort 1: MaveDB, three genes, grouped by protein position and by region
# --------------------------------------------------------------------------
def mavedb_pool(gene, mavedb_dir, pdb_dir):
    """Candidate pool, baseline and a knowledge ordering for one gene."""
    sys.path.insert(0, str(ROOT / "scripts"))
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "variant_bench", ROOT / "scripts" / "run_variant_mechanism_benchmark.py")
    bench = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bench)

    annotation = bench.ANNOTATION[gene]
    directory = pathlib.Path(mavedb_dir)
    activity = bench.load_scores(
        directory / f"mavedb_{annotation['activity_urn'].split(':')[-1]}"
                    "_scores.csv")
    abundance = bench.load_scores(
        directory / f"mavedb_{annotation['abundance_urn'].split(':')[-1]}"
                    "_scores.csv")
    shared = sorted(set(activity) & set(abundance))
    if len(shared) < 200:
        return {"verdict": "NOT_EVALUABLE",
                "reason": f"only {len(shared)} variants in both assays"}
    parsed = [activity[v][0] for v in shared]
    y = np.array([activity[v][1] for v in shared], dtype=float)
    x = np.array([abundance[v][1] for v in shared], dtype=float)
    positions = np.array([p[1] for p in parsed], dtype=int)
    keep = np.array([p[2] != "=" for p in parsed])
    parsed = [p for p, k in zip(parsed, keep) if k]
    y, x, positions = y[keep], x[keep], positions[keep]

    substitution, mechanism = VMS.build_design(parsed, annotation)
    blocks = [("abundance", x.reshape(-1, 1), ("protein_abundance",)),
              ("substitution", substitution,
               tuple(f"substitution_{i}" for i in range(substitution.shape[1]))),
              ("mechanism", mechanism,
               tuple(f"mechanism_{i}" for i in range(mechanism.shape[1])))]

    structural_names = ()
    centre = SA.CATALYTIC_CENTRES.get(gene)
    if pdb_dir and centre:
        path = pathlib.Path(pdb_dir) / f"{centre['pdb']}.pdb"
        if path.is_file():
            text = path.read_text(encoding="utf-8", errors="replace")
            residues, ligands = SA.parse_pdb(text, chain=centre.get("chain"))
            SA.verify_catalytic_geometry(residues, ligands, centre)
            substrate_points = (
                SA.polymer_chain_points(text, centre["substrate_chain"])
                if centre.get("substrate_chain") else None)
            distances = SA.residue_distances(residues, ligands, centre,
                                             substrate_points=substrate_points)
            structural = SA.build_structural_block(positions, distances)
            structural_names = tuple(f"structural_{n}"
                                     for n in SA.STRUCTURAL_NAMES)
            blocks.append(("structural", structural, structural_names))

    design = np.column_stack([block for _, block, _ in blocks])
    names = [n for _, _, these in blocks for n in these]
    # The baseline is protein abundance: the measurement a lab would make
    # anyway, and already the best single predictor on all three genes.
    baseline_columns = (names.index("protein_abundance"),)
    # Knowledge proposes mechanism and structure first, in the order the
    # architecture document ranks them, and never reads the outcome.
    knowledge_first = [i for i, n in enumerate(names)
                       if n.startswith(("mechanism_", "structural_"))]
    knowledge_order = np.array(
        knowledge_first + [i for i in range(len(names))
                           if i not in knowledge_first], dtype=int)
    return {"target": y, "design": design, "groups_position": positions,
            "names": names, "baseline_columns": baseline_columns,
            "knowledge_order": knowledge_order,
            "region_blocks": VMS.region_groups(positions,
                                               VMS.FROZEN["region_blocks"],
                                               annotation["length"]),
            "n_structural": len(structural_names),
            "role": annotation["role"]}


def run_mavedb(genes, mavedb_dir, pdb_dir, *, estimator="ridge"):
    out = {}
    for gene in genes:
        pool = mavedb_pool(gene, mavedb_dir, pdb_dir)
        if pool.get("verdict") == "NOT_EVALUABLE":
            out[gene] = pool
            continue
        print(f"  MaveDB {gene} [{pool['role']}] "
              f"n={pool['target'].size} p={len(pool['names'])}", flush=True)
        entry = {"role": pool["role"], "n": int(pool["target"].size),
                 "candidate_names": pool["names"],
                 "n_structural_columns": pool["n_structural"],
                 "baseline": "protein abundance, the best single predictor"}
        for label, groups in (
                ("grouped_by_position", pool["groups_position"]),
                ("grouped_by_contiguous_region", pool["region_blocks"])):
            result = SSS.screen(
                pool["target"], pool["design"], groups,
                baseline_columns=pool["baseline_columns"],
                knowledge_order=pool["knowledge_order"],
                feature_names=pool["names"], estimator=estimator)
            # The resampling unit is the POSITION in both cases: it is the unit
            # of dependence, not the fold label. Nineteen substitutions at one
            # residue are not nineteen observations.
            entry[label] = _report(result, pool["groups_position"],
                                   knowledge_order=pool["knowledge_order"])
            permuted = SSS.screen(
                pool["target"], pool["design"], groups,
                baseline_columns=pool["baseline_columns"],
                knowledge_order=pool["knowledge_order"],
                feature_names=pool["names"], estimator=estimator,
                permute_target=True, selectors=CONTROL_SELECTORS)
            entry[label]["permuted_outcome_control"] = {
                "mean_loss": permuted["mean_loss"],
                "best_selector_vs_everything": SSS.contrast(
                    permuted, "S4_shap_xgboost", "S0_everything",
                    groups=pool["groups_position"], family_size=FAMILY_SIZE),
                "selectors_run_under_permutation": list(CONTROL_SELECTORS),
                "what_must_happen": (
                    "with the outcome permuted no selector may beat using "
                    "everything by a margin-reaching amount. If one does, the "
                    "harness is reading its own folds."),
            }
        out[gene] = entry
    return out


# --------------------------------------------------------------------------
# cohort 2: HLC, measured in-vitro CYP activity, grouped by donor
# --------------------------------------------------------------------------
def hlc_pools(expression_zip, phenotype_zip, endpoints):
    from omics2pk import hlc_real_cyp_validation as HLC

    loaded = HLC.load_original_hlc(expression_zip, phenotype_zip)
    clinical, rna, rows = loaded["clinical"], loaded["rna"], loaded["gene_rows"]
    panel = np.asarray(sorted({j for g in HLC.GENES for j in rows[g]}), int)
    out, skipped = {}, {}
    for activity in endpoints:
        observed = loaded["outcomes"][activity]
        mask = (np.isfinite(observed) & (observed >= 0)
                & np.isfinite(clinical).all(axis=1)
                & np.isfinite(rna[:, panel]).all(axis=1))
        if mask.sum() < 300:
            skipped[activity] = f"only {int(mask.sum())} donors with complete rows"
            continue
        gene = ("CYP3A4" if activity.startswith("CYP3A4")
                else activity.split("_")[0])
        design = np.column_stack([clinical[mask], rna[mask][:, panel]])
        names = ["age", "sex"] + [f"probe_{loaded['feature_ids'][j]}"
                                  for j in panel]
        # Knowledge proposes the probes of the TARGET enzyme first, then the
        # rest of the mechanistically chosen panel. It never reads activity.
        target_probes = [2 + int(np.flatnonzero(panel == j)[0])
                         for j in rows[gene] if j in set(panel.tolist())]
        knowledge_order = np.array(
            target_probes + [i for i in range(len(names))
                             if i not in target_probes], dtype=int)
        out[activity] = {
            "target": observed[mask], "design": design, "names": names,
            # Age and sex: what is available with no assay at all.
            "baseline_columns": (0, 1),
            "knowledge_order": knowledge_order,
            # One row is one donor, so the donor IS the row here. The cluster
            # bootstrap and the row bootstrap coincide, and that is stated
            # rather than left for a reader to assume.
            "groups": np.arange(int(mask.sum())),
            "n_target_probes": len(target_probes),
            "gene": gene,
        }
    return out, skipped, int(panel.size)


def run_hlc(expression_zip, phenotype_zip, *, estimator="ridge"):
    endpoints = (*HLC_PRESPECIFIED, HLC_NEGATIVE_CONTROL,
                 HLC_SPECIFICITY_CONTROL)
    pools, skipped, panel_size = hlc_pools(expression_zip, phenotype_zip,
                                           endpoints)
    for activity in HLC_PRESPECIFIED:
        if activity not in pools:
            raise RunError(
                f"prespecified endpoint {activity} is absent from the cohort: "
                f"{skipped.get(activity, 'reason not recorded')}. A declared "
                "endpoint that silently vanishes would turn this run into a "
                "different experiment.")
    out = {"n_panel_probes": panel_size, "skipped": skipped,
           "endpoint_roles": {
               **{a: "prespecified" for a in HLC_PRESPECIFIED},
               HLC_NEGATIVE_CONTROL: "negative_control",
               HLC_SPECIFICITY_CONTROL: "substrate_specificity_control"},
           "endpoints": {}}
    for activity, pool in pools.items():
        print(f"  HLC {activity} n={pool['target'].size} "
              f"p={len(pool['names'])}", flush=True)
        result = SSS.screen(
            pool["target"], pool["design"], pool["groups"],
            baseline_columns=pool["baseline_columns"],
            knowledge_order=pool["knowledge_order"],
            feature_names=pool["names"], estimator=estimator)
        entry = _report(result, pool["groups"],
                        knowledge_order=pool["knowledge_order"])
        entry["gene"] = pool["gene"]
        entry["n_target_enzyme_probes"] = pool["n_target_probes"]
        entry["baseline"] = "age and sex, available without any assay"
        entry["one_row_is_one_donor"] = True
        permuted = SSS.screen(
            pool["target"], pool["design"], pool["groups"],
            baseline_columns=pool["baseline_columns"],
            knowledge_order=pool["knowledge_order"],
            feature_names=pool["names"], estimator=estimator,
            permute_target=True, selectors=CONTROL_SELECTORS)
        entry["permuted_outcome_control"] = {
            "mean_loss": permuted["mean_loss"],
            "best_selector_vs_everything": SSS.contrast(
                permuted, "S4_shap_xgboost", "S0_everything",
                groups=pool["groups"], family_size=FAMILY_SIZE)}
        out["endpoints"][activity] = entry
    return out


def hlc_fullarray_pools(expression_zip, phenotype_zip, endpoints):
    """The whole processed array as the candidate pool, with two filters only.

    Both filters are **outcome-blind and declared**: a probe is dropped when
    it is not finite for every retained donor (completeness), and the
    most-variable `PREFILTER_TOP_VARIANCE` of the survivors are kept inside
    each training fold for the arms that cannot otherwise run. Neither reads
    an activity value.
    """
    from omics2pk import hlc_real_cyp_validation as HLC

    loaded = HLC.load_original_hlc(expression_zip, phenotype_zip)
    clinical, rna, rows = loaded["clinical"], loaded["rna"], loaded["gene_rows"]
    mechanistic = np.asarray(sorted({j for g in HLC.GENES for j in rows[g]}),
                             int)
    out, skipped = {}, {}
    for activity in endpoints:
        observed = loaded["outcomes"][activity]
        donors = (np.isfinite(observed) & (observed >= 0)
                  & np.isfinite(clinical).all(axis=1))
        if donors.sum() < 300:
            skipped[activity] = f"only {int(donors.sum())} donors"
            continue
        complete = np.isfinite(rna[donors]).all(axis=0)
        # Keep every mechanistic probe that is complete, and say how many of
        # them the completeness filter removed rather than quietly shrinking
        # the pool the knowledge arm depends on.
        kept = np.flatnonzero(complete)
        mechanistic_kept = np.asarray([j for j in mechanistic if complete[j]],
                                      dtype=int)
        position = {int(c): i for i, c in enumerate(kept)}
        design = np.column_stack([clinical[donors], rna[np.ix_(donors, kept)]])
        names = ["age", "sex"] + [f"probe_{loaded['feature_ids'][j]}"
                                  for j in kept]
        pool = np.asarray([2 + position[int(j)] for j in mechanistic_kept],
                          dtype=int)
        gene = ("CYP3A4" if activity.startswith("CYP3A4")
                else activity.split("_")[0])
        target_probes = np.asarray(
            [2 + position[int(j)] for j in rows[gene] if complete[int(j)]],
            dtype=int)
        # Knowledge order: the target enzyme's probes, then the rest of the
        # mechanistic pool, then everything else. The POOL arm reads only the
        # first `pool_size` of it, which is the mechanistic set exactly.
        rest_of_pool = [int(c) for c in pool if c not in set(target_probes.tolist())]
        head = [int(c) for c in target_probes] + rest_of_pool
        knowledge_order = np.array(
            head + [i for i in range(len(names)) if i not in set(head)],
            dtype=int)
        out[activity] = {
            "target": observed[donors], "design": design, "names": names,
            "baseline_columns": (0, 1),
            "knowledge_order": knowledge_order,
            "groups": np.arange(int(donors.sum())),
            "pool_size": int(pool.size),
            "n_target_probes": int(target_probes.size),
            "n_candidates": len(names),
            "n_probes_dropped_as_incomplete": int((~complete).sum()),
            "n_mechanistic_probes_dropped_as_incomplete": int(
                mechanistic.size - mechanistic_kept.size),
            "gene": gene,
        }
    return out, skipped, int(mechanistic.size)


def matched_pool_null(pool, *, estimator="ridge"):
    """The mechanistic pool against `N_MATCHED_POOLS` random pools of its size.

    Every draw is matched on mean and log variance, excludes the real pool and
    the baseline columns, and is ordered by the same Spearman step. The
    statistic is where the mechanistic pool's out-of-fold loss falls in the
    distribution of matched-random-pool losses.

    Matching is outcome-blind -- mean, variance and correlation involve no
    activity -- so the draws are made once on the retained donors without
    leaking into any fold.
    """
    real_pool = np.asarray(
        pool["knowledge_order"], dtype=int)[:pool["pool_size"]]
    shared = dict(baseline_columns=pool["baseline_columns"],
                  feature_names=pool["names"], estimator=estimator,
                  pool_size=pool["pool_size"],
                  selectors=("S8_pool_then_univariate",))
    real = SSS.screen(pool["target"], pool["design"], pool["groups"],
                      knowledge_order=pool["knowledge_order"], **shared)
    observed = float(real["mean_loss"]["S8_pool_then_univariate"])
    draws, diagnostics = [], []
    for index in range(N_MATCHED_POOLS):
        drawn, report = SSS.matched_random_pool(
            pool["design"], real_pool,
            exclude=pool["baseline_columns"],
            seed=SSS.FROZEN["seed"] + 1000 + index)
        if drawn.size < real_pool.size:
            continue
        order = np.concatenate([
            drawn,
            np.asarray([i for i in range(len(pool["names"]))
                        if i not in set(drawn.tolist())], dtype=int)])
        run = SSS.screen(pool["target"], pool["design"], pool["groups"],
                         knowledge_order=order, **shared)
        draws.append(float(run["mean_loss"]["S8_pool_then_univariate"]))
        if index < 3:
            diagnostics.append(report)
    draws = np.asarray(draws, dtype=float)
    if draws.size < 10:
        return {"verdict": "NOT_EVALUABLE",
                "reason": f"only {int(draws.size)} matched pools were drawable"}
    better = int((draws <= observed).sum())
    return {
        "n_pools": int(draws.size),
        "mechanistic_pool_loss": observed,
        "matched_random_pool_loss_mean": float(draws.mean()),
        "matched_random_pool_loss_best": float(draws.min()),
        "matched_random_pool_loss_5th_percentile": float(
            np.percentile(draws, 5)),
        "n_random_pools_at_least_as_good": better,
        "empirical_p": float((1 + better) / (1 + draws.size)),
        "mechanistic_pool_beats_matched_chance": bool(
            observed < np.percentile(draws, 5)),
        "relative_gain_over_the_mean_random_pool": float(
            1.0 - observed / draws.mean()) if draws.mean() > 0 else float("nan"),
        "match_quality_first_three_draws": diagnostics,
        "what_this_isolates": (
            "the pool size, the probes' expression level and their spread are "
            "held fixed; only the biology used to choose them differs. A "
            "mechanistic pool that does not beat matched chance has reduced "
            "dimension, not applied knowledge."),
    }


def run_hlc_fullarray(expression_zip, phenotype_zip, *, estimator="ridge"):
    """The pre-registered high-dimension protocol. p >> n, for the first time."""
    endpoints = (*HLC_PRESPECIFIED, HLC_NEGATIVE_CONTROL,
                 HLC_SPECIFICITY_CONTROL)
    pools, skipped, mechanistic_size = hlc_fullarray_pools(
        expression_zip, phenotype_zip, endpoints)
    for activity in HLC_PRESPECIFIED:
        if activity not in pools:
            raise RunError(f"prespecified endpoint {activity} is absent: "
                           f"{skipped.get(activity, 'reason not recorded')}")
    out = {"protocol": "OMICS2PK_HIGH_DIMENSION_PREREGISTRATION_20261111.md",
           "mechanistic_panel_size": mechanistic_size,
           "prefilter_top_variance": PREFILTER_TOP_VARIANCE,
           "skipped": skipped, "endpoints": {}}
    for activity, pool in pools.items():
        print(f"  HLC full array {activity} n={pool['target'].size} "
              f"p={pool['n_candidates']} pool={pool['pool_size']}", flush=True)
        shared = dict(
            baseline_columns=pool["baseline_columns"],
            knowledge_order=pool["knowledge_order"],
            feature_names=pool["names"], estimator=estimator,
            pool_size=pool["pool_size"])
        filtered = SSS.screen(pool["target"], pool["design"], pool["groups"],
                              selectors=NEEDS_PREFILTER,
                              prefilter_top_variance=PREFILTER_TOP_VARIANCE,
                              **shared)
        whole = SSS.screen(pool["target"], pool["design"], pool["groups"],
                           selectors=CHEAP_AT_FULL_WIDTH, **shared)
        entry = {
            "gene": pool["gene"],
            "n_candidates": pool["n_candidates"],
            "n_probes_dropped_as_incomplete":
                pool["n_probes_dropped_as_incomplete"],
            "n_mechanistic_probes_dropped_as_incomplete":
                pool["n_mechanistic_probes_dropped_as_incomplete"],
            "pool_size": pool["pool_size"],
            "n_target_enzyme_probes": pool["n_target_probes"],
            "variance_prefiltered": _report(filtered, pool["groups"]),
            "whole_array_no_prefilter": _report(whole, pool["groups"]),
        }
        # THE CONTROL THAT ISOLATES THE KNOWLEDGE. S8 against S1 changes two
        # things at once -- the biology used to choose, and the dimension
        # searched. A random pool of the same size, matched on level and
        # spread, holds the dimension fixed and removes the biology, so what
        # is left between them is what the knowledge contributed.
        entry["matched_random_pools"] = matched_pool_null(
            pool, estimator=estimator)
        entry["non_inferiority"] = {
            "S8_vs_S4": SSS.non_inferiority(
                filtered, "S8_pool_then_univariate", "S4_shap_xgboost",
                groups=pool["groups"], family_size=FAMILY_SIZE),
            "S8_vs_S1_whole_array": SSS.non_inferiority(
                whole, "S8_pool_then_univariate", "S1_univariate_outcome",
                groups=pool["groups"], family_size=FAMILY_SIZE),
            "why_a_second_test": (
                "a non-significant difference is not equivalence. Concluding "
                "that the simple pipeline suffices requires excluding the "
                "possibility that it is worse by an amount that matters, "
                "which is a statement about the LOWER bound against minus "
                "the margin and a different test from the superiority "
                "ladder."),
        }
        # What the variance filter COSTS, paired on identical rows and folds.
        # Reported rather than assumed harmless.
        entry["cost_of_the_variance_prefilter"] = {
            arm: VMS.ablation(whole, filtered, arm, label="whole_array",
                              family_size=FAMILY_SIZE, groups=pool["groups"])
            for arm in ("S1_univariate_outcome", "S2_residual_on_baseline",
                        "S8_pool_then_univariate")
            if arm in whole["losses"] and arm in filtered["losses"]}
        permuted = SSS.screen(pool["target"], pool["design"], pool["groups"],
                              selectors=CONTROL_SELECTORS,
                              prefilter_top_variance=PREFILTER_TOP_VARIANCE,
                              permute_target=True, **shared)
        entry["permuted_outcome_control"] = {
            "mean_loss": permuted["mean_loss"],
            "n_folds_that_abstained": permuted["n_folds_that_abstained"],
            "best_selector_vs_everything": SSS.contrast(
                permuted, "S4_shap_xgboost", "S0_everything",
                groups=pool["groups"], family_size=FAMILY_SIZE)}
        out["endpoints"][activity] = entry
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mavedb-dir", default=None)
    parser.add_argument("--pdb-dir", default=None)
    parser.add_argument("--genes", nargs="*",
                        default=["NUDT15", "CYP2C9", "BRAF"])
    parser.add_argument("--hlc-expression", default=None)
    parser.add_argument("--hlc-phenotype", default=None)
    parser.add_argument("--hlc-fullarray", action="store_true",
                        help="run the pre-registered high-dimension protocol "
                             "on the whole processed array instead of the "
                             "mechanistically targeted panel")
    parser.add_argument("--estimator", default="ridge",
                        choices=("ridge", "xgboost"))
    parser.add_argument("--out", default=str(
        ROOT / "reports" / "data" / "supervised_screening_selection.json"))
    args = parser.parse_args(argv)

    for directory in (args.mavedb_dir, args.pdb_dir):
        if directory is None:
            continue
        resolved = pathlib.Path(directory).resolve()
        if ROOT == resolved or ROOT in resolved.parents:
            raise RunError("source files must stay outside the repository")

    payload = {
        "run": "supervised_screening_selection",
        "generated": "2026-11-11",
        "question": ("Does the methodology choose better which biological "
                     "information to use for predicting enzyme activity?"),
        "frozen": {k: (list(v) if isinstance(v, tuple) else v)
                   for k, v in SSS.FROZEN.items()},
        "selectors": list(SSS.SELECTORS),
        "family_size": FAMILY_SIZE,
        "primary_contrast": f"{PRIMARY[0]} versus {PRIMARY[1]}",
        "estimator": args.estimator,
        "positive_control": SSS.self_test(),
        "what_this_is_not": SSS.what_this_module_does_not_do(),
        "no_individual_level_data_in_this_artefact": True,
    }

    if args.mavedb_dir:
        print("MaveDB cohort ...", flush=True)
        payload["mavedb"] = run_mavedb(args.genes, args.mavedb_dir,
                                       args.pdb_dir, estimator=args.estimator)
    else:
        payload["mavedb"] = {"verdict": "NOT_EVALUABLE",
                             "reason": "--mavedb-dir was not given"}

    if args.hlc_expression and args.hlc_phenotype:
        if args.hlc_fullarray:
            print("HLC FULL ARRAY cohort ...", flush=True)
            payload["hlc_full_array"] = run_hlc_fullarray(
                args.hlc_expression, args.hlc_phenotype,
                estimator=args.estimator)
            payload["hlc"] = {
                "verdict": "NOT_RUN",
                "reason": "--hlc-fullarray selected the high-dimension protocol"}
        else:
            print("HLC cohort ...", flush=True)
            payload["hlc"] = run_hlc(args.hlc_expression, args.hlc_phenotype,
                                     estimator=args.estimator)
    else:
        payload["hlc"] = {
            "verdict": "NOT_EVALUABLE",
            "reason": ("--hlc-expression and --hlc-phenotype were not both "
                       "given. The HLC source is access-restricted and its "
                       "per-donor rows never enter this repository.")}

    payload["epa"] = {
        "verdict": "NOT_EVALUABLE",
        "reason": ("no EPA eight-CYP archive is present on this machine. "
                   "`epa_eight_cyp_sensitivity` and `epa_liver_cyp_validation` "
                   "exist as readers, but a reader without data is not a "
                   "cohort, and inventing one would break the brief."),
    }

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(payload, indent=2, ensure_ascii=False, default=float)
    out.write_text(text + "\n", encoding="utf-8")
    print(f"\nartefact: {out}  sha256 "
          f"{hashlib.sha256(text.encode('utf-8')).hexdigest()[:16]}")

    for cohort in ("mavedb", "hlc", "hlc_full_array"):
        if cohort not in payload:
            continue
        block = payload[cohort]
        if block.get("verdict") == "NOT_EVALUABLE":
            print(f"\n{cohort.upper()}: NOT_EVALUABLE - {block['reason']}")
            continue
        print(f"\n===== {cohort.upper()} =====")
        items = (block.get("endpoints") or block).items()
        for key, entry in items:
            if not isinstance(entry, dict) or entry.get("verdict"):
                continue
            for label in ("grouped_by_position",
                          "grouped_by_contiguous_region", None):
                cell = entry if label is None else entry.get(label)
                if not isinstance(cell, dict) or "contrasts" not in cell:
                    continue
                head = key if label is None else f"{key} / {label}"
                print(f"  --- {head}")
                print("      " + "  ".join(
                    f"{n.split('_')[0]}={v:.4f}"
                    for n, v in cell["mean_loss"].items()))
                print("      k chosen: " + "  ".join(
                    f"{n.split('_')[0]}={v}"
                    for n, v in cell["k_chosen_per_fold"].items()))
                for name, value in cell["contrasts"].items():
                    print(f"      {name:<46} "
                          f"{value['relative_gain']:+7.2%} "
                          f"[{value['interval'][0]:+.4f},"
                          f"{value['interval'][1]:+.4f}] {value['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
