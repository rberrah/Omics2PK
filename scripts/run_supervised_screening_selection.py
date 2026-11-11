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
                      "S5_knowledge_only", "S6_knowledge_then_shap")

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
                          ("S6_knowledge_then_shap", "S5_knowledge_only")):
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


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mavedb-dir", default=None)
    parser.add_argument("--pdb-dir", default=None)
    parser.add_argument("--genes", nargs="*",
                        default=["NUDT15", "CYP2C9", "BRAF"])
    parser.add_argument("--hlc-expression", default=None)
    parser.add_argument("--hlc-phenotype", default=None)
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

    for cohort in ("mavedb", "hlc"):
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
