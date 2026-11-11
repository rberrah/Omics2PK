#!/usr/bin/env python
"""Compute the corrected fusion-headroom diagnostics on ARIANA, Lee and HLC.

No fusion method is implemented and no gate is fitted. This runner collects
out-of-fold predictions from arms that **already exist** and reports four
ceilings, two recoverable quantities and a context probe, under the corrected
rules in `src/omics2pk/fusion_headroom_diagnostics.py`:

* the constant-weight ceiling bounds **stacking only**, not a contextual gate;
* the per-individual **convex** ceiling bounds every method in the family;
* an interval containing zero is a **non-detection**, never a refutation;
* the weights are fitted on one partition and scored on another, because
  fitting them on the same out-of-fold predictions that measure them leaks
  even when the base predictions are honestly out-of-fold.

Usage:
    python -I scripts/run_fusion_headroom_diagnostics.py \
        --ariana <dir>/Solanidine_metoprolol_ratios_PGx.xlsx \
        --lee-per-donor <outside-repo>/lee_jackson_cyp3a_panel_per_donor.csv \
        --out reports/data/fusion_headroom_diagnostics.json
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.stdout.reconfigure(encoding="utf-8")

from omics2pk import fusion_headroom_diagnostics as FHD  # noqa: E402
from omics2pk import plasma_biomarker_validation as PBV  # noqa: E402

PROTOCOL = ROOT / "docs" / "OMICS2PK_FUSING_REFUTATION_PROTOCOL_20261110.md"

#: The declared family for the Bonferroni level: one headroom gap per cohort,
#: times the two ceilings that get a verdict.
DECLARED_FAMILY_SIZE = 6

#: Ten HLC endpoints times two context sets. Declared before the run so
#: that a probe clearing only at a nominal level is labelled as such.
PROBE_FAMILY_SIZE = 20


class RunError(RuntimeError):
    pass


def _load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


# --------------------------------------------------------------------------
# the diagnostic block, identical for every cohort
# --------------------------------------------------------------------------
def diagnose(predictions, truth, context, *, arm_names, context_names,
             family_size=DECLARED_FAMILY_SIZE):
    matrix = np.asarray(predictions, dtype=float)
    target = np.asarray(truth, dtype=float)

    single_oracle = FHD.best_single_oracle(matrix, target)
    single_honest = FHD.best_single_honest(matrix, target)
    constant = FHD.constant_weight_oracle(matrix, target)
    selection = FHD.individual_selection_oracle(matrix, target)
    convex = FHD.individual_convex_oracle(matrix, target)
    stacked = FHD.nested_constant_stacking(matrix, target)

    reference = single_honest["per_individual_losses"]
    per_arm = np.abs(matrix - target[:, None])
    constant_weights = np.array([constant["weights"][f"arm_{i}"]
                                 for i in range(matrix.shape[1])])

    gaps = {
        # Bounds stacking only. A zero gap here refutes STACKING, not gates.
        "reference_minus_constant_ceiling": FHD.paired_gap(
            reference, np.abs(matrix @ constant_weights - target),
            family_size=family_size),
        # Bounds per-individual selection.
        "reference_minus_selection_ceiling": FHD.paired_gap(
            reference, per_arm.min(axis=1), family_size=family_size),
        # Bounds EVERYTHING, gates included. This is the family's ceiling.
        "reference_minus_convex_ceiling": FHD.paired_gap(
            reference, convex["per_individual_losses"],
            family_size=family_size),
        # What stacking actually recovers, honestly.
        "reference_minus_nested_stacking": FHD.paired_gap(
            reference, stacked["per_individual_losses"],
            family_size=family_size),
    }

    verdicts = {
        "stacking_family": FHD.verdict(
            gaps["reference_minus_constant_ceiling"],
            single_honest["mean_loss"]),
        "whole_fusion_family": FHD.verdict(
            gaps["reference_minus_convex_ceiling"],
            single_honest["mean_loss"]),
    }

    # Two context sets, because a gate may condition on anything observable
    # WITHOUT the outcome, and the arms' own disagreement is observable.
    # Testing only the clinical covariates would let a refutation of
    # contextual fusion rest on having looked in one place.
    spread = matrix.max(axis=1) - matrix.min(axis=1)
    endogenous = np.column_stack(
        [matrix, spread, matrix.mean(axis=1), matrix.std(axis=1)])
    endogenous_names = ([f"prediction::{name}" for name in arm_names]
                        + ["prediction_spread", "prediction_mean",
                           "prediction_sd"])

    probes = {}
    for label, block, names in (
            ("clinical_covariates", context, list(context_names)),
            ("arm_predictions_and_their_disagreement", endogenous,
             endogenous_names)):
        probe = FHD.context_informativeness(
            matrix, target, block, family_size=PROBE_FAMILY_SIZE)
        probe["context_variables"] = names
        probe["compared_arm_names"] = [arm_names[i]
                                       for i in probe["compared_arms"]]
        probes[label] = probe
    probe = probes["clinical_covariates"]

    errors = per_arm - per_arm.mean(axis=0)
    norms = np.sqrt((errors ** 2).sum(axis=0))
    safe = np.where(norms > 0.0, norms, 1.0)
    correlation = (errors / safe).T @ (errors / safe)
    off_diagonal = correlation[~np.eye(matrix.shape[1], dtype=bool)]

    return {
        "n_individuals": int(target.size),
        "arms": list(arm_names),
        "n_arms": int(matrix.shape[1]),
        "ceilings": {
            "best_single_oracle": {k: v for k, v in single_oracle.items()},
            "best_single_honest_training_choice": {
                "mean_loss": single_honest["mean_loss"],
                "arm_chosen_per_partition": [
                    arm_names[i] for i in
                    single_honest["arm_chosen_per_partition"]],
            },
            "constant_weight_oracle": {
                k: v for k, v in constant.items() if k != "weights"},
            "constant_weight_oracle_weights": {
                arm_names[i]: constant["weights"][f"arm_{i}"]
                for i in range(matrix.shape[1])},
            "individual_selection_oracle": {
                "mean_loss": selection["mean_loss"],
                "bounds": selection["bounds"],
                "best_arm_shares": {
                    arm_names[i]: selection["best_arm_shares"][f"arm_{i}"]
                    for i in range(matrix.shape[1])},
            },
            "individual_convex_oracle": {
                k: v for k, v in convex.items()
                if k != "per_individual_losses"},
        },
        "recoverable": {
            "nested_constant_stacking": {
                k: v for k, v in stacked.items()
                if k != "per_individual_losses"},
        },
        "gaps": gaps,
        "verdicts": verdicts,
        "context_probes": probes,
        "error_structure": {
            "median_pairwise_error_correlation": float(
                np.median(off_diagonal)),
            "max_pairwise_error_correlation": float(off_diagonal.max()),
            "how_to_read": (
                "combination gain is driven by decorrelation. A median above "
                "0.95 means the arms are one arm and no weighting can help."),
        },
    }


# --------------------------------------------------------------------------
# ARIANA
# --------------------------------------------------------------------------
def ariana_block(source):
    h1 = _load_module(ROOT / "scripts"
                      / "run_h1_ariana_incremental_omic_value.py", "h1_ariana")
    h1.assert_no_prohibited_predictor()
    frame, _ = h1.load_ariana(source)
    target_column = h1.COL_TROUGH
    observed = frame[target_column].to_numpy(dtype=float)
    if np.any(observed <= 0):
        raise RunError("a non-positive concentration cannot be log-transformed")
    keep = np.isfinite(observed)
    y_log = np.log(observed[keep])
    folds = h1.make_folds(int(keep.sum()), h1.OUTER_FOLDS, h1.OUTER_REPEATS,
                          h1.SEED)
    factories = h1.model_factories()
    factory, grid = factories[h1.PRIMARY_MODEL]

    # `out_of_fold_predictions` returns one complete out-of-fold assignment
    # PER REPEAT, and this project's declared loss averages the per-patient
    # ERROR across repeats rather than the predictions. Averaging predictions
    # would invent a prediction no arm made, so the diagnostic is computed
    # separately on each repeat and both are reported. Disagreement between
    # repeats is itself information about stability at n = 47.
    arm_names, per_arm_repeats = [], []
    for name, features in h1.ARMS.items():
        design = np.column_stack(
            [frame[f].to_numpy(dtype=float)[keep] for f in features])
        if not np.isfinite(design).all():
            continue
        predicted_log = np.asarray(
            h1.out_of_fold_predictions(design, y_log, folds, factory, grid),
            dtype=float)
        arm_names.append(name)
        per_arm_repeats.append(np.exp(predicted_log))

    context_names = h1.CLINICAL
    context = np.column_stack(
        [frame[f].to_numpy(dtype=float)[keep] for f in context_names])

    n_repeats = per_arm_repeats[0].shape[0]
    block = {
        "per_repeat": {
            f"repeat_{repeat}": diagnose(
                np.column_stack([arm[repeat] for arm in per_arm_repeats]),
                observed[keep], context, arm_names=arm_names,
                context_names=context_names)
            for repeat in range(n_repeats)},
        "n_repeats": int(n_repeats),
        "why_per_repeat": (
            "each repeat is one complete out-of-fold assignment. This "
            "project's loss averages per-patient ERRORS across repeats; "
            "averaging PREDICTIONS would create a prediction no arm made, so "
            "the ceilings are computed per repeat and every repeat is "
            "reported."),
        "n_individuals": int(keep.sum()),
        "arms": list(arm_names),
        "n_arms": len(arm_names),
    }
    block["cohort"] = "ARIANA n=47, observed pre-dose metoprolol concentration"
    block["loss_units"] = "ng/mL, absolute error"
    block["endpoint_caveat"] = (
        "an observed pre-dose concentration. NOT an AUC, NOT a clearance, NOT "
        "a steady-state exposure metric.")
    block["status_ceiling"] = "exploratory; this cohort has been analysed "\
        "repeatedly by this project"
    block["source_sha256"] = hashlib.sha256(
        pathlib.Path(source).read_bytes()).hexdigest()
    return block


# --------------------------------------------------------------------------
# Lee / Jackson
# --------------------------------------------------------------------------
def lee_block(per_donor):
    path = pathlib.Path(per_donor).resolve()
    if ROOT == path.parent or ROOT in path.parents:
        raise RunError("the per-donor table must stay outside this repository")
    rows, order = {}, []
    with open(path, newline="", encoding="utf-8") as handle:
        for record in csv.DictReader(handle):
            donor = record.pop("donor")
            parsed = {}
            for key, value in record.items():
                try:
                    parsed[key] = float(value)
                except (TypeError, ValueError):
                    parsed[key] = value
            parsed["cyp3a_total_abundance"] = (parsed["cyp3a4_abundance"]
                                               + parsed["cyp3a5_abundance"])
            parsed["is_male"] = 1.0 if str(parsed.get("sex", "")).upper() == "M" \
                else 0.0
            rows[donor] = parsed
            order.append(donor)

    endpoint = PBV.FROZEN["primary_endpoint"]
    arms = {**PBV.ARMS, **PBV.SENSITIVITY_ARMS}
    arm_names, columns = [], []
    target = PBV.log_target(rows, order, endpoint)
    for name, features in arms.items():
        predicted, _, _ = PBV.leave_one_donor_out(
            PBV.design(rows, order, features), target)
        arm_names.append(name)
        columns.append(predicted)

    context_names = ["age_years", "is_male", "bmi"]
    context = np.array([[rows[d][c] for c in context_names] for d in order],
                       dtype=float)
    block = diagnose(np.column_stack(columns), target, context,
                     arm_names=arm_names, context_names=context_names)
    block["cohort"] = ("Lee / Jackson n=15 hepatocyte donors, midazolam "
                       "1'-hydroxylation")
    block["loss_units"] = "log10 of pmol per million cells per minute"
    block["endpoint_caveat"] = (
        "an activity rate at a single substrate concentration. NOT an "
        "intrinsic clearance and NOT a pharmacokinetic quantity.")
    block["status_ceiling"] = "exploratory; n = 15"
    block["per_donor_table_sha256"] = hashlib.sha256(
        path.read_bytes()).hexdigest()
    return block


# --------------------------------------------------------------------------
# Human Liver Cohort — the only cohort where a gate could be estimated
# --------------------------------------------------------------------------
def hlc_block(expression_zip, phenotype_zip, endpoints=None):
    """Four arms on identical folds, per activity endpoint, n ~ 400.

    The designs, the folds, the model and the seed are taken unchanged from
    `hlc_real_cyp_validation.evaluate_fold_matched_cyp`. The only thing added
    is retaining the PREDICTIONS, which that function discards because it
    needs only the errors. No frozen artefact is touched and no model changes.
    """
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import Ridge
    from sklearn.model_selection import KFold
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    from omics2pk import hlc_real_cyp_validation as HLC

    loaded = HLC.load_original_hlc(expression_zip, phenotype_zip)
    clinical, rna, rows = loaded["clinical"], loaded["rna"], loaded["gene_rows"]
    panel = np.asarray(sorted({j for g in HLC.GENES for j in rows[g]}), int)
    wanted = list(endpoints) if endpoints else list(HLC.ACTIVITIES)

    out = {}
    for activity in wanted:
        observed = loaded["outcomes"][activity]
        mask = (np.isfinite(observed) & (observed >= 0)
                & np.isfinite(clinical).all(axis=1))
        if mask.sum() < 300:
            out[activity] = {"verdict": "NOT_EVALUABLE",
                             "reason": f"only {int(mask.sum())} donors"}
            continue
        gene = ("CYP3A4" if activity.startswith("CYP3A4")
                else activity.split("_")[0])
        target = observed[mask]
        designs = {
            "age_sex": clinical[mask],
            "target_CYP_RNA_only": rna[mask][:, rows[gene]],
            "age_sex_plus_target_CYP_RNA": np.column_stack(
                (clinical[mask], rna[mask][:, rows[gene]])),
            "age_sex_plus_targeted_panel": np.column_stack(
                (clinical[mask], rna[mask][:, panel])),
        }
        folds = list(KFold(5, shuffle=True,
                           random_state=HLC.SEED).split(target))
        arm_names, columns = [], []
        for name, design in designs.items():
            predicted = np.empty(target.size)
            for train, test in folds:
                model = make_pipeline(SimpleImputer(strategy="median"),
                                      StandardScaler(), Ridge(alpha=1.0))
                model.fit(design[train], np.log1p(target[train]))
                predicted[test] = np.maximum(
                    0.0, np.expm1(np.clip(model.predict(design[test]),
                                          -20.0, 20.0)))
            arm_names.append(name)
            columns.append(predicted)
        block = diagnose(np.column_stack(columns), target, clinical[mask],
                         arm_names=arm_names,
                         context_names=["age_years", "is_male"])
        block["cohort"] = f"HLC n={int(mask.sum())}, {activity}"
        block["loss_units"] = "original source activity units, absolute error"
        block["endpoint_caveat"] = (
            "measured in-vitro microsomal enzyme activity. NOT a clinical "
            "pharmacokinetic quantity.")
        block["context_caveat"] = (
            "the only clinical context this cohort publishes is age and sex. "
            "A null context probe here is a statement about those two "
            "variables, not about context in general.")
        out[activity] = block
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ariana")
    parser.add_argument("--lee-per-donor")
    parser.add_argument("--hlc-expression")
    parser.add_argument("--hlc-phenotype")
    parser.add_argument("--hlc-endpoints", nargs="*", default=None)
    parser.add_argument("--out",
                        default=str(ROOT / "reports" / "data"
                                    / "fusion_headroom_diagnostics.json"))
    args = parser.parse_args(argv)

    if not PROTOCOL.is_file():
        raise RunError(f"the protocol is missing: {PROTOCOL}")

    cohorts, skipped = {}, {}
    if args.ariana:
        cohorts["ariana"] = ariana_block(args.ariana)
    else:
        skipped["ariana"] = "NOT_EVALUABLE: --ariana not supplied"
    if args.lee_per_donor:
        cohorts["lee_jackson"] = lee_block(args.lee_per_donor)
    else:
        skipped["lee_jackson"] = "NOT_EVALUABLE: --lee-per-donor not supplied"
    if args.hlc_expression and args.hlc_phenotype:
        for activity, block in hlc_block(args.hlc_expression,
                                         args.hlc_phenotype,
                                         args.hlc_endpoints).items():
            if block.get("verdict") == "NOT_EVALUABLE":
                skipped[f"hlc::{activity}"] = block["reason"]
            else:
                cohorts[f"hlc::{activity}"] = block
        skipped["hlc_h2a_arms"] = (
            "The eight H2a arms are stored as LOSSES, not per-donor "
            "predictions, so their convex ceiling cannot be recomputed from "
            "the frozen artefacts. The four arms used here come from "
            "hlc_real_cyp_validation with its own designs, folds, model and "
            "seed unchanged. Extending the diagnostic to the H2a arms means "
            "re-running them to retain predictions, which is a separate job "
            "and still not a Fusing implementation.")
    else:
        skipped["hlc"] = ("NOT_EVALUABLE: --hlc-expression and "
                          "--hlc-phenotype not supplied")

    payload = {
        "run": "fusion_headroom_diagnostics",
        "generated": "2026-11-11",
        "protocol": PROTOCOL.name,
        "protocol_sha256": hashlib.sha256(PROTOCOL.read_bytes()).hexdigest(),
        "frozen": {k: (list(v) if isinstance(v, tuple) else v)
                   for k, v in FHD.FROZEN.items()},
        "declared_family_size": DECLARED_FAMILY_SIZE,
        "corrects": {
            "what_was_wrong": (
                "the protocol's Gate 1 used a single CONSTANT-weight ceiling "
                "and proposed to stop the programme when it showed no "
                "detectable headroom. A constant weight vector cannot express "
                "a gate, so that ceiling does not bound contextual fusion, "
                "and the rule would have discarded the hypothesis it was "
                "meant to test."),
            "what_replaces_it": (
                "four ceilings, of which the per-individual CONVEX one bounds "
                "the whole family; a verdict rule in which an interval "
                "containing zero is a non-detection and refutation requires "
                "the upper bound to fall below a declared minimum interesting "
                "gain; weights fitted on one partition and scored on another; "
                "and a context probe that can refute contextual fusion "
                "independently of any ceiling."),
            "proof_that_the_correction_matters": FHD.self_test(),
        },
        "cohorts": cohorts,
        "not_evaluated": skipped,
        "what_this_run_does_not_do": FHD.what_this_module_does_not_do(),
        "no_individual_level_data_in_this_artefact": True,
    }

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
                   encoding="utf-8")

    flat = {}
    for name, block in cohorts.items():
        if "per_repeat" in block:
            for label, inner in block["per_repeat"].items():
                merged = dict(inner)
                merged["cohort"] = f"{block['cohort']} [{label}]"
                flat[f"{name}::{label}"] = merged
        else:
            flat[name] = block
    for name, block in flat.items():
        print(f"\n===== {name}: {block['cohort']} "
              f"({block['n_arms']} arms, n={block['n_individuals']}) =====")
        ceilings = block["ceilings"]
        print(f"  best single, honest training choice   "
              f"{ceilings['best_single_honest_training_choice']['mean_loss']:.4f}")
        print(f"  best single, knowing the answers      "
              f"{ceilings['best_single_oracle']['mean_loss']:.4f}")
        print(f"  constant-weight ceiling (stacking)    "
              f"{ceilings['constant_weight_oracle']['mean_loss']:.4f}")
        print(f"  per-individual SELECTION ceiling      "
              f"{ceilings['individual_selection_oracle']['mean_loss']:.4f}")
        print(f"  per-individual CONVEX ceiling         "
              f"{ceilings['individual_convex_oracle']['mean_loss']:.4f}"
              f"   (bracketing rate "
              f"{ceilings['individual_convex_oracle']['bracketing_rate']:.3f})")
        print(f"  nested stacking, actually recovered   "
              f"{block['recoverable']['nested_constant_stacking']['mean_loss']:.4f}")
        for key, gap in block["gaps"].items():
            print(f"    {key:<42} {gap['gap']:+.4f} "
                  f"[{gap['interval'][0]:+.4f}, {gap['interval'][1]:+.4f}]")
        for key, cell in block["verdicts"].items():
            print(f"    VERDICT {key:<22} {cell['verdict']}")
        for label, probe in block["context_probes"].items():
            majority = max(probe["share_where_the_runner_up_wins"],
                           1.0 - probe["share_where_the_runner_up_wins"])
            print(f"  probe [{label}]: rho "
                  f"{probe['spearman_out_of_sample']:+.4f}, sign agreement "
                  f"{probe['sign_agreement']:.3f} vs majority "
                  f"{majority:.3f}")
        print(f"  median pairwise error correlation     "
              f"{block['error_structure']['median_pairwise_error_correlation']:.4f}")
    for name, reason in skipped.items():
        print(f"\n{name}: {reason}")
    print(f"\nartefact: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
