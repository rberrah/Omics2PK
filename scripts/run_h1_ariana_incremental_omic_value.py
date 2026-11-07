"""H1 — incremental value of omics for metoprolol trough prediction in ARIANA.

Executes the programme's single PRIMARY CONTRAST, frozen in
OMICS2PK_FIVE_OBJECTIVES_VALIDATION_PROTOCOL v1.0.0 §4:

    ARIANA, Ridge M3 versus Ridge M1, MAE, identical folds.

Ridge rather than the mechanistic model on purpose: Ridge estimates its
coefficients from the cohort and cannot borrow strength from literature-fixed
constants, so it is the harder test of whether the omic variables carry
incremental information.

Prespecification, stated exactly
--------------------------------
Frozen in protocol v1.0.0 (commit 7632e9e): the arm definitions, the primary
contrast and its loss, the loss hierarchy with its reversal rule, the repeated
subject-disjoint nested cross-validation scheme, the mandatory influence
reporting, the prohibition on using the OH-metoprolol/metoprolol ratio as a
predictor, and the exclusion of MAP-BE.

Operationalised at analysis time on 2026-11-07, NOT frozen in advance: fitting
on the log of the trough with back-transformation, the Ridge penalty grid, the
Elastic Net l1 ratio grid, the feature encodings, and the bootstrap resample
count. These are declared here before any result is reported.

Comparators: regularized Ridge, Elastic Net and XGBoost, as the protocol names,
with complexity capped for n=47. No substitution is in force.

What is forbidden here
----------------------
The ``lnMR OH-metoprolol/metoprolol`` column shares its denominator with the
target concentration. It is NEVER a predictor, in any arm. MAP-BE is excluded
by construction: one observation per patient makes conditioning on the target
and then predicting it circular. The endpoint is an observed pre-dose
concentration and is never called AUC or clearance.

Only aggregate statistics are printed; no individual prediction is emitted.

Usage
-----
    python scripts/run_h1_ariana_incremental_omic_value.py --source PATH
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

EXPECTED_SHA256 = "7397b13e43cd546b93f83a20b2dc84b6d49e97633173d56711c39973beedfe37"

# --- columns -------------------------------------------------------------
COL_TROUGH = "c(Metoprolol)[ng/ml]"
COL_TROUGH_WEIGHTED = "c(metoprolol)weighted [ng/ml]"
COL_PROHIBITED = "lnMR OH-metoprolol/metoprolol"
COL_OH = "lnMR OH-solanidine/solanidine"
COL_SSDA = "lnMR SSDA/solanidine"
COL_AS = "AS"
COL_INHIBITOR = "CYP2D6 inhibitor intake"
COL_DOSE = "Daily Dose Metoprolol"

CLINICAL = ["Age (years)", "male sex", "weight", "BMI", COL_DOSE]
PGX = [COL_AS, COL_INHIBITOR]
OMIC = [COL_OH, COL_SSDA]

#: Arm definitions, frozen in protocol v1.0.0 §3. The mechanistic arms M4/M5
#: and the hybrid HY arm are not computed here: they are the designated key
#: secondary and are run separately, so this module covers the arms the PRIMARY
#: contrast needs plus their statistical neighbours.
ARMS = {
    "M0_clinical": CLINICAL,
    "M1_clinical_pgx": CLINICAL + PGX,
    "M2_clinical_omic": CLINICAL + OMIC,
    "M3_clinical_pgx_omic": CLINICAL + PGX + OMIC,
    # --- PGx decomposition, added 2026-11-07 as a secondary analysis -------
    # The M1/M3 PGx block bundles the genotype-derived activity score with the
    # comedication flag. Bundled, an effect cannot be attributed to either, so
    # these arms separate them. Secondary and exploratory by construction.
    "M1a_clinical_genotype": CLINICAL + [COL_AS],
    "M1b_clinical_comedication": CLINICAL + [COL_INHIBITOR],
    "M3a_clinical_genotype_omic": CLINICAL + [COL_AS] + OMIC,
    "M3b_clinical_comedication_omic": CLINICAL + [COL_INHIBITOR] + OMIC,
}
PRIMARY_CONTRAST = ("M3_clinical_pgx_omic", "M1_clinical_pgx")
PRIMARY_MODEL = "ridge"
PRIMARY_LOSS = "mae"

# --- validation machinery, frozen -----------------------------------------
OUTER_FOLDS = 5
OUTER_REPEATS = 2
INNER_FOLDS = 3
SEED = 20261107
N_BOOTSTRAP = 2000

#: Exploratory success threshold, frozen in protocol v1.0.0 §4.
RELATIVE_IMPROVEMENT_THRESHOLD = 0.10
#: Trimming fraction for the robustness co-criterion, frozen in §2.1.
TRIM_FRACTION = 0.10


class AnalysisError(RuntimeError):
    pass


def sha256_of(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_source(path):
    actual = sha256_of(path)
    if actual != EXPECTED_SHA256:
        raise AnalysisError(
            f"checksum mismatch for {path}\n"
            f"  expected {EXPECTED_SHA256}\n  actual   {actual}"
        )
    return actual


def load_ariana(path):
    import pandas as pd

    checksum = verify_source(path)
    frame = pd.read_excel(path)
    frame.columns = [str(c).strip() for c in frame.columns]
    required = set(CLINICAL + PGX + OMIC + [COL_TROUGH, COL_TROUGH_WEIGHTED])
    missing = sorted(required - set(frame.columns))
    if missing:
        raise AnalysisError(f"source workbook is missing columns: {missing}")
    if len(frame) != 47:
        raise AnalysisError(f"expected 47 patients, found {len(frame)}")
    return frame, checksum


def assert_no_prohibited_predictor():
    """The target-contaminated ratio must not appear in any arm."""
    for arm, features in ARMS.items():
        if COL_PROHIBITED in features:
            raise AnalysisError(
                f"arm {arm!r} contains the prohibited predictor {COL_PROHIBITED!r}"
            )


# --- models ---------------------------------------------------------------
def model_factories():
    """Return {name: (estimator factory, inner-CV parameter grid)}."""
    from sklearn.linear_model import ElasticNet, Ridge
    from xgboost import XGBRegressor
    import numpy as np

    return {
        "ridge": (
            lambda **kw: Ridge(**kw),
            {"alpha": list(np.logspace(-3, 3, 25))},
        ),
        "elastic_net": (
            lambda **kw: ElasticNet(max_iter=50000, **kw),
            {
                "alpha": list(np.logspace(-3, 1, 17)),
                "l1_ratio": [0.1, 0.5, 0.9, 1.0],
            },
        ),
        # Complexity deliberately capped: n=47 cannot support a deep ensemble.
        "xgboost": (
            lambda **kw: XGBRegressor(
                random_state=SEED,
                n_jobs=1,
                tree_method="exact",
                objective="reg:squarederror",
                **kw,
            ),
            {
                "n_estimators": [50, 150],
                "max_depth": [1, 2],
                "learning_rate": [0.05, 0.1],
                "reg_lambda": [1.0, 10.0],
            },
        ),
    }


def _grid(parameter_grid):
    from itertools import product

    keys = sorted(parameter_grid)
    for values in product(*(parameter_grid[key] for key in keys)):
        yield dict(zip(keys, values))


def make_folds(n, n_splits, n_repeats, seed):
    """Subject-disjoint outer folds, identical for every arm and model."""
    from sklearn.model_selection import RepeatedKFold

    splitter = RepeatedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=seed)
    return [(train, test) for train, test in splitter.split(range(n))]


def out_of_fold_predictions(X, y_log, folds, factory, parameter_grid):
    """Nested CV. Returns a (n_repeats, n) array of log-scale predictions.

    Preprocessing and hyperparameter selection happen inside the training part
    of each outer fold only. No information crosses from test to train.
    """
    import numpy as np
    from sklearn.model_selection import KFold
    from sklearn.preprocessing import StandardScaler

    n = X.shape[0]
    folds_per_repeat = OUTER_FOLDS
    n_repeats = len(folds) // folds_per_repeat
    predictions = np.full((n_repeats, n), np.nan)

    for fold_index, (train_idx, test_idx) in enumerate(folds):
        repeat = fold_index // folds_per_repeat
        X_train, X_test = X[train_idx], X[test_idx]
        y_train = y_log[train_idx]

        # --- inner CV for hyperparameters, training fold only -------------
        inner = KFold(n_splits=INNER_FOLDS, shuffle=True, random_state=SEED + repeat)
        best_params, best_score = None, np.inf
        for params in _grid(parameter_grid):
            errors = []
            for inner_train, inner_val in inner.split(X_train):
                scaler = StandardScaler().fit(X_train[inner_train])
                model = factory(**params)
                model.fit(scaler.transform(X_train[inner_train]), y_train[inner_train])
                predicted = model.predict(scaler.transform(X_train[inner_val]))
                errors.append(np.mean((predicted - y_train[inner_val]) ** 2))
            score = float(np.mean(errors))
            if score < best_score:
                best_params, best_score = params, score

        # --- refit on the whole training fold -----------------------------
        scaler = StandardScaler().fit(X_train)
        model = factory(**best_params)
        model.fit(scaler.transform(X_train), y_train)
        predictions[repeat, test_idx] = model.predict(scaler.transform(X_test))

    if np.isnan(predictions).any():
        raise AnalysisError("some patients received no out-of-fold prediction")
    return predictions


# --- losses ---------------------------------------------------------------
def per_patient_absolute_error(y_true, predictions_log):
    """Mean across repeats of the per-patient absolute error, original scale."""
    import numpy as np

    predicted = np.exp(predictions_log)
    return np.mean(np.abs(predicted - y_true[None, :]), axis=0)


def per_patient_log_absolute_error(y_log, predictions_log):
    import numpy as np

    return np.mean(np.abs(predictions_log - y_log[None, :]), axis=0)


def trimmed_mean(values, fraction=TRIM_FRACTION):
    """Mean after removing the largest `fraction` of values.

    One-sided: the robustness question is whether a few extreme patients carry
    the result, so the largest errors are what gets trimmed.
    """
    import numpy as np

    array = np.sort(np.asarray(values, dtype=float))
    keep = max(1, int(round(array.size * (1.0 - fraction))))
    return float(array[:keep].mean())


def loss_summary(y_true, y_log, predictions_log):
    import numpy as np

    predicted = np.exp(predictions_log)
    mean_predicted = predicted.mean(axis=0)
    absolute = per_patient_absolute_error(y_true, predictions_log)
    log_absolute = per_patient_log_absolute_error(y_log, predictions_log)
    within_20 = np.mean(np.abs(mean_predicted - y_true) / y_true <= 0.20)
    return {
        "mae": float(absolute.mean()),
        "trimmed_mae_10pct": trimmed_mean(absolute),
        "log_mae": float(log_absolute.mean()),
        "rmse": float(np.sqrt(np.mean((mean_predicted - y_true) ** 2))),
        "bias": float(np.mean(mean_predicted - y_true)),
        "p20_fraction_within_20pct": float(within_20),
        "n_patients": int(y_true.size),
    }


def paired_contrast(errors_a, errors_b, n_boot=N_BOOTSTRAP, seed=SEED):
    """Paired difference a minus b, bootstrapped over PATIENTS.

    Patients are the independent unit. Fold repetition is never treated as
    additional independent subjects, per protocol §9.
    """
    import numpy as np

    errors_a = np.asarray(errors_a, dtype=float)
    errors_b = np.asarray(errors_b, dtype=float)
    difference = errors_a - errors_b
    point = float(difference.mean())

    rng = np.random.default_rng(seed)
    n = difference.size
    draws = np.empty(n_boot)
    for index in range(n_boot):
        picked = rng.integers(0, n, size=n)
        draws[index] = difference[picked].mean()
    low, high = np.percentile(draws, [2.5, 97.5])
    baseline = float(errors_b.mean())

    # Non-significant is not equivalent. Protocol section 10 forbids reading a
    # zero-crossing interval as evidence of no effect, so the honest statement
    # is whether the interval can exclude a benefit of the PRESPECIFIED size.
    best_case_relative = low / baseline if baseline else None
    excludes_prespecified_benefit = (
        None
        if best_case_relative is None
        else bool(best_case_relative > -RELATIVE_IMPROVEMENT_THRESHOLD)
    )
    return {
        "n_independent_patients": int(n),
        "mean_difference": point,
        "relative_change_vs_baseline": point / baseline if baseline else None,
        "descriptive_95_interval": [float(low), float(high)],
        "interval_excludes_zero": bool(low * high > 0),
        "direction": "improvement" if point < 0 else "degradation",
        "best_case_relative_in_interval": best_case_relative,
        "excludes_a_prespecified_size_benefit": excludes_prespecified_benefit,
        "how_to_read": (
            "interval_excludes_zero false means no effect was DETECTED, never "
            "that there is none. excludes_a_prespecified_size_benefit true "
            f"means a benefit of at least {RELATIVE_IMPROVEMENT_THRESHOLD:.0%} "
            "is incompatible with the interval; false means the design cannot "
            "rule one out, which is INCONCLUSIVE rather than NOT_SUPPORTED."
        ),
        "n_bootstrap": int(n_boot),
    }


def influence_report(errors_new, errors_baseline, y_true):
    """Mandatory influence reporting, protocol §2.2."""
    import numpy as np

    gain = np.asarray(errors_baseline, dtype=float) - np.asarray(errors_new, dtype=float)
    net = float(gain.sum())
    order = np.argsort(-gain)
    top5 = float(gain[order[:5]].sum())
    quartiles = np.quantile(y_true, [0.25, 0.5, 0.75])
    bucket = np.digitize(y_true, quartiles)
    return {
        "n_improved": int((gain > 0).sum()),
        "n_worsened": int((gain < 0).sum()),
        "net_gain": net,
        "share_of_net_gain_from_five_most_improved": (
            float(top5 / net) if net != 0 else None
        ),
        "mean_gain_by_endpoint_quartile": {
            f"Q{index + 1}": float(gain[bucket == index].mean())
            for index in range(4)
            if (bucket == index).any()
        },
        "reading": (
            "A net gain concentrated in a handful of patients, or present in only "
            "one endpoint quartile, is reported as such and caps the status at "
            "PARTIALLY_SUPPORTED per protocol section 2.1."
        ),
    }


def reversal_check(contrast_primary, contrast_trimmed, contrast_log):
    """Protocol §2.1 reversal rule."""
    improved = contrast_primary["mean_difference"] < 0
    trimmed_reversed = contrast_trimmed["mean_difference"] > 0
    log_reversed = contrast_log["mean_difference"] > 0
    both_reversed = trimmed_reversed and log_reversed
    return {
        "primary_improved": bool(improved),
        "trimmed_mae_reversed": bool(trimmed_reversed),
        "log_mae_reversed": bool(log_reversed),
        "status_capped_at_partially_supported": bool(improved and both_reversed),
        "rule": (
            "If the primary loss improves while BOTH robustness co-criteria "
            "reverse, the status caps at PARTIALLY_SUPPORTED and the influence "
            "analysis becomes mandatory."
        ),
    }


# --- orchestration --------------------------------------------------------
def run(source, target_column=COL_TROUGH, n_boot=N_BOOTSTRAP):
    import numpy as np

    assert_no_prohibited_predictor()
    frame, checksum = load_ariana(source)

    y_true = frame[target_column].to_numpy(dtype=float)
    if np.any(y_true <= 0):
        raise AnalysisError("non-positive concentrations cannot be log-transformed")
    y_log = np.log(y_true)

    folds = make_folds(len(frame), OUTER_FOLDS, OUTER_REPEATS, SEED)
    factories = model_factories()

    report = {
        "objective": "H1",
        "cohort": "ariana_metoprolol_47",
        "source_sha256": checksum,
        "target": target_column,
        "target_is": (
            "observed pre-dose metoprolol concentration. NOT an AUC, NOT a "
            "clearance, NOT a steady-state exposure metric."
        ),
        "prohibited_predictor_absent": COL_PROHIBITED,
        "map_be": "excluded by construction; one observation per patient",
        "status_ceiling": (
            "exploratory. This cohort has been analysed repeatedly by this "
            "project, so no ARIANA result may be reported as confirmatory."
        ),
        "validation": {
            "outer_folds": OUTER_FOLDS,
            "outer_repeats": OUTER_REPEATS,
            "inner_folds": INNER_FOLDS,
            "seed": SEED,
            "folds_identical_across_arms_and_models": True,
            "preprocessing": "StandardScaler fitted inside the training fold only",
            "hyperparameters": "selected by inner CV inside the training fold only",
        },
        "prespecification": {
            "frozen_in_v1_0_0": (
                "arm definitions, the primary contrast and its loss, the loss "
                "hierarchy and reversal rule, the CV scheme, influence reporting, "
                "the prohibited predictor, the MAP-BE exclusion."
            ),
            "operationalised_at_analysis_time": (
                "log-target fitting with back-transformation, the penalty grids, "
                "feature encodings, bootstrap count."
            ),
            "comparators": (
                "Ridge, Elastic Net and XGBoost as the protocol names, with "
                "complexity capped for n=47. No substitution in force."
            ),
        },
        "arms": {name: list(features) for name, features in ARMS.items()},
        "results": {},
    }

    # Out-of-fold predictions for every (model, arm), on identical folds.
    errors = {}
    for model_name, (factory, grid) in factories.items():
        report["results"][model_name] = {}
        for arm_name, features in ARMS.items():
            X = frame[features].to_numpy(dtype=float)
            predictions = out_of_fold_predictions(X, y_log, folds, factory, grid)
            report["results"][model_name][arm_name] = loss_summary(
                y_true, y_log, predictions
            )
            errors[(model_name, arm_name)] = {
                "absolute": per_patient_absolute_error(y_true, predictions),
                "log_absolute": per_patient_log_absolute_error(y_log, predictions),
            }

    # --- the primary contrast --------------------------------------------
    new_arm, baseline_arm = PRIMARY_CONTRAST
    new = errors[(PRIMARY_MODEL, new_arm)]
    baseline = errors[(PRIMARY_MODEL, baseline_arm)]

    primary = paired_contrast(new["absolute"], baseline["absolute"], n_boot=n_boot)
    log_contrast = paired_contrast(
        new["log_absolute"], baseline["log_absolute"], n_boot=n_boot
    )
    # Trimmed co-criterion: the contrast of trimmed losses, not a trimmed contrast.
    trimmed = {
        "mean_difference": trimmed_mean(new["absolute"]) - trimmed_mean(baseline["absolute"]),
        "new_trimmed_mae": trimmed_mean(new["absolute"]),
        "baseline_trimmed_mae": trimmed_mean(baseline["absolute"]),
    }

    relative = primary["relative_change_vs_baseline"]
    report["primary_contrast"] = {
        "definition": f"{PRIMARY_MODEL} {new_arm} versus {baseline_arm}, {PRIMARY_LOSS}",
        "why_ridge": (
            "Ridge estimates its coefficients from the cohort and cannot borrow "
            "strength from literature-fixed constants, so it is the harder test."
        ),
        "primary_loss": primary,
        "robustness_co_criteria": {"log_mae": log_contrast, "trimmed_mae_10pct": trimmed},
        "reversal_rule": reversal_check(
            primary,
            {"mean_difference": trimmed["mean_difference"]},
            log_contrast,
        ),
        "influence": influence_report(new["absolute"], baseline["absolute"], y_true),
        "exploratory_success_criterion": {
            "threshold": f"at least {RELATIVE_IMPROVEMENT_THRESHOLD:.0%} relative MAE "
                         "improvement with a descriptive interval excluding zero",
            "relative_improvement": -relative if relative is not None else None,
            "threshold_met": bool(
                relative is not None
                and -relative >= RELATIVE_IMPROVEMENT_THRESHOLD
                and primary["interval_excludes_zero"]
                and primary["direction"] == "improvement"
            ),
        },
    }

    # Secondary contrasts across every model, same folds, reported as secondary.
    report["secondary_contrasts"] = {}
    for model_name in factories:
        report["secondary_contrasts"][model_name] = {}
        for new_name, base_name in (
            ("M3_clinical_pgx_omic", "M1_clinical_pgx"),
            ("M2_clinical_omic", "M0_clinical"),
            ("M1_clinical_pgx", "M0_clinical"),
            ("M3_clinical_pgx_omic", "M2_clinical_omic"),
            # PGx decomposition: genotype and comedication separately.
            ("M1a_clinical_genotype", "M0_clinical"),
            ("M1b_clinical_comedication", "M0_clinical"),
            ("M3a_clinical_genotype_omic", "M2_clinical_omic"),
            ("M3b_clinical_comedication_omic", "M2_clinical_omic"),
            ("M3_clinical_pgx_omic", "M3a_clinical_genotype_omic"),
        ):
            report["secondary_contrasts"][model_name][f"{new_name}_vs_{base_name}"] = (
                paired_contrast(
                    errors[(model_name, new_name)]["absolute"],
                    errors[(model_name, base_name)]["absolute"],
                    n_boot=n_boot,
                )
            )
    report["multiplicity"] = (
        "One primary contrast, reported unadjusted. All other contrasts are "
        "secondary within the H1 family and require Holm adjustment before any "
        "claim. Fold repetition is not additional independent subjects."
    )
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--source", required=True, help="ARIANA workbook")
    parser.add_argument(
        "--target",
        default=COL_TROUGH,
        choices=[COL_TROUGH, COL_TROUGH_WEIGHTED],
        help="primary target, or the prespecified dose-weighted secondary",
    )
    parser.add_argument("--bootstrap", type=int, default=N_BOOTSTRAP)
    args = parser.parse_args(argv)
    try:
        report = run(args.source, args.target, args.bootstrap)
    except AnalysisError as exc:
        parser.exit(2, f"INVALID: {exc}\n")
    print(json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
