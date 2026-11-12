"""H1 key secondary and hybrid arm: ARIANA mechanistic M5 versus M4, plus HY.

What this runs
--------------
Protocol §4 designates one primary contrast — Ridge M3 versus Ridge M1, already
run by ``run_h1_ariana_incremental_omic_value.py`` — and one **key secondary**:
the mechanistic **M5 versus M4**, same folds, same losses. Protocol §4 also
names ``P3 vs P1`` and ``HY vs P1`` among the other secondaries. This runner
produces those, on the **same folds and the same loss hierarchy** as the
primary, by importing the primary runner's machinery rather than reimplementing
it.

The order is the protocol's and is not revisited in the light of the numbers:
the primary came first, the key secondary is the contrast reported here, and HY
is examined after it. No arm is promoted to primary because it did better.

What is fixed before the run
----------------------------
* **Primary backbone**: the extended-release two-compartment scenario, on the
  pharmaceutical argument that ARIANA's 47.5 and 95 mg unit doses are the
  German succinate extended-release presentations. The two immediate-release
  backbones run as structural sensitivity.
* **Primary omic panel**: 4-OH-solanidine alone, the mechanistically
  prespecified panel, because it has an authentic synthesised standard and SSDA
  does not. The prior ARIANA bundle promoted ``OH + SSDA`` to primary **on the
  observed loss**; that choice is not admissible as a primary here and runs as
  a declared sensitivity.
* **Transferred constants**: f_m = 0.80, gamma = 1, genotype/omic weight 0.5,
  activity-score decrement 0.54 per inhibitor. Each is swept afterwards, and
  the sweep is reported whether or not it favours the fixed value.
* **HY**: the published-coefficient PopPK structural prediction, plus a ridge
  model of its log residual fitted on the training patients of the outer fold,
  with the penalty chosen by inner cross-validation inside that fold and the
  clinical + PGx + omic design of the primary M3 arm. A variant whose residual
  stage sees only the omic block is reported as a secondary.

What this runner refuses to do
------------------------------
* No arm uses ``lnMR OH-metoprolol/metoprolol``: it shares its denominator
  with the target. A preflight check enforces it across every design.
* No MAP-Bayesian step: one concentration per subject makes conditioning on it
  and then predicting it circular.
* No NLME fit is claimed.
* No individual prediction is printed; only aggregates.
* Nothing here may be reported as confirmatory. ARIANA has been analysed
  repeatedly by this project and the status ceiling is ``exploratory``.

Usage
-----
    python scripts/run_h1_ariana_mechanistic_and_hybrid_arms.py \
        --source PATH --out reports/data/h1_ariana_mechanistic_hybrid.json
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from omics2pk import ariana_mechanistic_arms as MECH  # noqa: E402

PROTOCOL_VERSION = "1.4.7"
REPORT_VERSION = "1.0.0"


def load_primary_runner():
    """Import the primary H1 runner so the folds and losses are literally its own."""
    path = ROOT / "scripts" / "run_h1_ariana_incremental_omic_value.py"
    spec = importlib.util.spec_from_file_location("run_h1_primary", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


H1 = load_primary_runner()


def check_primary_machinery_unchanged():
    """The key secondary is only on the same folds if these are still the same."""
    expected = {"OUTER_FOLDS": 5, "OUTER_REPEATS": 2, "INNER_FOLDS": 3,
                "SEED": 20261107, "N_BOOTSTRAP": 2000,
                "RELATIVE_IMPROVEMENT_THRESHOLD": 0.10, "TRIM_FRACTION": 0.10}
    mismatched = {
        name: (value, getattr(H1, name))
        for name, value in expected.items() if getattr(H1, name) != value
    }
    if mismatched:
        raise H1.AnalysisError(
            "the primary runner's frozen validation constants have changed, so "
            f"these arms would not share its folds: {mismatched}"
        )


def fold_fingerprint(folds):
    """A hash of the fold membership, so fold identity is auditable from the report."""
    digest = hashlib.sha256()
    for train, test in folds:
        digest.update(b"|")
        digest.update(",".join(str(int(i)) for i in train).encode())
        digest.update(b";")
        digest.update(",".join(str(int(i)) for i in test).encode())
    return digest.hexdigest()


# --- designs ---------------------------------------------------------------
#: The residual stage of HY uses the primary M3 design.
HY_RESIDUAL_DESIGNS = {
    "clinical_pgx_omic": H1.CLINICAL + H1.PGX + H1.OMIC,
    "omic_only": list(H1.OMIC),
}
HY_PRIMARY_RESIDUAL_DESIGN = "clinical_pgx_omic"


def assert_no_prohibited_predictor():
    """``lnMR OH-metoprolol/metoprolol`` must appear in no design, anywhere."""
    for name, features in HY_RESIDUAL_DESIGNS.items():
        if H1.COL_PROHIBITED in features:
            raise H1.AnalysisError(
                f"HY residual design {name!r} contains the prohibited "
                f"predictor {H1.COL_PROHIBITED!r}"
            )
    for panel, columns in MECH.OMIC_PANELS.items():
        if H1.COL_PROHIBITED in columns:
            raise H1.AnalysisError(
                f"omic panel {panel!r} contains the prohibited predictor"
            )
    H1.assert_no_prohibited_predictor()


def build_arrays(frame, panel_columns):
    import numpy as np

    missing = [column for column in panel_columns if column not in frame.columns]
    if missing:
        raise H1.AnalysisError(f"omic panel columns absent from source: {missing}")
    return {
        "observed": frame[H1.COL_TROUGH].to_numpy(dtype=float),
        "activity_score": frame[H1.COL_AS].to_numpy(dtype=float),
        "inhibitor": frame[H1.COL_INHIBITOR].to_numpy(dtype=float),
        "weight": frame["weight"].to_numpy(dtype=float),
        "male": frame["male sex"].to_numpy(dtype=float),
        "daily_dose_mg": frame[H1.COL_DOSE].to_numpy(dtype=float),
        "ratios": np.column_stack(
            [frame[column].to_numpy(dtype=float) for column in panel_columns]
        ),
    }


def subset(arrays, index):
    return {
        key: (value[index] if key != "ratios" else value[index, :])
        for key, value in arrays.items()
    }


# --- out-of-fold predictions ----------------------------------------------
def mechanistic_out_of_fold(arrays, folds, backbone, *, gamma, blend_weight,
                            fm, delta_as, dose_equivalence,
                            arms=None):
    """Log-scale out-of-fold predictions for the mechanistic arms.

    Returns ``(predictions, trace)`` where ``predictions[arm]`` has the shape
    the primary runner's loss functions expect, ``(n_repeats, n_patients)``.
    """
    import numpy as np

    arms = dict(arms or MECH.MECHANISTIC_ARMS)
    n = arrays["observed"].size
    n_repeats = len(folds) // H1.OUTER_FOLDS
    predictions = {
        name: np.full((n_repeats, n), np.nan) for name in arms
    }
    locked, reference_sizes, relaxed = {name: [] for name in arms}, [], 0

    for fold_index, (train_index, test_index) in enumerate(folds):
        repeat = fold_index // H1.OUTER_FOLDS
        train = subset(arrays, train_index)
        test = subset(arrays, test_index)
        reference = MECH.omic_reference(
            train["activity_score"], train["inhibitor"], train["ratios"]
        )
        reference_sizes.append(reference["n_reference_patients"])
        relaxed += int(reference["inhibitor_free_requirement_relaxed"])
        shared = dict(gamma=gamma, blend_weight=blend_weight, fm=fm,
                      delta_as=delta_as, dose_equivalence=dose_equivalence)
        for name, arm in arms.items():
            baseline = MECH.lock_baseline_clearance(
                train, arm, backbone, reference["reference"], **shared
            )
            locked[name].append(baseline)
            predicted = MECH.structural_prediction(
                test, baseline, arm, backbone, reference["reference"], **shared
            )
            predictions[name][repeat, test_index] = np.log(predicted)

    for name, matrix in predictions.items():
        if np.isnan(matrix).any():
            raise H1.AnalysisError(
                f"arm {name!r} left a patient without an out-of-fold prediction"
            )
    trace = {
        "locked_baseline_clearance_L_per_h": {
            name: {
                "median": float(np.median(values)),
                "min": float(np.min(values)),
                "max": float(np.max(values)),
            }
            for name, values in locked.items()
        },
        "omic_reference_patients_per_fold": {
            "median": float(np.median(reference_sizes)),
            "min": int(np.min(reference_sizes)),
            "max": int(np.max(reference_sizes)),
        },
        "folds_where_the_inhibitor_free_requirement_was_relaxed": int(relaxed),
        "what_is_fitted_on_the_cohort": (
            "one parameter per arm per fold: the reference apparent clearance, "
            "by least squares on the training patients' log residuals, locked "
            "before the test fold is touched"
        ),
    }
    return predictions, trace


def poppk_out_of_fold(arrays, folds, backbone, *, delta_as):
    """Log-scale out-of-fold predictions for P1, P2 and P3."""
    import numpy as np

    n = arrays["observed"].size
    n_repeats = len(folds) // H1.OUTER_FOLDS
    names = ("P1_poppk_published_coefficients",
             "P2_poppk_estimated_coefficients",
             "P3_poppk_estimated_plus_omic")
    predictions = {name: np.full((n_repeats, n), np.nan) for name in names}
    estimated, failures, cross_checks = [], 0, []

    for fold_index, (train_index, test_index) in enumerate(folds):
        repeat = fold_index // H1.OUTER_FOLDS
        train = subset(arrays, train_index)
        test = subset(arrays, test_index)
        reference = MECH.omic_reference(
            train["activity_score"], train["inhibitor"], train["ratios"]
        )["reference"]
        design_train = MECH.poppk_design(train, reference, delta_as)
        design_test = MECH.poppk_design(test, reference, delta_as)

        # P1: published coefficients, one free parameter, locked exactly as the
        # mechanistic arms lock theirs so that the two are comparable.
        baseline = MECH.lock_baseline_clearance(
            train, "pgx", backbone, reference, gamma=MECH.OMIC_GAMMA,
            blend_weight=MECH.BLEND_WEIGHT, fm=MECH.FM_CYP2D6,
            delta_as=delta_as, dose_equivalence=1.0,
        )
        theta_published = {"log_cl0": float(np.log(baseline))}
        predictions[names[0]][repeat, test_index] = np.log(
            MECH.poppk_predict(test, design_test, theta_published, backbone,
                               published=True)
        )

        for name, with_omic in ((names[1], False), (names[2], True)):
            fit = MECH.poppk_fit(train, design_train, backbone, with_omic)
            failures += int(not fit["converged"])
            cross_checks.append(fit["cross_check_max_parameter_difference"])
            predictions[name][repeat, test_index] = np.log(
                MECH.poppk_predict(test, design_test, fit["theta"], backbone,
                                   published=False)
            )
            estimated.append({"arm": name, "repeat": repeat,
                              "fold": fold_index % H1.OUTER_FOLDS, **fit["theta"]})

    for name, matrix in predictions.items():
        if np.isnan(matrix).any():
            raise H1.AnalysisError(f"arm {name!r} left a patient unpredicted")

    def spread(arm, key):
        values = [row[key] for row in estimated if row["arm"] == arm and key in row]
        if not values:
            return None
        array = np.asarray(values, dtype=float)
        mean = float(array.mean())
        return {
            "median": float(np.median(array)),
            "q1": float(np.percentile(array, 25)),
            "q3": float(np.percentile(array, 75)),
            "cv_percent": (
                float(100.0 * array.std(ddof=1) / abs(mean)) if mean else None
            ),
        }

    trace = {
        "fits_failing_the_convergence_cross_check": int(failures),
        "convergence_cross_check": {
            "how": (
                "each fold's gradient-based fit is refitted with a "
                "derivative-free method from its own solution; the largest "
                "disagreement in any coefficient is recorded. A solver's own "
                "success flag is not used as evidence: the objective is "
                "differentiated by finite differences, so a tight gradient "
                "tolerance is unreachable at a perfectly good optimum."
            ),
            "max_parameter_disagreement_over_all_fits": (
                float(np.max(cross_checks)) if cross_checks else None
            ),
            "tolerance": 1e-4,
        },
        "estimated_coefficients_across_folds": {
            arm: {
                key: spread(arm, key)
                for key in ("log_cl0", "sex", "activity_score", "omic")
                if spread(arm, key) is not None
            }
            for arm in (names[1], names[2])
        },
        "published_values_for_comparison": {
            "sex": MECH.SEX_LOG_EFFECT,
            "activity_score": (
                "not a free power in the published arm: the mechanistic "
                f"relation (1 - f_m) + f_m * a_rel with f_m = {MECH.FM_CYP2D6}"
            ),
            "omic": (
                "not published as a clearance coefficient; its mechanistic "
                f"counterpart is f_m * gamma = {MECH.FM_CYP2D6 * MECH.OMIC_GAMMA}"
            ),
        },
        "how_to_read_a_wide_cv": (
            "A coefficient whose fold-to-fold CV exceeds 100 % is not "
            "identified at this sample size. That is an argument for fixing it "
            "from the literature, not a result about the covariate."
        ),
    }
    return predictions, trace


def hybrid_out_of_fold(arrays, frame, folds, backbone, *, delta_as,
                       design_name):
    """HY: the published-coefficient structural prediction plus a ridge residual."""
    import numpy as np
    from sklearn.linear_model import Ridge
    from sklearn.model_selection import KFold
    from sklearn.preprocessing import StandardScaler

    features = HY_RESIDUAL_DESIGNS[design_name]
    design_matrix = frame[features].to_numpy(dtype=float)
    alphas = list(np.logspace(-3, 3, 25))

    n = arrays["observed"].size
    n_repeats = len(folds) // H1.OUTER_FOLDS
    predictions = np.full((n_repeats, n), np.nan)
    chosen_alphas, residual_sd = [], []

    for fold_index, (train_index, test_index) in enumerate(folds):
        repeat = fold_index // H1.OUTER_FOLDS
        train = subset(arrays, train_index)
        test = subset(arrays, test_index)
        reference = MECH.omic_reference(
            train["activity_score"], train["inhibitor"], train["ratios"]
        )["reference"]

        baseline = MECH.lock_baseline_clearance(
            train, "pgx", backbone, reference, gamma=MECH.OMIC_GAMMA,
            blend_weight=MECH.BLEND_WEIGHT, fm=MECH.FM_CYP2D6,
            delta_as=delta_as, dose_equivalence=1.0,
        )
        structural_train = MECH.structural_prediction(
            train, baseline, "pgx", backbone, reference, gamma=MECH.OMIC_GAMMA,
            blend_weight=MECH.BLEND_WEIGHT, fm=MECH.FM_CYP2D6,
            delta_as=delta_as, dose_equivalence=1.0,
        )
        structural_test = MECH.structural_prediction(
            test, baseline, "pgx", backbone, reference, gamma=MECH.OMIC_GAMMA,
            blend_weight=MECH.BLEND_WEIGHT, fm=MECH.FM_CYP2D6,
            delta_as=delta_as, dose_equivalence=1.0,
        )
        residual = np.log(train["observed"]) - np.log(structural_train)
        residual_sd.append(float(residual.std(ddof=1)))

        X_train = design_matrix[train_index]
        X_test = design_matrix[test_index]

        # Penalty chosen by inner CV inside the training fold only, with the
        # same inner scheme the primary runner uses.
        inner = KFold(n_splits=H1.INNER_FOLDS, shuffle=True,
                      random_state=H1.SEED + repeat)
        best_alpha, best_score = None, np.inf
        for alpha in alphas:
            errors = []
            for inner_train, inner_validation in inner.split(X_train):
                scaler = StandardScaler().fit(X_train[inner_train])
                model = Ridge(alpha=alpha)
                model.fit(scaler.transform(X_train[inner_train]),
                          residual[inner_train])
                predicted = model.predict(
                    scaler.transform(X_train[inner_validation])
                )
                errors.append(
                    float(np.mean((predicted - residual[inner_validation]) ** 2))
                )
            score = float(np.mean(errors))
            if score < best_score:
                best_alpha, best_score = alpha, score
        chosen_alphas.append(float(best_alpha))

        scaler = StandardScaler().fit(X_train)
        model = Ridge(alpha=best_alpha)
        model.fit(scaler.transform(X_train), residual)
        correction = model.predict(scaler.transform(X_test))
        predictions[repeat, test_index] = np.log(structural_test) + correction

    if np.isnan(predictions).any():
        raise H1.AnalysisError("HY left a patient without a prediction")
    trace = {
        "residual_design": features,
        "structural_base": "P1, the published-coefficient PopPK arm",
        "residual_target": (
            "log observed minus log structural prediction, computed on the "
            "training patients of the outer fold only"
        ),
        "ridge_penalty_chosen_per_fold": {
            "median": float(np.median(chosen_alphas)),
            "min": float(np.min(chosen_alphas)),
            "max": float(np.max(chosen_alphas)),
        },
        "training_log_residual_sd_per_fold": {
            "median": float(np.median(residual_sd)),
            "min": float(np.min(residual_sd)),
            "max": float(np.max(residual_sd)),
        },
    }
    return predictions, trace


# --- losses and contrasts -------------------------------------------------
def arm_errors(observed, observed_log, predictions):
    return {
        "absolute": H1.per_patient_absolute_error(observed, predictions),
        "log_absolute": H1.per_patient_log_absolute_error(observed_log, predictions),
    }


def bonferroni_interval(errors_new, errors_baseline, n_boot, family_size,
                        alpha=0.05, seed=None):
    """The same paired bootstrap, read at a level corrected for the family.

    Protocol §9 requires every contrast outside the primary and the designated
    key secondary to be adjusted for multiplicity before any claim. Saying so
    in prose is not enough, so the corrected interval is computed: the draws
    are the primary runner's own, only the percentiles move.
    """
    import numpy as np

    difference = (
        np.asarray(errors_new, dtype=float)
        - np.asarray(errors_baseline, dtype=float)
    )
    rng = np.random.default_rng(H1.SEED if seed is None else seed)
    n = difference.size
    draws = np.empty(n_boot)
    for index in range(n_boot):
        draws[index] = difference[rng.integers(0, n, size=n)].mean()
    half = 100.0 * alpha / (2.0 * family_size)
    low, high = np.percentile(draws, [half, 100.0 - half])
    baseline = float(np.asarray(errors_baseline, dtype=float).mean())
    return {
        "family_size": int(family_size),
        "level": f"{100.0 * (1.0 - alpha / family_size):.2f} %",
        "interval": [float(low), float(high)],
        "interval_excludes_zero": bool(low * high > 0),
        "relative_interval_vs_baseline": (
            [float(low / baseline), float(high / baseline)] if baseline else None
        ),
        "method": (
            "Bonferroni across the contrasts of this family, applied to the "
            "percentiles of the same paired patient bootstrap"
        ),
    }


def full_contrast(errors_new, errors_baseline, observed, n_boot):
    """The protocol's loss hierarchy for one contrast, with its reversal rule."""
    primary = H1.paired_contrast(errors_new["absolute"],
                                 errors_baseline["absolute"], n_boot=n_boot)
    log_contrast = H1.paired_contrast(errors_new["log_absolute"],
                                      errors_baseline["log_absolute"],
                                      n_boot=n_boot)
    trimmed_new = H1.trimmed_mean(errors_new["absolute"])
    trimmed_baseline = H1.trimmed_mean(errors_baseline["absolute"])
    trimmed = {
        "mean_difference": trimmed_new - trimmed_baseline,
        "new_trimmed_mae": trimmed_new,
        "baseline_trimmed_mae": trimmed_baseline,
    }
    relative = primary["relative_change_vs_baseline"]
    return {
        "primary_loss": primary,
        "robustness_co_criteria": {"log_mae": log_contrast,
                                   "trimmed_mae_10pct": trimmed},
        "reversal_rule": H1.reversal_check(
            primary, {"mean_difference": trimmed["mean_difference"]}, log_contrast
        ),
        "influence": H1.influence_report(errors_new["absolute"],
                                         errors_baseline["absolute"], observed),
        "exploratory_success_criterion": {
            "threshold": (
                f"at least {H1.RELATIVE_IMPROVEMENT_THRESHOLD:.0%} relative MAE "
                "improvement with a descriptive 95 % interval excluding zero, "
                "and no reversal on both robustness co-criteria"
            ),
            "relative_improvement": -relative if relative is not None else None,
            "threshold_met": bool(
                relative is not None
                and -relative >= H1.RELATIVE_IMPROVEMENT_THRESHOLD
                and primary["interval_excludes_zero"]
                and primary["direction"] == "improvement"
            ),
        },
    }


# --- sensitivity ----------------------------------------------------------
#: Declared before the run. Each sweep varies one transferred constant and
#: leaves the others at their published value.
SENSITIVITY_SWEEPS = {
    "fraction_metabolised_by_cyp2d6": ("fm", (0.60, 0.70, 0.80, 0.90)),
    "omic_exponent_gamma": ("gamma", (0.25, 0.50, 0.75, 1.00)),
    "genotype_omic_weight": ("blend_weight", (0.00, 0.25, 0.50, 0.75, 1.00)),
    "activity_score_decrement_per_inhibitor": (
        "delta_as", (0.00, 0.54, 1.00, 3.00)),
    "salt_dose_equivalence": (
        "dose_equivalence", (1.0, MECH.SUCCINATE_TARTRATE_EQUIVALENCE)),
}


def key_secondary_only(arrays, folds, backbone, observed, observed_log, n_boot,
                       **overrides):
    """MAE of M4 and M5 and their contrast, for a sensitivity cell."""
    settings = dict(gamma=MECH.OMIC_GAMMA, blend_weight=MECH.BLEND_WEIGHT,
                    fm=MECH.FM_CYP2D6, delta_as=MECH.DELTA_AS_PER_INHIBITOR,
                    dose_equivalence=1.0)
    settings.update(overrides)
    new_arm, baseline_arm = MECH.KEY_SECONDARY_CONTRAST
    arms = {name: MECH.MECHANISTIC_ARMS[name] for name in (new_arm, baseline_arm)}
    predictions, _ = mechanistic_out_of_fold(arrays, folds, backbone,
                                             arms=arms, **settings)
    errors = {
        name: arm_errors(observed, observed_log, matrix)
        for name, matrix in predictions.items()
    }
    contrast = H1.paired_contrast(errors[new_arm]["absolute"],
                                 errors[baseline_arm]["absolute"], n_boot=n_boot)
    return {
        "M4_mae": float(errors[baseline_arm]["absolute"].mean()),
        "M5_mae": float(errors[new_arm]["absolute"].mean()),
        "M4_log_mae": float(errors[baseline_arm]["log_absolute"].mean()),
        "M5_log_mae": float(errors[new_arm]["log_absolute"].mean()),
        "relative_change": contrast["relative_change_vs_baseline"],
        "interval_excludes_zero": contrast["interval_excludes_zero"],
    }


# --- orchestration --------------------------------------------------------
def run(source, n_boot=None, backbone_name=MECH.PRIMARY_BACKBONE,
        panel_name=MECH.PRIMARY_OMIC_PANEL, sensitivity=True,
        sensitivity_bootstrap=200):
    import numpy as np

    check_primary_machinery_unchanged()
    assert_no_prohibited_predictor()
    n_boot = H1.N_BOOTSTRAP if n_boot is None else n_boot

    frame, checksum = H1.load_ariana(source)
    observed = frame[H1.COL_TROUGH].to_numpy(dtype=float)
    if np.any(observed <= 0):
        raise H1.AnalysisError("non-positive concentrations cannot be logged")
    observed_log = np.log(observed)

    folds = H1.make_folds(len(frame), H1.OUTER_FOLDS, H1.OUTER_REPEATS, H1.SEED)
    backbone = MECH.BACKBONES[backbone_name]
    panel_columns = MECH.OMIC_PANELS[panel_name]
    arrays = build_arrays(frame, panel_columns)

    report = {
        "objective": "H1 key secondary and hybrid arm",
        "report_version": REPORT_VERSION,
        "protocol_version": PROTOCOL_VERSION,
        "cohort": "ariana_metoprolol_47",
        "source_sha256": checksum,
        "target": H1.COL_TROUGH,
        "target_is": (
            "observed pre-dose metoprolol concentration, ng/mL. NOT an AUC, "
            "NOT a clearance, NOT a steady-state exposure metric."
        ),
        "prohibited_predictor_absent": H1.COL_PROHIBITED,
        "map_be": (
            "excluded by construction: one observation per patient makes "
            "conditioning on it and predicting it circular"
        ),
        "nlme": (
            "not claimed: between-subject and residual variability are not "
            "separately identifiable from one observation per subject"
        ),
        "status_ceiling": (
            "exploratory. This cohort has been analysed repeatedly by this "
            "project, so no ARIANA result may be reported as confirmatory, "
            "however favourable it is."
        ),
        "validation": {
            "outer_folds": H1.OUTER_FOLDS,
            "outer_repeats": H1.OUTER_REPEATS,
            "inner_folds": H1.INNER_FOLDS,
            "seed": H1.SEED,
            "folds_identical_to_the_primary_contrast": True,
            "fold_fingerprint_sha256": fold_fingerprint(folds),
            "how_identity_is_guaranteed": (
                "the fold list comes from the primary runner's own "
                "make_folds, imported rather than reimplemented, and its "
                "frozen constants are checked before the run"
            ),
            "bootstrap": {
                "n": int(n_boot),
                "unit": "patients; fold repetition is never additional subjects",
            },
        },
        "configuration": {
            "backbone": {"key": backbone_name, **{
                key: value for key, value in backbone.items()}},
            "omic_panel": {"key": panel_name, "columns": list(panel_columns)},
            "why_this_panel_is_primary": (
                "4-OH-solanidine has an authentic synthesised standard; SSDA "
                "is semi-quantitative without one. The prior ARIANA bundle "
                "promoted OH+SSDA to primary on the observed loss, which is a "
                "result-guided choice and is carried here as a sensitivity."
            ),
            "why_this_backbone_is_primary": (
                "47.5 and 95 mg unit doses are the German succinate "
                "extended-release presentations. The argument is "
                "pharmaceutical, never the observed loss; both "
                "immediate-release backbones run as structural sensitivity."
            ),
            "dosing_assumption": (
                "two equal daily doses, tau = 12 h, trough drawn at 12 h. The "
                "source documents neither clock time nor time since the last "
                "dose, so this is an assumption, not a measurement."
            ),
        },
        "literature_constants": [
            {"module": module, "fixed_coefficient": coefficient,
             "source": source_text, "estimated_on_ariana": False}
            for module, coefficient, source_text in MECH.LITERATURE_CONSTANTS
        ],
        "free_parameters_fitted_on_this_cohort": dict(MECH.FREE_PARAMETERS),
        "activity_score_alignment": MECH.activity_score_grid_audit(
            arrays["activity_score"]
        ),
        "arms": {},
        "losses": {},
    }

    # --- the arms ---------------------------------------------------------
    mechanistic_predictions, mechanistic_trace = mechanistic_out_of_fold(
        arrays, folds, backbone, gamma=MECH.OMIC_GAMMA,
        blend_weight=MECH.BLEND_WEIGHT, fm=MECH.FM_CYP2D6,
        delta_as=MECH.DELTA_AS_PER_INHIBITOR, dose_equivalence=1.0,
    )
    poppk_predictions, poppk_trace = poppk_out_of_fold(
        arrays, folds, backbone, delta_as=MECH.DELTA_AS_PER_INHIBITOR
    )
    hybrid_predictions, hybrid_trace = hybrid_out_of_fold(
        arrays, frame, folds, backbone,
        delta_as=MECH.DELTA_AS_PER_INHIBITOR,
        design_name=HY_PRIMARY_RESIDUAL_DESIGN,
    )
    hybrid_secondary, hybrid_secondary_trace = hybrid_out_of_fold(
        arrays, frame, folds, backbone,
        delta_as=MECH.DELTA_AS_PER_INHIBITOR, design_name="omic_only",
    )

    predictions = dict(mechanistic_predictions)
    predictions.update(poppk_predictions)
    predictions["HY_poppk_published_plus_ml_residual"] = hybrid_predictions
    predictions["HYc_residual_on_the_omic_block_only"] = hybrid_secondary

    report["arms"] = {
        "mechanistic": mechanistic_trace,
        "poppk_covariate": poppk_trace,
        "hybrid": hybrid_trace,
        "hybrid_secondary": hybrid_secondary_trace,
        "what_each_arm_is": {
            "M4c_backbone_only": (
                "structural control: no CYP2D6 information at all. Not a "
                "protocol arm; it exists so that M4's own content is readable."
            ),
            "M4_mechanistic_pgx": (
                "protocol M4: mechanistic structure plus genotype and "
                "phenoconversion, all coefficients literature-fixed"
            ),
            "M5c_mechanistic_omic_only": (
                "the omic module substituting for the genotype. Secondary."
            ),
            "M5_mechanistic_pgx_omic": (
                "protocol M5: M4 plus the omic module, geometric blend at the "
                "fixed weight 0.5"
            ),
            "P1_poppk_published_coefficients": (
                "protocol P1. At ARIANA's covariate set this is the SAME model "
                "as M4: identical structure, identical published coefficients, "
                "identical single free parameter. The report measures their "
                "agreement rather than asserting it."
            ),
            "P2_poppk_estimated_coefficients": (
                "protocol P2: the same covariates with the sex and CYP2D6 "
                "coefficients estimated on the training fold"
            ),
            "P3_poppk_estimated_plus_omic": (
                "protocol P3: P2 plus the omic as an estimated coefficient"
            ),
            "HY_poppk_published_plus_ml_residual": (
                "protocol HY: P1's structural prediction plus a ridge model of "
                "its log residual, fitted inside the training fold"
            ),
            "HYc_residual_on_the_omic_block_only": (
                "HY variant whose residual stage sees only the two solanidine "
                "ratios. Secondary, declared before the run."
            ),
        },
    }

    errors = {
        name: arm_errors(observed, observed_log, matrix)
        for name, matrix in predictions.items()
    }
    report["losses"] = {
        name: H1.loss_summary(observed, observed_log, matrix)
        for name, matrix in predictions.items()
    }

    # P1 is M4 under the protocol's own definitions: measure, do not assert.
    difference = np.abs(
        predictions["P1_poppk_published_coefficients"]
        - predictions["M4_mechanistic_pgx"]
    )
    report["p1_and_m4_are_the_same_model"] = {
        "claim": (
            "Protocol §3 lists P1 (PopPK with published covariate "
            "coefficients, 1 free parameter) and M4 (fixed-coefficient "
            "mechanistic model plus PGx) as separate arms. At ARIANA's "
            "covariate set they are the same model, so this is a property of "
            "the arm table and not a coding shortcut."
        ),
        "max_absolute_log_difference": float(difference.max()),
        "two_independent_code_paths": True,
        "consequence": (
            "'HY versus P1' is therefore also 'HY versus M4', and P3 versus P1 "
            "is the fix-versus-estimate comparison at identical structure."
        ),
    }

    # --- the key secondary contrast ---------------------------------------
    new_arm, baseline_arm = MECH.KEY_SECONDARY_CONTRAST
    report["key_secondary_contrast"] = {
        "definition": (
            f"{new_arm} versus {baseline_arm}, mean absolute error on the "
            "observed trough, identical folds to the primary contrast"
        ),
        "designated_in": "protocol §4, key secondary",
        "reported": "unadjusted and labelled secondary, per protocol §9",
        **full_contrast(errors[new_arm], errors[baseline_arm], observed, n_boot),
    }

    # --- the other secondaries protocol §4 names --------------------------
    report["other_secondary_contrasts"] = {}
    other_secondaries = {
        "P3_vs_P1": ("P3_poppk_estimated_plus_omic",
                     "P1_poppk_published_coefficients"),
        "HY_vs_P1": ("HY_poppk_published_plus_ml_residual",
                     "P1_poppk_published_coefficients"),
        "P2_vs_P1": ("P2_poppk_estimated_coefficients",
                     "P1_poppk_published_coefficients"),
        "M4_vs_backbone_only": ("M4_mechanistic_pgx", "M4c_backbone_only"),
        "M5c_omic_only_vs_M4": ("M5c_mechanistic_omic_only",
                                "M4_mechanistic_pgx"),
        "HYc_omic_residual_vs_P1": ("HYc_residual_on_the_omic_block_only",
                                    "P1_poppk_published_coefficients"),
        "HY_vs_M5": ("HY_poppk_published_plus_ml_residual",
                     "M5_mechanistic_pgx_omic"),
    }
    for label, (new_name, base_name) in other_secondaries.items():
        contrast = full_contrast(errors[new_name], errors[base_name], observed,
                                 n_boot)
        contrast["multiplicity_adjusted"] = bonferroni_interval(
            errors[new_name]["absolute"], errors[base_name]["absolute"],
            n_boot, len(other_secondaries),
        )
        report["other_secondary_contrasts"][label] = contrast
    report["other_secondary_contrasts"]["HY_vs_M5"]["not_a_fair_arm_comparison"] = (
        "HY's residual stage uses the primary contrast's M3 design, which "
        f"contains both solanidine ratios, while M5 uses the primary panel "
        f"{list(panel_columns)}. The two arms therefore see different amounts "
        "of omic information and this contrast is descriptive only. It is "
        "reported because the design was declared before the run, not because "
        "it answers which arm is better."
    )

    # --- sensitivity ------------------------------------------------------
    if sensitivity:
        report["structural_sensitivity"] = {}
        for name in MECH.BACKBONES:
            report["structural_sensitivity"][name] = key_secondary_only(
                arrays, folds, MECH.BACKBONES[name], observed, observed_log,
                sensitivity_bootstrap,
            )
        report["omic_panel_sensitivity"] = {}
        for name, columns in MECH.OMIC_PANELS.items():
            report["omic_panel_sensitivity"][name] = key_secondary_only(
                build_arrays(frame, columns), folds, backbone, observed,
                observed_log, sensitivity_bootstrap,
            )
        report["transferred_constant_sensitivity"] = {}
        for label, (keyword, values) in SENSITIVITY_SWEEPS.items():
            report["transferred_constant_sensitivity"][label] = {
                f"{value:g}": key_secondary_only(
                    arrays, folds, backbone, observed, observed_log,
                    sensitivity_bootstrap, **{keyword: value}
                )
                for value in values
            }
        report["how_to_read_the_sensitivity"] = (
            "These sweeps use a reduced bootstrap and are descriptive only. "
            "Their purpose is to show whether the key secondary's direction "
            "survives the transfer assumptions, not to choose a value: the "
            "published values were fixed before the sweep and are not changed "
            "after it."
        )

    report["multiplicity"] = {
        "rule": (
            "One primary contrast for the whole programme, already reported "
            "elsewhere. The key secondary here is reported unadjusted and "
            "labelled secondary, as protocol §9 prescribes; every other "
            "contrast is exploratory and carries a Bonferroni-adjusted "
            "interval alongside its nominal one. Fold repetition is never "
            "additional independent subjects."
        ),
        "other_secondary_family_size": len(other_secondaries),
        "which_other_secondaries_survive_adjustment": sorted(
            label for label, contrast in report["other_secondary_contrasts"].items()
            if contrast["multiplicity_adjusted"]["interval_excludes_zero"]
        ),
        "which_lose_their_interval_on_adjustment": sorted(
            label for label, contrast in report["other_secondary_contrasts"].items()
            if contrast["primary_loss"]["interval_excludes_zero"]
            and not contrast["multiplicity_adjusted"]["interval_excludes_zero"]
        ),
    }
    report["what_this_cannot_establish"] = [
        "No replication cohort exists for the solanidine axis in public data, "
        "so no result here can be upgraded beyond exploratory.",
        "The galenic form and the exact time since the last dose are not "
        "documented in the source; every prediction assumes tau = 12 h.",
        "SSDA is semi-quantitative without an authentic analytical standard, "
        "so any panel containing it carries an unquantified identification "
        "uncertainty.",
        "The omic module's validation against the OH-metoprolol phenotype is "
        "internal to these 47 patients, one laboratory and one assay.",
        f"{report['activity_score_alignment']['n_patients_interpolated']} of "
        f"{observed.size} patients receive an interpolated relative activity "
        "that the published table does not contain.",
    ]
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--source", required=True, help="ARIANA workbook")
    parser.add_argument("--backbone", default=MECH.PRIMARY_BACKBONE,
                        choices=sorted(MECH.BACKBONES))
    parser.add_argument("--omic-panel", default=MECH.PRIMARY_OMIC_PANEL,
                        choices=sorted(MECH.OMIC_PANELS))
    parser.add_argument("--bootstrap", type=int, default=H1.N_BOOTSTRAP)
    parser.add_argument("--sensitivity-bootstrap", type=int, default=200)
    parser.add_argument("--no-sensitivity", action="store_true")
    parser.add_argument("--out", help="write the report here instead of stdout")
    args = parser.parse_args(argv)
    try:
        report = run(args.source, args.bootstrap, args.backbone,
                     args.omic_panel, not args.no_sensitivity,
                     args.sensitivity_bootstrap)
    except (H1.AnalysisError, MECH.MechanisticError) as exc:
        parser.exit(2, f"INVALID: {exc}\n")
    text = json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True)
    if args.out:
        Path(args.out).write_text(text + "\n", encoding="utf-8")
        print(f"written: {args.out}")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
