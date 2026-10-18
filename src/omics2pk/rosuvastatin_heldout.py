"""Private, exploratory E3 rosuvastatin n=40 held-human-out benchmark.

Only operates on a locally supplied, author-verified, one-row-per-human
predose molecular CSV with genuinely quantified AUC and Cmax. A structural
CSV match is NOT scientific authentication. No per-person predictions,
identifiers, names or feature values are emitted. No model/dose deployment.
"""
from __future__ import annotations

from math import exp, fsum, isfinite, log, log1p, sqrt
from pathlib import Path
from random import Random
from statistics import mean, median
from typing import Mapping, Sequence

from .rosuvastatin_original_local import (
    EMPTY_TOKENS, PK_OUTCOME_TOKEN, EXPECTED_HUMANS,
    SourceGateError, original_csv, qualify_combined_source,
)

ATTESTATIONS = frozenset({
    "original_author_csv_authenticated",
    "same_biological_subject_ids_verified",
    "all_selected_omic_signals_predose_verified",
    "auc_cmax_assay_units_verified",
})
MAX_FEATURES = 2500


class BenchmarkGateError(ValueError):
    """A human/source/methodological gate, never an original cell value."""


def _fail(message: str):
    raise BenchmarkGateError(message)


def _numeric(value, *, allow_missing: bool, positive: bool) -> float | None:
    raw = str(value).strip()
    if raw.upper() in EMPTY_TOKENS:
        if allow_missing:
            return None
        _fail("A required original covariate or outcome is missing")
    try:
        number = float(raw)
    except (TypeError, ValueError, OverflowError):
        _fail("Non-numeric source covariate or outcome in selected analysis columns")
    if not isfinite(number) or (positive and number <= 0) or (not positive and number < 0):
        _fail("Non-finite or out-of-domain source covariate or outcome")
    return number


def _median_impute_log(train: Sequence[float | None], test: float | None):
    known = [v for v in train if v is not None]
    if len(known) < max(4, int(0.8 * len(train))):
        _fail("Selected omic feature lacks >=80% of original training-subject measurements")
    mid = median(known)
    return [log1p(v if v is not None else mid) for v in train], log1p(test if test is not None else mid)


def _cholesky_solve(matrix: list[list[float]], rhs: list[float]) -> list[float]:
    """Deterministic strictly positive definite ridge solve, stdlib only."""
    n = len(rhs)
    lower = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            term = matrix[i][j] - fsum(lower[i][k] * lower[j][k] for k in range(j))
            if i == j:
                if term <= 1e-10 or not isfinite(term):
                    _fail("Numerically unstable original-source ridge kernel")
                lower[i][j] = sqrt(term)
            else:
                lower[i][j] = term / lower[j][j]
    tmp = [0.0] * n
    for i in range(n):
        tmp[i] = (rhs[i] - fsum(lower[i][k] * tmp[k] for k in range(i))) / lower[i][i]
    out = [0.0] * n
    for i in range(n - 1, -1, -1):
        out[i] = (tmp[i] - fsum(lower[k][i] * out[k] for k in range(i + 1, n))) / lower[i][i]
    return out


def _predict_ridge(
    rows: list[dict[str, str]],
    held: int,
    target: str,
    clinical: Sequence[str],
    omics: Sequence[str],
) -> float:
    train = [r for i, r in enumerate(rows) if i != held]
    held_row = rows[held]
    target_log = [log(_numeric(r[target], allow_missing=False, positive=True)) for r in train]
    y_center = mean(target_log)
    columns: list[list[float]] = []
    test_values: list[float] = []

    # Fit ALL preprocessing from training persons only, including imputation.
    for col in clinical:
        values = [_numeric(r[col], allow_missing=False, positive=False) for r in train]
        held_val = _numeric(held_row[col], allow_missing=False, positive=False)
        columns.append(values)
        test_values.append(held_val)
    for col in omics:
        values = [_numeric(r[col], allow_missing=True, positive=False) for r in train]
        held_val = _numeric(held_row[col], allow_missing=True, positive=False)
        transformed, test_transformed = _median_impute_log(values, held_val)
        columns.append(transformed)
        test_values.append(test_transformed)

    if not columns:
        return exp(y_center)  # independent training-only constant baseline
    n = len(train)
    ztrain = [[0.0] * len(columns) for _ in range(n)]
    ztest = [0.0] * len(columns)
    for j, col in enumerate(columns):
        center = mean(col)
        scale = sqrt(fsum((v - center) ** 2 for v in col) / n)
        if scale < 1e-12:
            scale = 1.0
        for i in range(n):
            ztrain[i][j] = (col[i] - center) / scale
        ztest[j] = (test_values[j] - center) / scale
    # Fixed, declared alpha = number of predictors. No hyperparameter search,
    # supervised variable screen, look at held-out targets or nested-CV claim.
    penalty = float(len(columns))
    kernel = [
        [fsum(a * b for a, b in zip(ztrain[i], ztrain[j])) +
         (penalty if i == j else 0.0) for j in range(n)]
        for i in range(n)
    ]
    dual = _cholesky_solve(kernel, [v - y_center for v in target_log])
    pred_log = y_center + fsum(
        dual[i] * fsum(a * b for a, b in zip(ztrain[i], ztest))
        for i in range(n)
    )
    if not isfinite(pred_log) or abs(pred_log) > 30:
        _fail("Exploratory held-out prediction exceeds safe numeric range; no clipping")
    return exp(pred_log)


def _metrics(truth: list[float], predictions: list[float]) -> dict:
    errors = [p - y for p, y in zip(predictions, truth)]
    return {
        "mae": round(mean(abs(e) for e in errors), 6),
        "rmse": round(sqrt(mean(e * e for e in errors)), 6),
        "bias": round(mean(errors), 6),
    }


def evaluate_original_combined(
    path: str | Path,
    *,
    person_id_col: str,
    auc_col: str,
    cmax_col: str,
    omic_columns: Sequence[str],
    clinical_columns: Sequence[str] = (),
    expected_sha256: str,
    attestations: Mapping[str, bool],
    expected_subjects: int = EXPECTED_HUMANS,
    seed: int = 20261018,
    bootstrap_replicates: int = 2000,
) -> dict:
    """LODO original-subject regression; explicit source attestations, no row export.

    Exploratory fixed-penalty ridge. Input features MUST be locked before
    target inspection and sources externally checked by the researcher.
    """
    if set(attestations) != ATTESTATIONS or not all(v is True for v in attestations.values()):
        _fail("All four documented independent human source attestations are mandatory")
    if not isinstance(expected_sha256, str) or len(expected_sha256) != 64 or any(
        c not in "0123456789abcdef" for c in expected_sha256
    ):
        _fail("Pinned lowercase original CSV SHA256 must be supplied")
    if not isinstance(seed, int) or not 100 <= bootstrap_replicates <= 10000:
        _fail("Invalid fixed uncertainty resampling configuration")
    if expected_subjects != 40:
        _fail("This protocol is predeclared for 40 independent rosuvastatin volunteers only")
    names = [*clinical_columns, *omic_columns]
    if len(omic_columns) < 2 or len(names) > MAX_FEATURES or len(names) != len(set(names)):
        _fail("Require >=2 distinct explicitly selected omic features and no duplicate columns")
    if set(names) & {person_id_col, auc_col, cmax_col}:
        _fail("Patient ID or PK outcome appears in candidate predictors")
    if any(PK_OUTCOME_TOKEN.search(col) for col in names):
        _fail("Post-dose PK descriptor cannot be used as a predictor")
    if not all(isinstance(col, str) and col.strip() for col in names):
        _fail("Missing explicit predictor column name")
    try:
        structural = qualify_combined_source(
            path, person_id_col=person_id_col, auc_col=auc_col,
            cmax_col=cmax_col, expected_subjects=expected_subjects,
        )
        cols, source_rows, digest = original_csv(path)
    except SourceGateError as error:
        _fail(str(error))
    if digest != expected_sha256:
        _fail("Actual private original source SHA256 differs from approved frozen manifest")
    if any(col not in cols for col in names):
        _fail("A prespecified clinical or molecular feature is not in the original CSV")
    # No identifiers or source values leave function; same-person order is stable
    # across source CSV row shuffles, so scores cannot change with file order.
    rows = sorted(source_rows, key=lambda row: row[person_id_col])
    if len(rows) != expected_subjects:
        _fail("Biological subject cardinality changed")
    summary = {}
    for target_name, target in (("AUC", auc_col), ("Cmax", cmax_col)):
        truth = [_numeric(row[target], allow_missing=False, positive=True) for row in rows]
        baseline = [
            _predict_ridge(rows, i, target, clinical_columns, ())
            for i in range(len(rows))
        ]
        plus_omic = [
            _predict_ridge(rows, i, target, clinical_columns, omic_columns)
            for i in range(len(rows))
        ]
        diff = [abs(p - t) - abs(b - t) for p, b, t in zip(plus_omic, baseline, truth)]
        rng = Random(seed + (0 if target_name == "AUC" else 1))
        boots = sorted(mean(diff[rng.randrange(len(diff))] for _ in diff)
                       for _ in range(bootstrap_replicates))
        lo = boots[int(0.025 * bootstrap_replicates)]
        hi = boots[min(bootstrap_replicates - 1, int(0.975 * bootstrap_replicates))]
        summary[target_name] = {
            "clinical_or_constant_baseline": _metrics(truth, baseline),
            "clinical_plus_predose_omics": _metrics(truth, plus_omic),
            "delta_mae_omics_minus_baseline": round(mean(diff), 6),
            "exploratory_paired_bootstrap_95pct_delta_mae": [round(lo, 6), round(hi, 6)],
        }
    return {
        "status": "PRIVATE_SOURCE_ATTESTED_EXPLORATORY_E3_LODO_NOT_EXTERNAL_VALIDATION",
        "study": "Silveira 2021; MTBLS3128; rosuvastatin",
        "source_sha256": digest,
        "original_structural_status": structural["source_status"],
        "source_attestation_is_human_supplied_not_software_proven": True,
        "n_independent_heldout_humans": expected_subjects,
        "n_heldout_subject_folds": expected_subjects,
        "n_selected_clinical_covariates": len(clinical_columns),
        "n_selected_predose_omic_features": len(omic_columns),
        "reference_is_clinical": bool(clinical_columns),
        "method": "fixed-penalty ridge(log PK), stdlib dual; train-only log1p/median/scale; grouped LODO",
        "hyperparameter_search_performed": False,
        "input_features_prespecification_independently_certified": False,
        "bootstrap_resamples": bootstrap_replicates,
        "outcomes": summary,
        "outputs_include_original_rows_or_individual_predictions": False,
        "interpretation": "Exploratory internal n=40 validation only; no external clinical validation, NLME/PBPK or dosing claim",
    }
