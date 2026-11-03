"""Leakage-controlled research benchmark for an *observed* enzyme activity phenotype.

The goal is to compare genotype, raw matched omics, and source-curated functional
omics on the same donors. No allele activity percentages, clinical dose, or
PopPK claim. External functional curation is NOT inferred from data correlation.

This standard-library implementation intentionally limits input dimensions and
fits one fixed ridge regression in log(activity) per held-out donor fold.
"""
from __future__ import annotations

from datetime import datetime
from math import exp, isfinite, log, sqrt
from random import Random
from statistics import mean
from typing import Mapping, Sequence

from omics2pk.functional_knowledge import assess_cyp3a5_evidence


class FunctionalBenchmarkError(ValueError):
    pass


_ALLOWED = {"transcriptomics", "proteomics", "metabolomics"}
_CURATED_KINDS = {
    "transcriptomics": {"canonical_splice_junction"},
    "proteomics": {"cyp3a5_full_length"},
    "metabolomics": {"metabolite_abundance", "tracer_flux"},
}


def _require(ok: bool, why: str) -> None:
    if not ok:
        raise FunctionalBenchmarkError(why)


def _date(value: object) -> datetime:
    _require(isinstance(value, str) and bool(value.strip()), "Original absolute timestamp required")
    try:
        value = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        raise FunctionalBenchmarkError("Invalid timestamp") from None
    _require(value.tzinfo is not None and value.utcoffset() is not None,
             "Timezone-aware timestamp required")
    return value


def _positive(value: object, name: str) -> float:
    _require(not isinstance(value, bool), f"{name}: numeric observation required")
    try:
        x = float(value)
    except (ValueError, TypeError, OverflowError):
        raise FunctionalBenchmarkError(f"{name}: numeric observation required") from None
    _require(isfinite(x) and x > 0, f"{name}: positive finite observation required")
    return x


def _solve(matrix: list[list[float]], rhs: list[float]) -> list[float]:
    """Small positive-definite ridge system solved with partial pivoting."""
    n = len(rhs)
    aug = [row[:] + [rhs[i]] for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = max(range(j, n), key=lambda i: abs(aug[i][j]))
        _require(abs(aug[pivot][j]) > 1e-12, "Singular ridge system")
        aug[pivot], aug[j] = aug[j], aug[pivot]
        s = aug[j][j]
        for k in range(j, n + 1):
            aug[j][k] /= s
        for i in range(n):
            if i != j:
                a = aug[i][j]
                for k in range(j, n + 1):
                    aug[i][k] -= a * aug[j][k]
    return [row[-1] for row in aug]


def _ridge_predict(train_x: list[list[float]], train_y: list[float],
                   test_x: list[list[float]]) -> list[float]:
    """Fit standardizer, response mean and fixed alpha=1 inside train ONLY."""
    n, p = len(train_x), len(train_x[0])
    _require(n >= 8 and p <= 42 and all(len(x) == p for x in train_x + test_x),
             "Unmatched or too-large train-only regression")
    avg = [mean(row[j] for row in train_x) for j in range(p)]
    sd = [sqrt(mean((row[j] - avg[j]) ** 2 for row in train_x)) for j in range(p)]
    sd = [v if v > 1e-12 else 1.0 for v in sd]
    z = [[(row[j] - avg[j]) / sd[j] for j in range(p)] for row in train_x]
    ym = mean(train_y)
    gram = [[sum(row[j] * row[k] for row in z) + (1.0 if j == k else 0.0)
             for k in range(p)] for j in range(p)]
    cross = [sum(z[i][j] * (train_y[i] - ym) for i in range(n)) for j in range(p)]
    coef = _solve(gram, cross)
    out = [exp(max(-30.0, min(30.0, ym + sum(
        (row[j] - avg[j]) / sd[j] * coef[j] for j in range(p)))))
           for row in test_x]
    _require(all(isfinite(v) and v > 0 for v in out), "Invalid held-out activity")
    return out


def _source_curated(omic: str, item: Mapping[str, object]) -> bool:
    """Strict whitelist: no automated biological association discovery."""
    kind = item.get("target")
    _require(kind in ({"total_gene", "canonical_splice_junction", "aberrant_splice_junction"}
                      if omic == "transcriptomics" else
                      {"pan_cyp3a", "cyp3a5_full_length", "cyp3a5_unresolved_peptide"}
                      if omic == "proteomics" else
                      {"metabolite_abundance", "tracer_flux"}), "Unrecognized assay kind")
    validated = item.get("assay_validated")
    _require(isinstance(validated, bool), "Explicit assay validation flag required")
    if not validated or kind not in _CURATED_KINDS[omic]:
        return False
    reference = item.get("external_evidence")
    _require(isinstance(reference, str) and reference.strip(),
             "Externally curated annotation requires source evidence reference")
    if omic == "metabolomics":
        enzymes = item.get("candidate_enzymes")
        _require(isinstance(enzymes, list) and enzymes
                 and all(isinstance(v, str) and v.strip() for v in enzymes),
                 "Metabolic biomarker needs its candidate enzyme list")
        # A metabolite shared by CYP3A4 and CYP3A5 is NOT CYP3A5-specific.
        # Still eligible for the exploratory (uncurated) metabolomic arm.
        return enzymes == ["CYP3A5"] and item.get("validated_as_cyp3a5_specific") is True
    return True


def benchmark_cyp3a5_activity(
    *, donor_ids: Sequence[str], diplotypes: Sequence[str],
    genotype_call_statuses: Sequence[str], measured_activity: Sequence[float],
    activity_assay_at: Sequence[str], activity_unit: str,
    activity_assay_validated: bool, modalities: Mapping[str, Mapping[str, object]],
    source_kind: str, source_sha256: str, n_folds: int = 5,
    max_lookback_days: int = 30,
) -> dict[str, object]:
    """Held-donor compare PGx vs independent omic layers and curated subset.

    modalities[name]: feature_ids, values (donor x feature), collected_at
    (donor timestamps), tissue, preprocessing_scope='per_sample' or
    'external_frozen', annotations {feature_id: {target, assay_validated,
    external_evidence, ...}}. All source rows are aligned by explicit donor IDs,
    never matched automatically by row position; caller must establish original
    donor identity for each assay before passing this aligned input.

    The outcome must be specifically measured CYP3A5 enzyme activity, NOT CYP3A
    shared probe activity, protein abundance, or an already fitted PK target.
    """
    ids = tuple(donor_ids)
    n = len(ids)
    _require(30 <= n <= 400 and all(isinstance(v, str) and v.strip() for v in ids)
             and len(set(ids)) == n, "30-400 unique actual biological donors required")
    _require(all(len(v) == n for v in (diplotypes, genotype_call_statuses,
                                      measured_activity, activity_assay_at)),
             "Every cohort vector must match original donors")
    _require(source_kind in {"synthetic", "original_source_declared"},
             "Synthetic vs original measurements must be explicit")
    _require(isinstance(source_sha256, str) and len(source_sha256) == 64
             and all(v in "0123456789abcdef" for v in source_sha256),
             "Exact source SHA-256 required")
    _require(activity_assay_validated is True and isinstance(activity_unit, str)
             and bool(activity_unit.strip()),
             "Measured, validated CYP3A5-specific activity and units required")
    _require(isinstance(n_folds, int) and not isinstance(n_folds, bool)
             and 3 <= n_folds <= 10 and n - (n + n_folds - 1) // n_folds >= 8,
             "Invalid donor-heldout fold count")
    _require(isinstance(max_lookback_days, int) and 0 <= max_lookback_days <= 365,
             "Invalid prespecified pre-activity lookback")
    _require(isinstance(modalities, Mapping) and 1 <= len(modalities) <= 3
             and set(modalities) <= _ALLOWED, "One to three supported omics needed")

    y = [_positive(v, "measured CYP3A5 activity") for v in measured_activity]
    dates = [_date(v) for v in activity_assay_at]
    genotype = []
    for dip, status in zip(diplotypes, genotype_call_statuses):
        result = assess_cyp3a5_evidence(diplotype=dip, call_status=status)
        label = result["diplotype_classification"]
        _require(label is not None, "Incomplete CYP3A5 diplotype cannot be used in PGx reference")
        genotype.append([float(label == "intermediate_expresser"),
                         float(label == "expresser")])
    _require(len({tuple(v) for v in genotype}) >= 2, "No genotype contrast in cohort")

    values: dict[str, list[list[float]]] = {}
    knowledge_columns: dict[str, list[int]] = {}
    n_total = 0
    for omic, spec in sorted(modalities.items()):
        _require(isinstance(spec, Mapping), "Invalid modality record")
        features = tuple(spec.get("feature_ids", ()))
        _require(1 <= len(features) <= 12 and len(set(features)) == len(features)
                 and all(isinstance(f, str) and f.strip() for f in features),
                 "Need 1-12 original unique feature IDs per modality")
        _require(tuple(spec.get("donor_ids", ())) == ids,
                 "Each modality must independently attest the same donor ID order")
        _require(isinstance(spec.get("tissue"), str) and spec["tissue"].strip().lower()
                 in {"liver", "hepatic", "liver_microsomes"},
                 "Validated hepatic CYP3A5 activity requires matched liver tissue")
        _require(spec.get("preprocessing_scope") in {"per_sample", "external_frozen"},
                 "Cross-donor preprocessing must be externally frozen before heldout CV")
        times = spec.get("collected_at", ())
        _require(len(times) == n, "Omic specimen dates must match all donors")
        for t, assay_t in zip(times, dates):
            difference = (assay_t - _date(t)).total_seconds() / 86400
            _require(0 <= difference <= max_lookback_days,
                     "Post-outcome or stale molecular features prohibited")
        matrix = spec.get("values", ())
        _require(len(matrix) == n, "Molecular observations must match donors")
        arr: list[list[float]] = []
        for row in matrix:
            _require(len(row) == len(features), "Molecular feature orientation mismatch")
            converted: list[float] = []
            for value in row:
                _require(not isinstance(value, bool), "Omic boolean numeric value invalid")
                try:
                    x = float(value)
                except (ValueError, TypeError, OverflowError):
                    raise FunctionalBenchmarkError("Omic measurement nonnumeric") from None
                _require(isfinite(x), "Missing or infinite omic measurement")
                converted.append(x)
            arr.append(converted)
        annotations = spec.get("annotations")
        _require(isinstance(annotations, Mapping) and set(annotations) == set(features),
                 "Every actual omic feature needs a source-attested annotation")
        curated = []
        for j, f in enumerate(features):
            a = annotations[f]
            _require(isinstance(a, Mapping), "Invalid source annotation")
            if _source_curated(omic, a):
                curated.append(j)
        values[omic] = arr
        knowledge_columns[omic] = curated
        n_total += len(features)
    _require(n_total <= 32, "Too many dimensions for small-donor benchmark")

    order = sorted(range(n), key=lambda i: ids[i])
    shuffled = list(order)
    Random(20261103).shuffle(shuffled)
    groups = [shuffled[k::n_folds] for k in range(n_folds)]
    names = ["genotype_only"] + [f"genotype_plus_{name}" for name in sorted(values)]
    if len(values) > 1:
        names.append("genotype_plus_all_omics")
    has_curated = any(knowledge_columns.values())
    if has_curated:
        names.append("genotype_plus_curated_omics")
    error: dict[str, list[float]] = {name: [] for name in names}
    full_predictions: dict[str, dict[int, float]] = {name: {} for name in names}
    for test in groups:
        test_set = set(test)
        train = [i for i in order if i not in test_set]
        train_y = [log(y[i]) for i in train]
        for name in names:
            def x(i: int) -> list[float]:
                row = list(genotype[i])
                for omic in sorted(values):
                    if name == "genotype_plus_all_omics" or name == f"genotype_plus_{omic}":
                        row.extend(values[omic][i])
                    elif name == "genotype_plus_curated_omics":
                        row.extend(values[omic][i][j] for j in knowledge_columns[omic])
                return row
            predictions = _ridge_predict([x(i) for i in train], train_y,
                                         [x(i) for i in test])
            for i, p in zip(test, predictions):
                full_predictions[name][i] = p
    for name in names:
        _require(len(full_predictions[name]) == n, "Incomplete held-donor predictions")
        error[name] = [abs(y[i] - full_predictions[name][i]) for i in order]
    out = {key: round(mean(v), 8) for key, v in error.items()}
    baseline = out["genotype_only"]
    differences = {name: round(val - baseline, 8)
                   for name, val in out.items() if name != "genotype_only"}
    # Important: paired MAE deltas are descriptive ONLY: no refitted bootstrap,
    # no hypothesis test, no 90% interval claims from tiny E2 datasets.
    return {
        "status": ("SYNTHETIC_SOFTWARE_QUALIFICATION_NO_BIOLOGICAL_CLAIM"
                   if source_kind == "synthetic" else
                   "ORIGINAL_E2_HELD_DONOR_EXPLORATORY_NOT_CLINICAL_PK"),
        "endpoint": "observed_CYP3A5_specific_enzyme_activity",
        "endpoint_unit": activity_unit,
        "n_independent_donors": n,
        "n_outer_donor_heldout_folds": n_folds,
        "fixed_estimator": "Ridge alpha=1 log(activity) train-only standardization",
        "n_source_features": {k: len(v[0]) for k, v in values.items()},
        "n_knowledge_curated_features": {k: len(knowledge_columns[k]) for k in values},
        "MAE_by_arm_observed_activity_units": out,
        "delta_MAE_vs_genotype_only": differences,
        "no_arbitrary_genotype_activity_percentages": True,
        "no_post_landmark_features_per_declared_timestamps": True,
        "feature_provenance_independently_verified": False,
        "no_individual_predictions_or_source_donor_ids_exported": True,
        "uncertainty_calibrated": False,
        "in_vivo_PK_prediction_demonstrated": False,
    }
