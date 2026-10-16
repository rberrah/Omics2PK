"""Audit the actual MultiOmics DOWNSTREAM_EXCHANGE_V1 for E2→E3→E4 PK.

This is a read-only, source-hash-verified compatibility check. It emits ONLY
aggregate QC, not sample IDs, dates, raw values, genomic features or patients.
PK outcomes and covariate scores are deliberately not invented or imputed.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import json
from pathlib import Path

from .exchange import ExchangeError, _csv, _file, inspect_exchange


_ALLOWED_UNITS = {
    "h": 1.0, "hour": 1.0, "hours": 1.0,
    "min": 1/60, "minute": 1/60, "minutes": 1/60,
    "d": 24.0, "day": 24.0, "days": 24.0,
}
_MODALITIES = ("transcriptomics", "proteomics", "metabolomics")


def _temporal_hour(raw: str, unit: str) -> float | None:
    """V1 numeric relative time only; categorical T0/T12 is NOT a unit."""
    from math import isfinite
    if not raw or not unit or unit.strip().lower() not in _ALLOWED_UNITS:
        return None
    try:
        v = float(raw)
    except (ValueError, TypeError, OverflowError):
        return None
    if not isfinite(v) or v < 0:
        return None
    return v * _ALLOWED_UNITS[unit.strip().lower()]


def inspect_pk_readiness(directory: str | Path) -> dict:
    """Read native exchange, then inspect longitudinal source metadata.

    A "PASS" for source file hashes is NOT molecular-to-PK modeling readiness.
    An explicit shared time origin is required later to join clinical doses
    and molecular samples: matching numeric hours alone does not prove it.
    """
    path = Path(directory).resolve()
    base = inspect_exchange(path)  # source SHA-256 and matrix identity QC
    manifest = json.loads((path / "manifest.json").read_text(encoding="utf-8"))
    headers, sample_rows = _csv(_file(path, manifest["files"]["samples"]),
                               ["assay_id", "sample_id", "subject_id", "omic"])
    has_collection_time = "collection_time" in headers
    has_time_unit = "collection_time_unit" in headers
    has_source_timepoint = "timepoint" in headers
    has_tissue = "tissue" in headers

    sample_to_owner: dict[str, str] = {}
    sample_to_tissue: dict[str, str] = {}
    sample_to_time: dict[str, float] = {}
    assay_counts_by_sample = Counter()
    unique_sample_by_omic_subject = defaultdict(set)
    n_invalid_time = n_missing_time = 0
    by_subject_visits = defaultdict(set)
    distinct_subjects = {r["subject_id"] for r in sample_rows}
    for r in sample_rows:
        sample = r["sample_id"]
        subject = r["subject_id"]
        assay_counts_by_sample[(r["omic"],sample)] += 1
        if sample in sample_to_owner and sample_to_owner[sample] != subject:
            raise ExchangeError("One specimen belongs to multiple subjects")
        sample_to_owner[sample] = subject
        if has_tissue:
            tissue = r.get("tissue", "").strip()
            if sample in sample_to_tissue and sample_to_tissue[sample] and tissue and sample_to_tissue[sample] != tissue:
                raise ExchangeError("Contradictory tissues for a shared specimen")
            if tissue:
                sample_to_tissue[sample] = tissue
        unique_sample_by_omic_subject[(r["omic"], subject)].add(sample)
        raw = r.get("collection_time", "").strip() if has_collection_time else ""
        unit = r.get("collection_time_unit", "").strip() if has_time_unit else ""
        if not raw:
            n_missing_time += 1
            continue
        hour = _temporal_hour(raw, unit)
        if hour is None:
            n_invalid_time += 1
            continue
        if sample in sample_to_time and abs(sample_to_time[sample] - hour) > 1e-8:
            raise ExchangeError("Shared specimen has conflicting quantitative collection times")
        sample_to_time[sample] = hour
        by_subject_visits[subject].add(hour)
    temporal_eligible = n_missing_time == 0 and n_invalid_time == 0 and (
        has_collection_time and has_time_unit
    )
    technical_replicate_groups = sum(n - 1 for n in assay_counts_by_sample.values() if n > 1)
    subjects_with_repeat_assays = len({
        subject for (_,subject), specimens in unique_sample_by_omic_subject.items()
        if len(specimens) > 1
    })
    longitudinal_subjects = sum(len(visits) > 1 for visits in by_subject_visits.values())

    needs = []
    if not temporal_eligible:
        needs.append("NUMERIC_COLLECTION_TIME_AND_EXPLICIT_UNIT_PER_ASSAY")
    if not has_source_timepoint and not temporal_eligible:
        needs.append("PRESERVE_CATEGORICAL_TIMEPOINT_FROM_ORIGINAL_MULTIOMICS_SOURCE")
    if not any(r.get("tissue", "").strip() for r in sample_rows):
        needs.append("EXPLICIT_TISSUE_AND_SPECIMEN_MATRIX_PROVENANCE")
    if manifest["software"].get("commit", "unknown") in ("unknown", "", None):
        needs.append("REAL_PRODUCER_COMMIT_AND_RUNTIME_PROVENANCE")
    if any(x["value_type"] in ("source_values_unspecified", "raw_counts")
           for x in manifest.get("matrices", [])):
        needs.append("DOCUMENT_PER_LAYER_FIT_READY_NORMALIZATION_NOT_RAW_COUNTS")
    # These are NOT present in the V1 export; explicitly require separate
    # study-specific PK/clinical inputs or a strictly versioned sidecar.
    needs.extend((
        "SUBJECT_MATCHED_ORIGINAL_PK_DOSE_CONCENTRATION_EVENTS",
        "REPRODUCIBLE_NCA_ENDPOINTS_WITH_TIMEWINDOW_AND_PK_UNITS",
        "INDIVIDUAL_PREDECLARED_BIOLOGICAL_SCORES_WITH_ASSAY_ID_AND_FOLD_PROVENANCE",
        "DECLARED_SHARED_MOLECULAR_PK_TIME_ORIGIN",
    ))
    return {
        "status": "E2_SOURCE_EXCHANGE_QC_ONLY_NOT_PK_READY" if needs else "PK_READINESS_CONDITIONAL",
        "source_exchange_schema_version": base["schema_version"],
        "n_independent_biological_subjects": len(distinct_subjects),
        "n_original_omic_assays": len(sample_rows),
        "n_original_unique_specimens": len(sample_to_owner),
        "n_subjects_with_repeated_same_modality_specimens": subjects_with_repeat_assays,
        "n_subjects_with_two_or_more_distinct_numeric_collection_times": longitudinal_subjects,
        "n_extra_technical_repeat_assays_not_independent_biology": technical_replicate_groups,
        "n_omic_assays_missing_collection_time": n_missing_time,
        "n_omic_assays_with_nonconvertible_time_or_unit": n_invalid_time,
        "temporal_numeric_metadata_ready": temporal_eligible,
        "has_categorical_timepoint_column": has_source_timepoint,
        "source_biological_unit": manifest["study"].get("biological_unit"),
        "source_validation_status": base["validation"],
        "matrix_value_types": {m["omic"]: m["value_type"] for m in manifest["matrices"]},
        "has_structured_person_level_NCA_targets": False,
        "has_original_PK_dose_concentration_history": False,
        "has_prospective_subject_assay_level_model_scores": False,
        "missing_for_E3_and_E4": needs,
        "clinical_or_PK_effect_validated": False,
        "no_patient_or_assay_identifiers_exported": True,
        "note": "Collection_time is only meaningful with an explicit shared PK time origin. Group contrast effect/q values cannot substitute individual fold-frozen biological covariates.",
    }
