"""Aggregate-only readiness assessment for MultiOmics → original observed PK.

This is a *source availability* gate, not proof of a causal PK predictor.
Only original individual patient/specimen keys are checked in memory; never
return names, pseudonyms, sample identifiers, PK values, or omic values.
No optional ML dependencies, and absolutely no network / data uploads.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime
from itertools import combinations
import json
from pathlib import Path

from .exchange import inspect_exchange, _csv, _file


class PKReadinessError(ValueError):
    pass


OMICS = ("transcriptomics", "proteomics", "metabolomics")
REQUIRED_FLAGS = ("has_clinical", "has_measured_pk", "has_pgx", "has_dose_timing")


def _require(value, message):
    if not value:
        raise PKReadinessError(message)


def _datetime(value):
    if not isinstance(value, str) or not value:
        return None
    try:
        dt = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return None
    return dt if dt.tzinfo is not None and dt.utcoffset() is not None else None


def assess_pk_readiness(exchange_dir, availability=None, *, lookback_days=30):
    """Audit original biological-unit source coverage without calculating PK.

    The optional caller-supplied *private* availability rows are metadata ONLY:
    subject_id, has_clinical, has_measured_pk, has_pgx, has_dose_timing
    (EXACT true/false booleans); optional assay_ids_by_omic and
    sample_ids_by_omic identity dictionaries, collected_at_by_omic of
    source-authenticated timezone-aware ISO collection instants, and
    timezone-aware index_dose_at.

    A person is temporally *qualified* for a layer only when the source
    samples.csv itself has an absolute collection_time and the caller-provided
    collection time agrees exactly; the same specimen must precede index dose.
    Private sidecar alone cannot manufacture pre-dose collection evidence.

    All counts are per UNIQUE biological person, not assay rows. Repeated
    specimens are separately counted and excluded from 1-assay eligibility.
    No independent clinical PK/genotype authenticity is inferred.
    """
    _require(isinstance(lookback_days, (float, int)) and
             not isinstance(lookback_days, bool) and
             0 < lookback_days <= 365, "Invalid prespecified predose lookback")
    root = Path(exchange_dir).resolve()
    report = inspect_exchange(root)
    mf = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    _require(mf["study"]["biological_unit"] == "participant",
             "Clinical patient-level PK needs an original participant biological unit")
    _, samples = _csv(_file(root, mf["files"]["samples"]),
                      ["subject_id", "sample_id", "assay_id", "omic"])
    patients = {r["subject_id"] for r in samples}
    assays = defaultdict(lambda: defaultdict(list))
    for r in samples:
        assays[r["omic"]][r["subject_id"]].append(r)
    layers = sorted(assays)
    _require(set(layers) <= set(OMICS), "Unknown research molecular layer")
    layer_people = {omic: {pid for pid, records in assays[omic].items()
                          if len(records) == 1} for omic in layers}
    result = {
        "status": "MULTIOMICS_ORIGINAL_PK_AVAILABILITY_AUDIT_NOT_EFFECT_VALIDATION",
        "source_validation": report["validation"],
        "n_original_biological_people": len(patients),
        "n_omic_assay_rows": len(samples),
        "n_unique_one_assay_people_by_modality": {
            omic: len(layer_people[omic]) for omic in layers},
        "n_people_with_repeated_assays_by_modality": {
            omic: sum(len(v) > 1 for v in assays[omic].values())
            for omic in layers},
        "n_intersecting_biological_people_by_modalities": {},
        "molecular_pk_reference_available": False,
        "private_original_pk_and_pgx_identity_unverified_by_reader": True,
        "external_pk_validation": False,
        "no_individual_ids_or_source_values_exported": True,
    }
    for count in range(1, len(layers) + 1):
        for subset in combinations(layers, count):
            key = "+".join(subset)
            common = set.intersection(*(layer_people[o] for o in subset))
            result["n_intersecting_biological_people_by_modalities"][key] = len(common)
    if availability is None:
        result["availability_sidecar_supplied"] = False
        result["n_people_with_original_measured_pk_declared"] = None
        result["n_clinical_pk_pgx_reference_declared"] = None
        result["n_predose_qualified_by_modalities"] = None
        result["blocking_requirements"] = [
            "same-person measured PK availability metadata",
            "authentic CYP/UGT/transporter PGx availability metadata",
            "original specimen/assay and absolute molecular-versus-dose time crosswalk",
        ]
        return result

    _require(isinstance(availability, (list, tuple)),
             "Private original clinical availability must be metadata records")
    private = {}
    for row in availability:
        _require(isinstance(row, dict) and isinstance(row.get("subject_id"), str)
                 and row["subject_id"] and row["subject_id"] not in private,
                 "Original biological people cannot be missing or repeated")
        for flag in REQUIRED_FLAGS:
            _require(type(row.get(flag)) is bool,
                     "Clinical PK/PGx availability must use explicit booleans")
        private[row["subject_id"]] = row
    _require(set(private) <= patients,
             "Private original clinical availability includes unknown participants")
    result["availability_sidecar_supplied"] = True
    result["n_private_clinical_people_declared"] = len(private)
    result["n_people_with_original_measured_pk_declared"] = sum(
        r["has_measured_pk"] for r in private.values())
    clinical_pk = {p for p, row in private.items()
                   if row["has_clinical"] and row["has_measured_pk"]
                   and row["has_dose_timing"]}
    reference = {p for p in clinical_pk if private[p]["has_pgx"]}
    result["n_clinical_pk_without_full_pgx_declared"] = len(clinical_pk - reference)
    result["n_clinical_pk_pgx_reference_declared"] = len(reference)
    result["molecular_pk_reference_available"] = bool(reference)
    joined = defaultdict(set)
    predose = defaultdict(set)
    for omic in layers:
        for pid in reference & layer_people[omic]:
            row = private[pid]
            specimen = assays[omic][pid][0]
            am = row.get("assay_ids_by_omic")
            sm = row.get("sample_ids_by_omic")
            if (not isinstance(am, dict) or not isinstance(sm, dict)
                or am.get(omic) != specimen["assay_id"]
                or sm.get(omic) != specimen["sample_id"]):
                continue
            joined[omic].add(pid)
            cm = row.get("collected_at_by_omic")
            sample_time = _datetime(specimen.get("collection_time", ""))
            declared = _datetime(cm.get(omic)) if isinstance(cm, dict) else None
            dose = _datetime(row.get("index_dose_at"))
            if sample_time is None or declared is None or dose is None:
                continue
            if sample_time != declared:
                continue
            lag = (dose - sample_time).total_seconds() / 86400
            if 0 <= lag <= lookback_days:
                predose[omic].add(pid)
    result["n_reference_and_verified_original_specimen_join_by_modality"] = {
        omic: len(joined[omic]) for omic in layers}
    result["n_reference_and_exact_predose_clock_by_modality"] = {
        omic: len(predose[omic]) for omic in layers}
    result["n_predose_qualified_by_modalities"] = {
        "+".join(subset): len(set.intersection(*(predose[o] for o in subset)))
        for count in range(1, len(layers) + 1)
        for subset in combinations(layers, count)
    }
    result["n_clinical_reference_missing_exact_original_omic_clock"] = {
        omic: len((reference & layer_people[omic]) - predose[omic])
        for omic in layers
    }
    result["lookback_days"] = lookback_days
    result["clinical_pk_and_pgx_values_not_loaded_or_independently_verified"] = True
    result["blocking_requirements"] = (
        [] if any(predose.values()) else [
            "No molecular cohort has source-and-private verified original predose timestamps",
            "Check actual original source sample_id and assay_id crosswalk",
        ])
    return result
