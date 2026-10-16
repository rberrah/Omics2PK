"""Fail-closed, local-only clinical Omics -> PK study linkage metadata contract.

Input rows contain real dates/subject identifiers and must NEVER be logged,
uploaded to GitHub Actions, or exported by the QC summary. This does NOT
estimate a PopPK model, adjust for clinical/PGx confounders, or imply causality.
"""
from __future__ import annotations

from collections import Counter
from datetime import datetime
from math import isfinite
from typing import Mapping, Sequence


class PKLinkError(ValueError):
    """Unverifiable identity, dose/concentration linkage or temporal alignment."""


def _s(value: object, field: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise PKLinkError(f"{field} must be a non-empty string")
    return value.strip()


def _dt(value: object, field: str) -> datetime:
    raw = _s(value, field)
    try:
        timestamp = datetime.fromisoformat(raw.replace("Z", "+00:00"))
    except ValueError:
        raise PKLinkError(f"{field} is not an ISO8601 timestamp") from None
    if timestamp.tzinfo is None or timestamp.utcoffset() is None:
        raise PKLinkError(f"{field} must explicitly contain a time zone")
    return timestamp


def _pos(value: object, field: str) -> float:
    if isinstance(value, bool):
        raise PKLinkError(f"{field} must be a real positive number")
    try:
        num = float(value)
    except (TypeError, ValueError):
        raise PKLinkError(f"{field} must be a real positive number") from None
    if not isfinite(num) or num <= 0.0:
        raise PKLinkError(f"{field} must be finite and strictly positive")
    return num


def qualify_clinical_link(
    omics: Sequence[Mapping[str, object]],
    doses: Sequence[Mapping[str, object]],
    concentrations: Sequence[Mapping[str, object]],
    *,
    drug: str,
    max_omics_lag_hours: float,
    selected_modality: str,
    selection_rule: str = "most_recent_before_dose",
) -> dict:
    """Validate identifiers and a BEFORE-DOSE molecular covariate for each PK sample.

    A concentration must link to an explicitly identified dose, never simply
    "the nearest dose". This is intentionally strict: sparse longitudinal
    concentration profiles and actual dosing histories can require a richer
    PK event contract after an independent protocol review.
    """
    target = _s(drug, "drug").casefold()
    if selected_modality not in {"RNA", "PROTEIN", "METABOLITE", "MULTIOMICS"}:
        raise PKLinkError("Modality must be explicitly RNA, PROTEIN, METABOLITE or MULTIOMICS")
    if selection_rule != "most_recent_before_dose":
        raise PKLinkError("Unregistered temporal selection rule")
    max_lag = _pos(max_omics_lag_hours, "max_omics_lag_hours")
    if max_lag > 24 * 365:
        raise PKLinkError("Unreasonably large molecular-state lookback window")

    samples = {}
    sample_groups = {}
    for row in omics:
        sid = _s(row.get("sample_id"), "sample_id")
        person = _s(row.get("subject_id"), "subject_id")
        when = _dt(row.get("collected_at"), "collected_at")
        modality = _s(row.get("modality"), "modality").upper()
        if modality not in {"RNA", "PROTEIN", "METABOLITE"}:
            raise PKLinkError("Unsupported individual assay type")
        if sid in samples:
            raise PKLinkError("Duplicate omics sample ID")
        # No modality conversion: one RNA measurement is not multiomics.
        samples[sid] = (person, when, modality)
        sample_groups.setdefault((person, modality), []).append((sid, when))
    if not samples:
        raise PKLinkError("No actual omics sample metadata supplied")

    drug_doses = {}
    dose_subjects = set()
    for row in doses:
        dose_id = _s(row.get("dose_id"), "dose_id")
        person = _s(row.get("subject_id"), "subject_id")
        administered = _dt(row.get("administered_at"), "administered_at")
        agent = _s(row.get("drug"), "drug").casefold()
        _pos(row.get("amount"), "amount")
        _s(row.get("dose_unit"), "dose_unit")
        _s(row.get("route"), "route")
        if dose_id in drug_doses:
            raise PKLinkError("Duplicate dose event ID")
        if agent == target:
            drug_doses[dose_id] = (person, administered)
            dose_subjects.add(person)
    if not drug_doses:
        raise PKLinkError("No observed target-drug dose records")

    seen = set()
    matched = 0
    intervals = []
    lag_hours = []
    for row in concentrations:
        sample_id = _s(row.get("observation_id"), "observation_id")
        if sample_id in seen:
            raise PKLinkError("Duplicate concentration observation ID")
        seen.add(sample_id)
        person = _s(row.get("subject_id"), "subject_id")
        agent = _s(row.get("drug"), "drug").casefold()
        when = _dt(row.get("sampled_at"), "sampled_at")
        _pos(row.get("concentration"), "concentration")
        _s(row.get("concentration_unit"), "concentration_unit")
        if agent != target:
            continue
        dose_id = _s(row.get("dose_id"), "dose_id")
        if dose_id not in drug_doses:
            raise PKLinkError("Concentration references unknown/not-target dose ID")
        dose_person, dose_time = drug_doses[dose_id]
        if person != dose_person:
            raise PKLinkError("Concentration and dose subject IDs disagree")
        if when <= dose_time:
            raise PKLinkError("Concentration was not collected after its dose")
        if selected_modality == "MULTIOMICS":
            required = {"RNA", "PROTEIN", "METABOLITE"}
        else:
            required = {selected_modality}
        for assay in required:
            earlier = [(sid, t) for sid, t in sample_groups.get((person, assay), [])
                       if t <= dose_time and 0 <= (dose_time - t).total_seconds() / 3600 <= max_lag]
            if not earlier:
                raise PKLinkError("No pre-dose omics assay within frozen time window")
            # Tie at the most recent sampling timestamp cannot be disambiguated
            # without explicit biological replicate/preprocessing provenance.
            most_recent = max(t for _, t in earlier)
            if sum(t == most_recent for _, t in earlier) != 1:
                raise PKLinkError("Ambiguous most-recent omics assay timestamp")
            lag_hours.append((dose_time - most_recent).total_seconds() / 3600)
        intervals.append((when - dose_time).total_seconds() / 3600)
        matched += 1
    if matched == 0:
        raise PKLinkError("No measured concentration linked to this drug")
    if not dose_subjects:
        raise PKLinkError("Target-drug subjects cannot be established")
    # This summary contains NO patient or assay IDs, dates, concentrations,
    # gene names or sample-level values; it is not a qualified model result.
    return {
        "status": "TEMPORAL_IDENTITY_GATE_ONLY_NOT_POPPK_VALIDATION",
        "drug": target,
        "modality": selected_modality,
        "selection_rule": selection_rule,
        "max_omics_lag_hours": max_lag,
        "n_omics_samples": len(samples),
        "n_dose_events_target": len(drug_doses),
        "n_dose_subjects_target": len(dose_subjects),
        "n_target_concentrations_linked": matched,
        "n_matched_assay_links": len(lag_hours),
        "max_actual_assay_lag_hours": round(max(lag_hours), 5),
        "n_independent_subjects_at_most": len(dose_subjects),
        "no_patient_identifiers_exported": True,
        "disclaimer": "Input doses may be incomplete; additional doses/infusions and rich PK records require explicit protocol qualification before fitting.",
    }
