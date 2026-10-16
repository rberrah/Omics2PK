"""Observed and simulated longitudinal dose/concentration event contract (PK_LONG_V1).

A PK concentration depends on the WHOLE preceding dose history, not simply an
arbitrary nearest dose. This independent data/temporal qualification layer
does NOT fit an NLME model or authorize clinical dosing.

All subject-level records stay in the caller's trusted environment; only
aggregate counts are returned. Use study-relative elapsed hours or an
externally vetted clock-time conversion. Never confuse EVID with a phenotype.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from math import isfinite
from typing import Mapping, Sequence


class LongitudinalPKError(ValueError):
    """The original timed PK course cannot be interpreted unambiguously."""


def _text(v: object, field: str) -> str:
    if not isinstance(v, str) or not v.strip():
        raise LongitudinalPKError(f"{field} must be a non-empty string")
    return v.strip()


def _number(v: object, field: str, *, positive: bool = False) -> float:
    if isinstance(v, bool):
        raise LongitudinalPKError(f"{field} cannot be boolean")
    try:
        result = float(v)
    except (TypeError, ValueError, OverflowError):
        raise LongitudinalPKError(f"{field} must be a finite numeric value") from None
    if not isfinite(result) or (positive and result <= 0) or (not positive and result < 0):
        raise LongitudinalPKError(f"{field} must be finite and {'positive' if positive else 'nonnegative'}")
    return result


def _kind(row: Mapping) -> tuple[str, str, str, float]:
    subject = _text(row.get("subject_id"), "subject_id")
    kind = _text(row.get("event_type"), "event_type").lower()
    if kind not in ("dose", "observation"):
        raise LongitudinalPKError("Only dose and observation events accepted")
    source = _text(row.get("data_origin"), "data_origin").lower()
    if source not in ("observed", "simulated"):
        raise LongitudinalPKError("Every PK record needs observed/simulated provenance")
    time = _number(row.get("time_h"), "time_h")
    return subject, kind, source, time


def qualify_longitudinal_pk(
    events: Sequence[Mapping[str, object]],
    *,
    drug: str,
    time_origin: str,
    dosing_history_status: str,
    dataset_kind: str,
    require_observed_concentrations: bool = True,
) -> dict:
    """Strict data-only checks; zero BLQ concentration is NEVER invented as zero DV.

    Event values:
      dose: subject_id, time_h, event_type="dose", data_origin,
            amount_mg>0, route, formulation, dose_id
      observation: ... event_type="observation", concentration_mg_l>=0 OR
            is_blq=True and quantifiable_value absent, analyte, specimen,
            sample_id, sampling_phase ("predose" for time tied with dose).
    Required provenance: dataset_kind=single_dose|multiple_dose,
    dosing_history_status=complete|partial|unknown.

    Unknown/partial dosing histories may be loaded for research data QC,
    but must NOT be passed to mechanistic estimation as if complete.
    """
    _text(drug, "drug")
    _text(time_origin, "time_origin")
    if dataset_kind not in {"single_dose", "multiple_dose"}:
        raise LongitudinalPKError("Dataset kind must explicitly be single/multiple dose")
    if dosing_history_status not in {"complete", "partial", "unknown"}:
        raise LongitudinalPKError("Dosing history needs explicit complete/partial/unknown flag")
    if not events:
        raise LongitudinalPKError("No original longitudinal PK events")

    subjects = set()
    times_by_subject = defaultdict(list)
    doses_by_subject = defaultdict(list)
    obs_by_subject = defaultdict(list)
    source_counts = Counter()
    phase_counts = Counter()
    dose_ids = set()
    sample_ids = set()
    analytes = set()
    specimens = set()
    seen_physical_samples = set()
    n_blq = 0
    for record in events:
        if not isinstance(record, Mapping):
            raise LongitudinalPKError("Expected explicit event mappings")
        person, kind, source, time = _kind(record)
        agent = _text(record.get("drug"), "drug")
        if agent.casefold() != drug.casefold():
            raise LongitudinalPKError("Multiple drugs require separate registered PK data streams")
        subjects.add(person)
        source_counts[(kind, source)] += 1
        times_by_subject[person].append(time)
        if kind == "dose":
            did = _text(record.get("dose_id"), "dose_id")
            if did in dose_ids:
                raise LongitudinalPKError("Duplicate dose event identifier")
            dose_ids.add(did)
            _number(record.get("amount_mg"), "amount_mg", positive=True)
            _text(record.get("route"), "route")
            _text(record.get("formulation"), "formulation")
            if "concentration_mg_l" in record or "is_blq" in record:
                raise LongitudinalPKError("PK dose events may not contain a DV")
            doses_by_subject[person].append((time, did, source))
        else:
            sid = _text(record.get("sample_id"), "sample_id")
            if sid in sample_ids:
                raise LongitudinalPKError("Duplicate concentration sample identifier")
            sample_ids.add(sid)
            analyte = _text(record.get("analyte"), "analyte").casefold()
            specimen = _text(record.get("specimen"), "specimen").casefold()
            analytes.add(analyte)
            specimens.add(specimen)
            if (person, time, analyte) in seen_physical_samples:
                raise LongitudinalPKError("Two observations at same patient/time/analyte need replicate contract")
            seen_physical_samples.add((person, time, analyte))
            phase = _text(record.get("sampling_phase", "postdose"), "sampling_phase")
            if phase not in {"predose", "postdose"}:
                raise LongitudinalPKError("Invalid simultaneous-dose sampling phase")
            phase_counts[phase] += 1
            blq = record.get("is_blq", False)
            if type(blq) is not bool:
                raise LongitudinalPKError("is_blq must be boolean")
            if blq:
                if record.get("concentration_mg_l") is not None:
                    raise LongitudinalPKError("BLQ must not be replaced with zero or a fake measured concentration")
                _number(record.get("lloq_mg_l"), "lloq_mg_l", positive=True)
                n_blq += 1
            else:
                _number(record.get("concentration_mg_l"), "concentration_mg_l")
            obs_by_subject[person].append((time, phase, source))

    n_observations = sum(map(len, obs_by_subject.values()))
    n_measured = source_counts[("observation", "observed")]
    if require_observed_concentrations and not n_measured:
        raise LongitudinalPKError("No genuinely measured original concentrations")
    if not all(doses_by_subject.get(p) for p in obs_by_subject):
        raise LongitudinalPKError("Concentration records without any recorded dose for same patient")
    if dataset_kind == "single_dose" and any(len(doses_by_subject[p]) != 1 for p in subjects):
        raise LongitudinalPKError("Single-dose source cannot have several recorded dose events")
    if dataset_kind == "multiple_dose" and not any(len(doses_by_subject[p]) >= 2 for p in subjects):
        raise LongitudinalPKError("Multiple-dose classification unsupported by dose event history")

    for p, samples in obs_by_subject.items():
        dose_times = [t for t, _, _ in doses_by_subject[p]]
        first = min(dose_times)
        for t, phase, _ in samples:
            if t < first:
                raise LongitudinalPKError("Pre-first-dose baseline not declared as a separate protocol")
            if t == first and phase != "predose":
                raise LongitudinalPKError("Time-zero concentration must be explicitly marked predose")
            if phase == "predose" and t not in dose_times:
                raise LongitudinalPKError("Predose concentration must coincide with a registered dose")
            if t in dose_times and phase == "postdose":
                raise LongitudinalPKError("Same-time dose/concentration ordering is ambiguous without predose mark")

    repeated_subjects = sum(len(obs_by_subject[p]) > 1 for p in subjects)
    return {
        "status": "QUALIFIED_LONGITUDINAL_PK_EVENT_SCHEMA_ONLY_NOT_A_POPPK_FIT",
        "dataset_kind": dataset_kind,
        "source_time_origin": time_origin,
        "dose_history_status": dosing_history_status,
        "dosing_history_suitable_for_mechanistic_fit": dosing_history_status == "complete",
        "n_independent_subjects": len(subjects),
        "n_subjects_with_repeated_concentration_observations": repeated_subjects,
        "n_dose_events": sum(len(v) for v in doses_by_subject.values()),
        "n_concentration_observations": n_observations,
        "n_observed_concentrations": n_measured,
        "n_simulated_concentrations": source_counts[("observation", "simulated")],
        "n_observed_dose_events": source_counts[("dose", "observed")],
        "n_simulated_dose_events": source_counts[("dose", "simulated")],
        "n_below_quantification_limit": n_blq,
        "n_predose_observations": phase_counts["predose"],
        "n_distinct_analytes": len(analytes),
        "n_distinct_specimen_types": len(specimens),
        "no_subject_ids_or_PK_values_exported": True,
        "model_validation": False,
        "warning": "Mixed simulated/observed records are enumerated separately. PK depends on all earlier doses, not the nearest dose; exposure/NLME requires regimen and censoring qualification.",
    }


def align_molecular_states_to_landmarks(
    omics: Sequence[Mapping[str, object]],
    landmarks: Sequence[Mapping[str, object]],
    *,
    modality: str,
    max_lookback_h: float,
) -> dict:
    """Temporally qualify a longitudinal molecular predictor without peeking ahead.

    A repeated assay collected after the first dose MAY be used for a later
    prediction landmark if causally available by then, unlike baseline-only
    heuristics. Exactly same-time omics requires a separately frozen sampling
    order and currently ABSTAINS.
    """
    if modality not in {"RNA", "PROTEIN", "METABOLITE"}:
        raise LongitudinalPKError("Individual modality must be prespecified")
    lookback = _number(max_lookback_h, "max_lookback_h", positive=True)
    by_subject = defaultdict(list)
    seen = set()
    for m in omics:
        sid = _text(m.get("sample_id"), "sample_id")
        if sid in seen:
            raise LongitudinalPKError("Duplicate molecular specimen/sample ID")
        seen.add(sid)
        subject = _text(m.get("subject_id"), "subject_id")
        t = _number(m.get("time_h"), "omics time_h")
        assay = _text(m.get("modality"), "modality")
        if assay == modality:
            by_subject[subject].append(t)

    if not landmarks:
        raise LongitudinalPKError("At least one prediction landmark is required")
    n_aligned = 0
    n_ambiguous = 0
    n_missing = 0
    n_future_discarded = 0
    for landmark in landmarks:
        person = _text(landmark.get("subject_id"), "landmark subject_id")
        t = _number(landmark.get("time_h"), "landmark time_h")
        times = by_subject.get(person, ())
        n_future_discarded += sum(assay_t > t for assay_t in times)
        candidates = [k for k in times if 0 <= t - k <= lookback and k < t]
        if not candidates:
            n_missing += 1
            continue
        last = max(candidates)
        if candidates.count(last) > 1 or any(k == t for k in times):
            n_ambiguous += 1
            continue
        n_aligned += 1
    return {
        "status": "LANDMARK_QC_ONLY_NO_COVARIATE_VALUES_EXPORTED",
        "modality": modality,
        "n_prediction_landmarks": len(landmarks),
        "n_landmarks_with_preceding_molecular_measurement": n_aligned,
        "n_landmarks_missing_prior_omic": n_missing,
        "n_landmarks_with_ambiguous_simultaneous_assay": n_ambiguous,
        "n_future_omic_samples_excluded_by_landmark": n_future_discarded,
        "no_future_molecular_state_used": True,
        "note": "A prior molecular state may be time-varying. Missing or simultaneous samples must be handled by frozen protocol; patient-level train/test splits belong outside this QC output.",
    }
