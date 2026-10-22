"""Fail-closed longitudinal pre-dose trough / omics identity-time gate.

Designed for studies such as He et al. 2022 (DOI 10.3390/molecules27113517)
where tacrolimus C0 is sampled before the morning dose and molecular material is
collected longitudinally. This validates source chronology only; it does not
fit PK, infer missing dose times, or treat concurrent metabolomics as a
prospective predictor.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime
from math import isfinite
from typing import Mapping, Sequence


class PredoseTroughLinkError(ValueError):
    pass


def _text(value: object, field: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise PredoseTroughLinkError(f"{field} must be a non-empty string")
    return value.strip()


def _time(value: object, field: str) -> datetime:
    raw = _text(value, field)
    try:
        out = datetime.fromisoformat(raw.replace("Z", "+00:00"))
    except ValueError:
        raise PredoseTroughLinkError(f"{field} must be ISO8601") from None
    if out.tzinfo is None or out.utcoffset() is None:
        raise PredoseTroughLinkError(f"{field} must include a time zone")
    return out


def _positive(value: object, field: str) -> float:
    if isinstance(value, bool):
        raise PredoseTroughLinkError(f"{field} must be positive")
    try:
        out = float(value)
    except (TypeError, ValueError):
        raise PredoseTroughLinkError(f"{field} must be positive") from None
    if not isfinite(out) or out <= 0:
        raise PredoseTroughLinkError(f"{field} must be finite and positive")
    return out


def qualify_predose_trough_longitudinal_link(
    omics: Sequence[Mapping[str, object]],
    doses: Sequence[Mapping[str, object]],
    troughs: Sequence[Mapping[str, object]],
    *,
    drug: str,
    visit_order: Sequence[str],
    max_omics_trough_gap_hours: float = 1.0,
) -> dict:
    """Validate same-person pre-dose C0 and longitudinal molecular sampling.

    Every C0 must explicitly reference the preceding administered dose and the
    next scheduled/administered dose. The observed trough time must be strictly
    after the preceding dose and strictly before the next dose. This prevents a
    pre-dose C0 from being incorrectly linked as a post-dose concentration.

    One same-visit molecular sample is required for each modeled trough. The
    function also counts earlier-visit molecular -> later-visit trough pairs,
    which are the only links eligible for a genuinely prospective longitudinal
    benchmark. No patient identifiers, timestamps, concentrations or doses are
    returned.
    """
    target = _text(drug, "drug").casefold()
    visits = tuple(_text(v, "visit_order") for v in visit_order)
    if len(visits) < 2 or len(set(visits)) != len(visits):
        raise PredoseTroughLinkError("visit_order must contain at least two unique visits")
    rank = {v: i for i, v in enumerate(visits)}
    max_gap = _positive(max_omics_trough_gap_hours, "max_omics_trough_gap_hours")
    if max_gap > 24:
        raise PredoseTroughLinkError("omics/trough matching window cannot exceed 24 h")

    dose_by_id = {}
    for row in doses:
        dose_id = _text(row.get("dose_id"), "dose_id")
        if dose_id in dose_by_id:
            raise PredoseTroughLinkError("Duplicate dose_id")
        subject = _text(row.get("subject_id"), "subject_id")
        agent = _text(row.get("drug"), "drug").casefold()
        when = _time(row.get("administered_at"), "administered_at")
        _positive(row.get("amount"), "amount")
        _text(row.get("dose_unit"), "dose_unit")
        _text(row.get("route"), "route")
        if agent == target:
            dose_by_id[dose_id] = (subject, when)
    if not dose_by_id:
        raise PredoseTroughLinkError("No target-drug dose events")

    omics_by_subject_visit = defaultdict(list)
    omics_ids = set()
    for row in omics:
        sample_id = _text(row.get("sample_id"), "sample_id")
        if sample_id in omics_ids:
            raise PredoseTroughLinkError("Duplicate omics sample_id")
        omics_ids.add(sample_id)
        subject = _text(row.get("subject_id"), "subject_id")
        visit = _text(row.get("visit"), "visit")
        if visit not in rank:
            raise PredoseTroughLinkError("Omics visit is outside frozen visit_order")
        modality = _text(row.get("modality"), "modality").upper()
        if modality != "METABOLITE":
            raise PredoseTroughLinkError("This source gate requires METABOLITE samples")
        when = _time(row.get("collected_at"), "collected_at")
        omics_by_subject_visit[(subject, visit)].append((sample_id, when))
    if not omics_ids:
        raise PredoseTroughLinkError("No molecular samples")

    seen_observations = set()
    seen_subject_visit = set()
    linked_visits_by_subject = defaultdict(set)
    matched = 0
    visit_counts = Counter()
    sample_gaps = []
    interval_hours = []
    concentration_units = set()

    for row in troughs:
        agent = _text(row.get("drug"), "drug").casefold()
        if agent != target:
            continue
        observation = _text(row.get("observation_id"), "observation_id")
        if observation in seen_observations:
            raise PredoseTroughLinkError("Duplicate trough observation_id")
        seen_observations.add(observation)
        subject = _text(row.get("subject_id"), "subject_id")
        visit = _text(row.get("visit"), "visit")
        if visit not in rank:
            raise PredoseTroughLinkError("Trough visit is outside frozen visit_order")
        key = (subject, visit)
        if key in seen_subject_visit:
            raise PredoseTroughLinkError("More than one modeled trough for subject/visit")
        seen_subject_visit.add(key)

        sampled = _time(row.get("sampled_at"), "sampled_at")
        _positive(row.get("concentration"), "concentration")
        concentration_units.add(_text(row.get("concentration_unit"), "concentration_unit"))

        previous_id = _text(row.get("previous_dose_id"), "previous_dose_id")
        next_id = _text(row.get("next_dose_id"), "next_dose_id")
        if previous_id == next_id:
            raise PredoseTroughLinkError("Previous and next dose IDs must differ")
        if previous_id not in dose_by_id or next_id not in dose_by_id:
            raise PredoseTroughLinkError("Trough references unknown target-drug dose")
        previous_subject, previous_time = dose_by_id[previous_id]
        next_subject, next_time = dose_by_id[next_id]
        if subject != previous_subject or subject != next_subject:
            raise PredoseTroughLinkError("Trough and bounding doses disagree on subject")
        if not (previous_time < sampled < next_time):
            raise PredoseTroughLinkError(
                "Pre-dose trough must fall strictly between previous and next dose"
            )

        candidates = omics_by_subject_visit.get(key, [])
        if len(candidates) != 1:
            raise PredoseTroughLinkError(
                "Require exactly one same-subject same-visit metabolomics sample"
            )
        _, omics_time = candidates[0]
        gap = abs((sampled - omics_time).total_seconds()) / 3600
        if gap > max_gap:
            raise PredoseTroughLinkError(
                "Metabolomics and trough samples exceed frozen matching window"
            )

        matched += 1
        visit_counts[visit] += 1
        linked_visits_by_subject[subject].add(visit)
        sample_gaps.append(gap)
        interval_hours.append((sampled - previous_time).total_seconds() / 3600)

    if matched == 0:
        raise PredoseTroughLinkError("No linked target-drug trough observations")
    if len(concentration_units) != 1:
        raise PredoseTroughLinkError("Trough concentration units must be consistent")

    prospective_pairs = 0
    subjects_with_repeated = 0
    for subject, subject_visits in linked_visits_by_subject.items():
        ordered = sorted((rank[v], v) for v in subject_visits)
        if len(ordered) >= 2:
            subjects_with_repeated += 1
        for i, _ in ordered:
            prospective_pairs += sum(j > i for j, _ in ordered)

    return {
        "status": "PREDOSE_TROUGH_LONGITUDINAL_IDENTITY_TIME_GATE_ONLY",
        "drug": target,
        "frozen_visit_order": list(visits),
        "n_unique_subjects_linked": len(linked_visits_by_subject),
        "n_linked_same_visit_trough_omics_pairs": matched,
        "linked_visit_counts": {v: visit_counts[v] for v in visits},
        "n_subjects_with_at_least_two_linked_visits": subjects_with_repeated,
        "n_earlier_omics_to_later_trough_pairs_available": prospective_pairs,
        "max_omics_trough_gap_hours_allowed": max_gap,
        "max_observed_omics_trough_gap_hours": round(max(sample_gaps), 6),
        "median_or_summary_dose_not_accepted_as_event_history": True,
        "trough_bounded_by_exact_previous_and_next_dose_events": True,
        "same_visit_molecular_state_is_concurrent_exposure_sensing_not_prospective": True,
        "prospective_use_requires_earlier_visit_omics_and_later_visit_trough": True,
        "future_model_split_must_group_by_whole_subject": True,
        "concentration_unit_consistent": True,
        "individual_ids_times_doses_concentrations_exported": False,
        "not_dense_PK_or_AUC": True,
        "not_PopPK_validation": True,
    }
