"""Fail-closed source gate for irregular repeated tacrolimus TDM + metabolomics.

Primary target: Zhu et al. 2022, Front Pharmacol 13:1045843
(DOI 10.3389/fphar.2022.1045843), where residual plasma from repeated
clinical tacrolimus TDM samples was used for untargeted metabolomics.

This module validates source identity and chronology only. It does not fit a
model, infer missing dose events, or treat contemporaneous metabolomics as a
prospective predictor.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime
from math import isfinite
from typing import Mapping, Sequence


class RepeatedTDMSourceError(ValueError):
    pass


def _text(value: object, field: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise RepeatedTDMSourceError(f"{field} must be a non-empty string")
    return value.strip()


def _time(value: object, field: str) -> datetime:
    raw=_text(value,field)
    try:
        out=datetime.fromisoformat(raw.replace("Z","+00:00"))
    except ValueError:
        raise RepeatedTDMSourceError(f"{field} must be ISO8601") from None
    if out.tzinfo is None or out.utcoffset() is None:
        raise RepeatedTDMSourceError(f"{field} must include a time zone")
    return out


def _positive(value: object, field: str) -> float:
    if isinstance(value,bool):
        raise RepeatedTDMSourceError(f"{field} must be positive")
    try:
        out=float(value)
    except (TypeError,ValueError):
        raise RepeatedTDMSourceError(f"{field} must be positive") from None
    if not isfinite(out) or out<=0:
        raise RepeatedTDMSourceError(f"{field} must be finite and positive")
    return out


def qualify_irregular_repeated_tdm_metabolomics(
    omics: Sequence[Mapping[str,object]],
    doses: Sequence[Mapping[str,object]],
    troughs: Sequence[Mapping[str,object]],
    *,
    drug: str="tacrolimus",
    max_omics_trough_gap_hours: float=1.0,
) -> dict:
    """Validate exact same-biospecimen repeated pre-dose TDM/metabolomics links.

    Required omics fields:
      sample_id, biospecimen_id, subject_id, modality, collected_at

    Required dose fields:
      dose_id, subject_id, drug, administered_at, amount, dose_unit, route

    Required trough fields:
      observation_id, biospecimen_id, subject_id, drug,
      previous_dose_id, next_dose_id, sampled_at,
      concentration, concentration_unit,
      daily_dose, daily_dose_unit, weight_kg, postoperative_day

    The same biospecimen_id must link metabolomics and C0. Every trough must be
    strictly bounded by explicit previous/next administered tacrolimus doses.
    No fixed D7/M1/M3 visit structure is assumed: longitudinal order is derived
    from actual timestamps, suitable for irregular 3-6-sample patient series.
    """
    target=_text(drug,"drug").casefold()
    max_gap=_positive(max_omics_trough_gap_hours,"max_omics_trough_gap_hours")
    if max_gap>24:
        raise RepeatedTDMSourceError("omics/trough matching window cannot exceed 24 h")

    dose_by_id={}
    for row in doses:
        dose_id=_text(row.get("dose_id"),"dose_id")
        if dose_id in dose_by_id:
            raise RepeatedTDMSourceError("Duplicate dose_id")
        subject=_text(row.get("subject_id"),"subject_id")
        agent=_text(row.get("drug"),"drug").casefold()
        when=_time(row.get("administered_at"),"administered_at")
        amount=_positive(row.get("amount"),"amount")
        unit=_text(row.get("dose_unit"),"dose_unit")
        route=_text(row.get("route"),"route")
        if agent==target:
            dose_by_id[dose_id]=(subject,when,amount,unit,route)
    if not dose_by_id:
        raise RepeatedTDMSourceError("No target-drug dose events")

    omics_by_biospecimen={}
    sample_ids=set()
    for row in omics:
        sample_id=_text(row.get("sample_id"),"sample_id")
        biospecimen=_text(row.get("biospecimen_id"),"biospecimen_id")
        subject=_text(row.get("subject_id"),"subject_id")
        modality=_text(row.get("modality"),"modality").upper()
        when=_time(row.get("collected_at"),"collected_at")
        if sample_id in sample_ids:
            raise RepeatedTDMSourceError("Duplicate omics sample_id")
        sample_ids.add(sample_id)
        if biospecimen in omics_by_biospecimen:
            raise RepeatedTDMSourceError("More than one molecular sample per biospecimen")
        if modality!="METABOLITE":
            raise RepeatedTDMSourceError("This gate requires METABOLITE samples")
        omics_by_biospecimen[biospecimen]=(sample_id,subject,when)
    if not omics_by_biospecimen:
        raise RepeatedTDMSourceError("No molecular samples")

    seen_obs=set()
    seen_biospecimens=set()
    concentration_units=set()
    daily_dose_units=set()
    matched_by_subject=defaultdict(list)
    gaps=[]
    postoperative_days=[]
    previous_intervals=[]
    next_intervals=[]

    for row in troughs:
        agent=_text(row.get("drug"),"drug").casefold()
        if agent!=target:
            continue
        observation=_text(row.get("observation_id"),"observation_id")
        biospecimen=_text(row.get("biospecimen_id"),"biospecimen_id")
        subject=_text(row.get("subject_id"),"subject_id")
        if observation in seen_obs:
            raise RepeatedTDMSourceError("Duplicate trough observation_id")
        if biospecimen in seen_biospecimens:
            raise RepeatedTDMSourceError("Duplicate modeled trough biospecimen_id")
        seen_obs.add(observation); seen_biospecimens.add(biospecimen)

        sampled=_time(row.get("sampled_at"),"sampled_at")
        _positive(row.get("concentration"),"concentration")
        concentration_units.add(_text(row.get("concentration_unit"),"concentration_unit"))
        _positive(row.get("daily_dose"),"daily_dose")
        daily_dose_units.add(_text(row.get("daily_dose_unit"),"daily_dose_unit"))
        _positive(row.get("weight_kg"),"weight_kg")
        postop=_positive(row.get("postoperative_day"),"postoperative_day")

        prev_id=_text(row.get("previous_dose_id"),"previous_dose_id")
        next_id=_text(row.get("next_dose_id"),"next_dose_id")
        if prev_id==next_id:
            raise RepeatedTDMSourceError("Previous and next dose IDs must differ")
        if prev_id not in dose_by_id or next_id not in dose_by_id:
            raise RepeatedTDMSourceError("Trough references unknown target-drug dose")
        prev_subject,prev_time,_,_,_=dose_by_id[prev_id]
        next_subject,next_time,_,_,_=dose_by_id[next_id]
        if subject!=prev_subject or subject!=next_subject:
            raise RepeatedTDMSourceError("Trough and bounding doses disagree on subject")
        if not (prev_time<sampled<next_time):
            raise RepeatedTDMSourceError(
                "Pre-dose trough must fall strictly between previous and next dose")

        if biospecimen not in omics_by_biospecimen:
            raise RepeatedTDMSourceError(
                "Trough lacks exact same-biospecimen metabolomics link")
        _,omics_subject,omics_time=omics_by_biospecimen[biospecimen]
        if omics_subject!=subject:
            raise RepeatedTDMSourceError(
                "Metabolomics and trough biospecimen disagree on subject")
        gap=abs((sampled-omics_time).total_seconds())/3600
        if gap>max_gap:
            raise RepeatedTDMSourceError(
                "Metabolomics and trough exceed frozen matching window")

        matched_by_subject[subject].append((sampled,biospecimen,postop))
        gaps.append(gap); postoperative_days.append(postop)
        previous_intervals.append((sampled-prev_time).total_seconds()/3600)
        next_intervals.append((next_time-sampled).total_seconds()/3600)

    n_matched=sum(len(v) for v in matched_by_subject.values())
    if n_matched==0:
        raise RepeatedTDMSourceError("No linked target-drug trough observations")
    if len(concentration_units)!=1:
        raise RepeatedTDMSourceError("Trough concentration units must be consistent")
    if len(daily_dose_units)!=1:
        raise RepeatedTDMSourceError("Daily dose units must be consistent")

    repeated=0
    prospective_pairs=0
    per_subject_counts=Counter()
    for subject,events in matched_by_subject.items():
        ordered=sorted(events,key=lambda x:x[0])
        times=[x[0] for x in ordered]
        if len(set(times))!=len(times):
            raise RepeatedTDMSourceError(
                "Repeated subject has duplicate trough timestamps; chronology ambiguous")
        n=len(ordered)
        per_subject_counts[n]+=1
        if n>=2:
            repeated+=1
            prospective_pairs+=n*(n-1)//2

    return {
        "status":"IRREGULAR_REPEATED_PREDOSE_TDM_METABOLOMICS_GATE_ONLY",
        "drug":target,
        "n_unique_subjects_linked":len(matched_by_subject),
        "n_linked_same_biospecimen_c0_omics_pairs":n_matched,
        "n_subjects_with_at_least_two_linked_samples":repeated,
        "n_earlier_omics_to_later_c0_pairs_available":prospective_pairs,
        "samples_per_subject_distribution":{
            str(k):per_subject_counts[k] for k in sorted(per_subject_counts)
        },
        "max_omics_trough_gap_hours_allowed":max_gap,
        "max_observed_omics_trough_gap_hours":round(max(gaps),6),
        "min_postoperative_day":round(min(postoperative_days),6),
        "max_postoperative_day":round(max(postoperative_days),6),
        "same_biospecimen_for_omics_and_c0_required":True,
        "trough_bounded_by_exact_previous_and_next_dose_events":True,
        "daily_dose_and_weight_available_for_c0d_reconstruction":True,
        "same_sample_metabolomics_to_c0_is_contemporaneous_not_prospective":True,
        "prospective_use_requires_earlier_sample_omics_and_later_c0":True,
        "future_model_split_must_group_by_whole_subject":True,
        "fixed_nominal_visit_schedule_assumed":False,
        "individual_ids_times_doses_concentrations_exported":False,
        "not_dense_pk_or_auc":True,
        "not_popPK_validation":True,
    }
