"""Research-only temporal assay→clinical PK or assay→PD gate for Omics2PK.

The actual upstream MultiOmics validated SHA-pinned DOWNSTREAM_EXCHANGE_V1
remains the molecular contract. This module controls downstream clinical
event linking using genuine biological subject IDs and ORIGINAL timestamps.
Does NOT load or emit a patient, molecular abundance or clinical response.

Do NOT interpret true drug exposure or a post-treatment omics sample as
prospective predictors of that SAME drug's earlier exposure; no circular
'drug measured by mass spec predicts drug concentration' use.
"""
from __future__ import annotations

from datetime import datetime
from math import isfinite

from .pk_link_contract import PKLinkError,qualify_clinical_link


class OmicsPDLinkError(ValueError):
    pass


def _label(v):
    if not isinstance(v,str) or not v.strip() or len(v)>200:
        raise OmicsPDLinkError("Missing original source key/label")
    return v.strip()


def _date(v):
    text=_label(v)
    try: dt=datetime.fromisoformat(text.replace("Z","+00:00"))
    except ValueError:raise OmicsPDLinkError("Bad absolute ISO timestamp") from None
    if dt.tzinfo is None or dt.utcoffset() is None:
        raise OmicsPDLinkError("Timezone required for clinical and omics events")
    return dt


def qualify_temporal_pd_link(
    original_molecular_samples,
    original_treatment_events,
    original_response_events,
    *,
    endpoint:str,
    min_followup_days:float=1.,
    max_molecular_lookback_days:float=30.,
    required_modalities=("metabolomics",),
):
    """Validate pre-treatment molecular assays and measured POST-treatment PD.

    For observational/experimental studies the endpoint may be death,
    rejection or continuous clinical response, but this gate does NOT imply
    randomization or causally individualized treatment benefit.
    """
    target=_label(endpoint)
    if not (isfinite(float(min_followup_days)) and min_followup_days>0
            and isfinite(float(max_molecular_lookback_days))
            and 0<max_molecular_lookback_days<=365):
        raise OmicsPDLinkError("Unreasonable research PD followup/prospectivity window")
    mods=tuple(_label(m).lower() for m in required_modalities)
    if not mods or len(set(mods))!=len(mods) or not set(mods)<={"transcriptomics","proteomics","metabolomics","genomics"}:
        raise OmicsPDLinkError("Invalid original modality selection")
    assays={}
    assay_ids=set()
    for row in original_molecular_samples:
        sample=_label(row.get("sample_id"))
        if sample in assay_ids:raise OmicsPDLinkError("Repeated source molecular assay ID")
        assay_ids.add(sample)
        person=_label(row.get("subject_id"))
        mode=_label(row.get("omic")).lower()
        when=_date(row.get("collected_at"))
        tissue=_label(row.get("tissue"))
        assays.setdefault((person,mode),[]).append((when,tissue))
    arms={}
    intervention_ids=set()
    for row in original_treatment_events:
        person=_label(row.get("subject_id"))
        if person in arms:raise OmicsPDLinkError("Same biological subject has ambiguous intervention start")
        intervention=_label(row.get("intervention"))
        start=_date(row.get("started_at"))
        allocation=_label(row.get("allocation")).lower()
        if allocation not in {"randomized","observational","not_recorded"}:
            raise OmicsPDLinkError("Do not assume trial randomization")
        arms[person]=(intervention,start,allocation)
        intervention_ids.add(person)
    if not arms:
        raise OmicsPDLinkError("No authentic treatment start linked to original subjects")
    matched=set()
    response_ids=set()
    observed=0
    for row in original_response_events:
        key=_label(row.get("outcome_event_id"))
        if key in response_ids:raise OmicsPDLinkError("Repeated original endpoint observation ID")
        response_ids.add(key)
        person=_label(row.get("subject_id"))
        if _label(row.get("endpoint"))!=target:continue
        if person not in arms:
            raise OmicsPDLinkError("Endpoint has no same-person intervention record")
        intervention,start,_=arms[person]
        observed_date=_date(row.get("assessed_at"))
        days=(observed_date-start).total_seconds()/86400.
        if days<min_followup_days:
            raise OmicsPDLinkError("Endpoint observed before prescribed prospective followup")
        _label(row.get("observed_value"))
        for m in mods:
            candidates=[when for when,tissue in assays.get((person,m),[])
                        if when<=start and (start-when).total_seconds()/86400 <= max_molecular_lookback_days]
            if len(candidates)!=1:
                raise OmicsPDLinkError("No unambiguous prospective same-person predose molecular layer")
        observed+=1
        matched.add(person)
    if not observed:raise OmicsPDLinkError("No measured prospective original PD outcome")
    return {
        "status":"ORIGINAL_SAME_PERSON_PROSPECTIVE_MOLECULAR_TO_PD_LINK_ONLY",
        "target_PD_endpoint":target,
        "n_original_intervention_subjects":len(arms),
        "n_original_molecular_assay_records":len(assay_ids),
        "n_original_prospective_endpoint_observations":observed,
        "n_original_unique_subjects_with_PD_and_omics":len(matched),
        "prespecified_pre_treatment_omic_modalities":list(mods),
        "observational_association_not_causal_individual_PD_treatment_effect":True,
        "original_person_identifiers_dose_outcomes_never_emitted":True,
    }


def qualify_omics2pk_or_pd_staging(
    *,
    direction,
    original_omics_samples,
    original_doses_or_treatments,
    original_observations,
    settings,
):
    if direction=="PK":
        # Source notes from PKLink existing original/source timing contract.
        return qualify_clinical_link(
            original_omics_samples,original_doses_or_treatments,
            original_observations,**settings)
    if direction=="PD":
        return qualify_temporal_pd_link(
            original_omics_samples,original_doses_or_treatments,
            original_observations,**settings)
    raise OmicsPDLinkError("Declared downstream must be PK or PD")
