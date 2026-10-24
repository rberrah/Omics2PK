"""Privacy-preserving aggregate feasibility gate for BIOMARGIN -> ABIS/TDM.

This module deliberately accepts aggregate counts and governance attestations
only. It is intended for a first data-custodian feasibility exchange before any
patient-level BIOMARGIN molecular data or ABIS therapeutic-drug-monitoring data
are moved into an Omics2PK analysis environment.

Public evidence establishes that:
- BIOMARGIN collected molecular samples in four European kidney-transplant
  centres, including CHU Limoges;
- ABIS is a real clinical therapeutic-drug-monitoring / Bayesian dose-adaptation
  database used by transplant centres;
- REXETRIS links ABIS with CRISTAL and SNDS under authorized pseudonymized data
  processing.

None of that publicly proves a BIOMARGIN <-> ABIS same-patient crosswalk.
"""
from __future__ import annotations

from typing import Mapping


class BiomarginAbisFeasibilityError(ValueError):
    pass


COUNT_FIELDS = (
    "biomargin_molecular_patients",
    "exact_biomargin_abis_crosswalk_patients",
    "tacrolimus_tdm_patients",
    "tacrolimus_timed_dose_concentration_patients",
    "tacrolimus_preoutcome_omics_patients",
    "tacrolimus_cyp3a5_patients",
    "tacrolimus_repeated_omics_patients",
    "mpa_tdm_patients",
    "mpa_timed_dose_concentration_patients",
    "mpa_preoutcome_omics_patients",
    "independent_site_or_time_holdout_patients",
)

ATTESTATION_FIELDS = (
    "deterministic_same_person_crosswalk_confirmed",
    "authorized_linkage_and_research_use_confirmed",
    "omics_collection_time_available",
    "no_direct_identifiers_in_aggregate_exchange",
)


def _count(value: object, field: str) -> int:
    if isinstance(value,bool):
        raise BiomarginAbisFeasibilityError(f"{field} must be an integer count")
    try:
        numeric=float(value)
    except (TypeError,ValueError):
        raise BiomarginAbisFeasibilityError(f"{field} must be an integer count") from None
    if not numeric.is_integer() or numeric<0:
        raise BiomarginAbisFeasibilityError(f"{field} must be a nonnegative integer")
    return int(numeric)


def _attestation(value: object, field: str) -> bool:
    if not isinstance(value,bool):
        raise BiomarginAbisFeasibilityError(f"{field} must be boolean")
    return value


def qualify_biomargin_abis_aggregate_feasibility(
    counts: Mapping[str,object],
    attestations: Mapping[str,object],
) -> dict:
    """Evaluate linkage feasibility without receiving any patient identifiers.

    Counts must refer to unique biological patients after exact deterministic
    source-system matching. They are not sample/assay/visit counts.

    The function does not assert that a cohort is statistically powered. It
    only determines whether the minimum source linkage needed for future
    Omics2PK work exists at all.
    """
    missing=[x for x in COUNT_FIELDS if x not in counts]
    extra=[x for x in counts if x not in COUNT_FIELDS]
    if missing or extra:
        raise BiomarginAbisFeasibilityError(
            "Counts must use the frozen aggregate feasibility schema")
    missing_att=[x for x in ATTESTATION_FIELDS if x not in attestations]
    extra_att=[x for x in attestations if x not in ATTESTATION_FIELDS]
    if missing_att or extra_att:
        raise BiomarginAbisFeasibilityError(
            "Attestations must use the frozen aggregate feasibility schema")

    c={field:_count(counts[field],field) for field in COUNT_FIELDS}
    a={field:_attestation(attestations[field],field) for field in ATTESTATION_FIELDS}

    total=c["biomargin_molecular_patients"]
    matched=c["exact_biomargin_abis_crosswalk_patients"]
    if matched>total:
        raise BiomarginAbisFeasibilityError(
            "Exact ABIS crosswalk cannot exceed BIOMARGIN molecular patients")

    def hierarchy(values, label):
        if any(right>left for left,right in zip(values,values[1:])):
            raise BiomarginAbisFeasibilityError(
                f"{label} counts violate nested-patient hierarchy")

    hierarchy((
        matched,
        c["tacrolimus_tdm_patients"],
        c["tacrolimus_timed_dose_concentration_patients"],
        c["tacrolimus_preoutcome_omics_patients"],
        c["tacrolimus_repeated_omics_patients"],
    ),"Tacrolimus")
    hierarchy((
        matched,
        c["mpa_tdm_patients"],
        c["mpa_timed_dose_concentration_patients"],
        c["mpa_preoutcome_omics_patients"],
    ),"MPA")

    if c["tacrolimus_cyp3a5_patients"]>c["tacrolimus_timed_dose_concentration_patients"]:
        raise BiomarginAbisFeasibilityError(
            "Tacrolimus CYP3A5 count cannot exceed timed-PK matched patients")
    if c["independent_site_or_time_holdout_patients"]>matched:
        raise BiomarginAbisFeasibilityError(
            "Holdout-eligible count cannot exceed exact crosswalk patients")

    if matched>0 and not a["deterministic_same_person_crosswalk_confirmed"]:
        raise BiomarginAbisFeasibilityError(
            "Positive crosswalk count requires deterministic crosswalk confirmation")
    if matched>0 and not a["authorized_linkage_and_research_use_confirmed"]:
        raise BiomarginAbisFeasibilityError(
            "Positive crosswalk count requires linkage/research authorization")
    if (
        c["tacrolimus_preoutcome_omics_patients"]>0
        or c["mpa_preoutcome_omics_patients"]>0
    ) and not a["omics_collection_time_available"]:
        raise BiomarginAbisFeasibilityError(
            "Pre-outcome omics count requires molecular collection timing")
    if not a["no_direct_identifiers_in_aggregate_exchange"]:
        raise BiomarginAbisFeasibilityError(
            "Aggregate feasibility exchange must exclude direct identifiers")

    tac_pairing=(
        a["deterministic_same_person_crosswalk_confirmed"]
        and a["authorized_linkage_and_research_use_confirmed"]
        and c["tacrolimus_timed_dose_concentration_patients"]>0
    )
    tac_preoutcome=tac_pairing and c["tacrolimus_preoutcome_omics_patients"]>0
    tac_longitudinal=tac_preoutcome and c["tacrolimus_repeated_omics_patients"]>0
    mpa_pairing=(
        a["deterministic_same_person_crosswalk_confirmed"]
        and a["authorized_linkage_and_research_use_confirmed"]
        and c["mpa_timed_dose_concentration_patients"]>0
    )
    mpa_preoutcome=mpa_pairing and c["mpa_preoutcome_omics_patients"]>0

    if matched==0:
        status="NO_BIOMARGIN_ABIS_PATIENT_CROSSWALK"
    elif not (tac_pairing or mpa_pairing):
        status="CROSSWALK_EXISTS_BUT_NO_TIMED_DRUG_PK_PAIRING"
    elif tac_preoutcome or mpa_preoutcome:
        status="PAIRED_OMICS_TDM_FEASIBILITY_CONFIRMED_AT_AGGREGATE_LEVEL"
    else:
        status="TIMED_TDM_PAIRING_EXISTS_BUT_PREOUTCOME_OMICS_NOT_CONFIRMED"

    return {
        "status":status,
        "counts":c,
        "tacrolimus_timed_pairing_exists":tac_pairing,
        "tacrolimus_preoutcome_omics_pairing_exists":tac_preoutcome,
        "tacrolimus_prospective_longitudinal_candidate_exists":tac_longitudinal,
        "mpa_timed_pairing_exists":mpa_pairing,
        "mpa_preoutcome_omics_pairing_exists":mpa_preoutcome,
        "cyp3a5_available_with_tacrolimus_timed_pk":
            c["tacrolimus_cyp3a5_patients"]>0,
        "independent_site_or_time_holdout_candidate_exists":
            c["independent_site_or_time_holdout_patients"]>0,
        "aggregate_counts_are_not_statistical_power_validation":True,
        "patient_level_data_requested_by_this_gate":False,
        "biomargin_itself_contains_pk":False,
        "abis_crosswalk_publicly_proven_by_this_gate":False,
    }
