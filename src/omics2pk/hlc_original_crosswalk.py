"""A strictly structural same-original-donor HLC GEO↔Synapse source gate.

One original GSM microarray is NOT by itself a biological donor. This gate
requires an independently sourced codebook crosswalk and an actual measured
original enzyme assay phenotype per liver, never sorting rows or assuming
427 expression + 398 activity observations are the same humans.

This only returns privacy-safe aggregate candidate overlap. It NEVER makes
a statistical prediction or certifies that a caller's crosswalk source file
is a real author-issued biological donor codebook.
"""
from __future__ import annotations

from collections import Counter
from math import isfinite
import re

GSM=re.compile(r"^GSM\d{5,9}$")
HASH=re.compile(r"^[0-9a-f]{64}$")
REQUIRED_CROSSWALK=("gsm_id","donor_id","codebook_sha256")
REQUIRED_ACTIVITY=("donor_id","assay_id","measured_value","assay_units","original_entity_id")


class HumanLiverCrosswalkError(ValueError):
    pass


def _text(value):
    return isinstance(value,str) and bool(value.strip()) and len(value)<150


def audit_explicit_hlc_crosswalk(
    public_gsm_ids:list[str],
    original_codebook_rows:list[dict],
    original_phenotype_rows:list[dict],
    *,
    assay_id:str,
)->dict:
    if len(public_gsm_ids)!=427 or len(set(public_gsm_ids))!=427 or not all(
        GSM.fullmatch(x) for x in public_gsm_ids
    ):
        raise HumanLiverCrosswalkError("Require exactly 427 original GEO GSE9588 unique GSM IDs")
    if not _text(assay_id):
        raise HumanLiverCrosswalkError("Specify actual source enzyme probe assay")
    if not original_codebook_rows:
        raise HumanLiverCrosswalkError("Original documented donor crosswalk is missing")
    if not original_phenotype_rows:
        raise HumanLiverCrosswalkError("Original measured CYP activity dataset is missing")
    gsm_set=set(public_gsm_ids)
    assigned={}
    reverse=set()
    codebook_hashes=set()
    for row in original_codebook_rows:
        if not all(k in row for k in REQUIRED_CROSSWALK):
            raise HumanLiverCrosswalkError("Original donor key metadata missing")
        gsm=row["gsm_id"];donor=row["donor_id"];digest=row["codebook_sha256"]
        if not GSM.fullmatch(gsm) or gsm not in gsm_set or not _text(donor) or not HASH.fullmatch(digest):
            raise HumanLiverCrosswalkError("Original donor GSM/source hash inconsistent")
        if gsm in assigned or donor in reverse:
            raise HumanLiverCrosswalkError("Repeated GEO sample or biological donor mapping")
        assigned[gsm]=donor;reverse.add(donor);codebook_hashes.add(digest)
    if len(codebook_hashes)!=1:
        raise HumanLiverCrosswalkError("Crosswalk must have one original source document checksum")
    seen=set()
    targets={}
    excluded=Counter()
    unit_types=set()
    for row in original_phenotype_rows:
        if not all(k in row for k in REQUIRED_ACTIVITY):
            raise HumanLiverCrosswalkError("Original enzymatic assay metadata missing")
        if row["original_entity_id"]!="syn4629":
            raise HumanLiverCrosswalkError("Phenotype must be original Synapse syn4629 source")
        if row["assay_id"]!=assay_id:
            continue
        donor=row["donor_id"]
        if not _text(donor) or donor in seen:
            raise HumanLiverCrosswalkError("Repeated or missing original assay biological donor")
        seen.add(donor)
        if not _text(row["assay_units"]):
            raise HumanLiverCrosswalkError("Source enzyme activity assay physical unit missing")
        unit_types.add(row["assay_units"].strip())
        try:measurement=float(row["measured_value"])
        except (ValueError,OverflowError,TypeError):
            excluded["non_numeric"]+=1
            continue
        if not isfinite(measurement) or measurement<0:
            excluded["invalid_or_censored_without_LOQ_definition"]+=1
            continue
        targets[donor]=measurement
    if len(unit_types)!=1:
        raise HumanLiverCrosswalkError("Original enzyme assay physical units inconsistent")
    matched=sum(donor in targets for donor in reverse)
    return {
        "status":"PROVISIONAL_ORIGINAL_HLC_SOURCE_MATCH_COUNT_ONLY",
        "GEO_GSE9588_unique_original_GSM_samples":len(public_gsm_ids),
        "source_crosswalk_unique_original_GSM_donor_pairs":len(assigned),
        "source_phenotype_subjects_with_numeric_selected_enzyme_activity":len(targets),
        "source_original_GSM_donor_matches_with_assay":matched,
        "source_original_activity_units_consistent":True,
        "n_unquantifiable_or_invalid_assays":sum(excluded.values()),
        "original_source_codebook_SHA256_present":True,
        "original_source_codebook_proven_authentic_by_this_function":False,
        "qualified_independent_biological_donor_numeric_model":False,
        "requires_independent_primary_source_codebook_and_assay_unit_review":True,
        "original_donor_keys_and_probe_activities_exported":False,
    }
