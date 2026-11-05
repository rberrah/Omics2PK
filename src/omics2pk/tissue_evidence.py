"""Aggregate-only assay specimen versus PK target-organ provenance.

Blood-derived RNA/protein/metabolite may be a useful INDIRECT biomarker, but
is not an experimentally validated liver/intestine enzyme activity by default.
No individual records/identifiers, private values, source location names or
genotype are returned. The audit does not prohibit a blood-derived statistical
predictor; it prohibits false mechanical *tissue equivalence* claims.
"""
from __future__ import annotations

BLOOD={"whole_blood","blood","peripheral_blood","pbmc","plasma","serum","extracellular_vesicles"}
LIVER={"liver","liver_biopsy","liver_microsomes","liver_s9","hepatic_tissue","hepatocytes"}
INTESTINE={"intestine","small_intestine","enterocytes","ileum","jejunum"}
VALID_DOSE={"pre_dose","post_dose","at_dose","not_applicable","unknown"}
VALID_ASSAY={"bulk","single_cell","single_nucleus","sorted_cells","cell_line","cell_free","unknown"}

class TissueEvidenceError(ValueError):
    pass

def canonical(v):
    return str(v or "").strip().lower().replace("-","_").replace(" ","_")

def audit_tissue_source_provenance(samples):
    if not isinstance(samples,(list,tuple)) or not samples:
        raise TissueEvidenceError("Require source assays before tissue qualification")
    n=len(samples)
    declared_tissue=0
    declared_specimen=0
    declared_cellular=0
    validated_organ_proxy=0  # CANNOT be proved by source labels
    cross_organ_candidates=0
    postdose=0
    predose=0
    mixed_bulk=0
    for sample in samples:
        if not isinstance(sample,dict):
            raise TissueEvidenceError("Sample metadata must be original records")
        tissue=canonical(sample.get("tissue"))
        specimen=canonical(sample.get("specimen_type"))
        matrix=canonical(sample.get("sample_matrix"))
        cellular=canonical(sample.get("cell_type_or_bulk"))
        organ=canonical(sample.get("pharmacology_target_organ"))
        relative=canonical(sample.get("relative_to_dose"))
        declared_tissue+=bool(tissue)
        declared_specimen+=bool(specimen)
        declared_cellular+=bool(cellular)
        if cellular and cellular not in VALID_ASSAY:
            raise TissueEvidenceError("Invalid original cell type / bulk declaration")
        if relative and relative not in VALID_DOSE:
            raise TissueEvidenceError("Invalid original relative-to-dose declaration")
        if relative=="post_dose":postdose+=1
        if relative=="pre_dose":predose+=1
        source_is_blood=any(v in BLOOD for v in (tissue,specimen,matrix))
        target_is_metabolizing_organ=organ in LIVER|INTESTINE
        cross_organ_candidates+=bool(source_is_blood and target_is_metabolizing_organ)
        mixed_bulk+=bool(cellular=="bulk" and canonical(sample.get("cell_composition_status")) not in ("measured","estimated"))
    return {
        "status":("SOURCE_TISSUE_AND_SPECIMEN_NOT_FULLY_DECLARED"
                  if declared_tissue<n or declared_specimen<n else
                  "ORIGINAL_TISSUE_DECLARED_BUT_CROSS_ORGAN_FUNCTION_NOT_PROVED"),
        "n_omic_assays":n,
        "n_original_tissue_declared":int(declared_tissue),
        "n_original_specimen_type_declared":int(declared_specimen),
        "n_cell_assay_mode_declared":int(declared_cellular),
        "n_bulk_assays_without_measured_or_estimated_cell_composition":int(mixed_bulk),
        "n_original_predose_declared":int(predose),
        "n_original_postdose_declared":int(postdose),
        "n_blood_vs_hepatic_or_intestinal_target_candidates_unvalidated":int(cross_organ_candidates),
        "n_cross_organ_biomarker_proxies_functionally_validated_by_metadata_alone":validated_organ_proxy,
        "blood_PBMC_and_plasma_are_distinct_sources":True,
        "whole_blood_omics_can_be_qualified_direct_hepatic_CYP3A_activity":False,
        "statistical_use_of_blood_as_indirect_Pk_predictor_is_not_prohibited":True,
        "source_omics_measurements_and_patient_ids_not_returned":True,
    }
