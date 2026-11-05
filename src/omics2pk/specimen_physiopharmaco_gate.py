"""Public-data-only anatomical sampling-site × physiology qualification.

The *original assayed specimen* (plasma, blood cells, liver tissue,
isolated hepatocytes, microsomes, cultures) is distinct from the
*biological target organ* in which a drug is cleared. A metabolite
in plasma is an integrated systemic proxy, not a liver enzyme biopsy.
A blood RNA count is not a CYP3A4/5 liver RNA count.

This is a pre-outcome study/assay eligibility & interpretation gate.
It never trains a model, infers an unmeasured layer, predicts an
enzyme's activity from blood, or claims Omics2PK benefit vs full PGx.
"""
from __future__ import annotations
from collections.abc import Mapping

class PhysiopharmacoSpecimenError(ValueError):
    pass

PROCESS={
    "tacrolimus":{
        "hepatic_CYP3A_metabolism":{"liver"},
        "intestinal_CYP3A_first_pass":{"small_intestine"},
        "whole_blood_erythrocyte_partition":{"blood"},
    },
    "metformin":{
        "renal_transporter_secretion":{"kidney"},
        "intestinal_transport_and_response":{"small_intestine"},
    },
    "imatinib":{
        "hepatic_CYP_metabolism":{"liver"},
    },
    "doxorubicin":{
        "cancer_cellular_PD":{"cancer_cells"},
    },
}
# Explicit direct molecular source: an assay in an organ where the process
# happens may still be ex vivo; direct in-tissue is NOT clinical drug PK.
SOURCES={
    "HLC_LIVER_RNA":{
        "original_sample":"bulk_liver_tissue",
        "original_organ":"liver","omic":"transcriptomics",
        "method":"Agilent_liver_microarray_gene_probe",
        "collection":"ex_vivo_baseline",
        "source_grade":"OWNER_AUTHORIZED_REPLAY_NOT_OPEN_RAW",
        "has_direct_patient_drug_PK":False,
        "has_full_PGx_panel":False,
        "original_reference":"Yang_2010_sage_syn4499",
    },
    "HLC_LIVER_MICROSOME_CYP_PROBE":{
        "original_sample":"liver_microsomes",
        "original_organ":"liver","omic":"functional_activity",
        "method":"CYP_microsomal_probe_substrate_activity",
        "collection":"ex_vivo_baseline",
        "source_grade":"OWNER_AUTHORIZED_REPLAY_NOT_OPEN_RAW",
        "has_direct_patient_drug_PK":False,
        "has_full_PGx_panel":False,
        "original_reference":"Yang_2010_sage_syn4499",
    },
    "LATHAM_PHH_CYP_PROTEIN":{
        "original_sample":"primary_hepatocytes",
        "original_organ":"liver","omic":"proteomics",
        "method":"QTAP_targeted_CYP_peptides",
        "collection":"ex_vivo_baseline",
        "source_grade":"ORIGINAL_PUBLIC_DONOR_DATA_REPLAYED",
        "has_direct_patient_drug_PK":False,
        "has_full_PGx_panel":False,
        "original_reference":"Latham_2025_PMC12990092_Table4",
    },
    "LATHAM_PHH_SHARED_PROBE":{
        "original_sample":"primary_hepatocytes",
        "original_organ":"liver","omic":"functional_activity",
        "method":"CYP2C8_and_shared_CYP3A_probe",
        "collection":"ex_vivo_baseline",
        "source_grade":"ORIGINAL_PUBLIC_DONOR_DATA_REPLAYED",
        "has_direct_patient_drug_PK":False,
        "has_full_PGx_panel":False,
        "original_reference":"Latham_2025_PMC12990092_Table4",
    },
    "ROTROFF_PLASMA_GCTOF_A":{
        "original_sample":"plasma",
        "original_organ":"blood","omic":"metabolomics",
        "method":"GC_TOF_untargeted",
        "collection":"untreated_predose",
        "source_grade":"ORIGINAL_PUBLIC_DONOR_DATA_REPLAYED",
        "has_direct_patient_drug_PK":True,
        "has_full_PGx_panel":False,
        "original_reference":"Rotroff_2016_PMC4906013",
    },
    "ROTROFF_PLASMA_GCTOF_B":{
        "original_sample":"plasma",
        "original_organ":"blood","omic":"metabolomics",
        "method":"GC_TOF_untargeted",
        "collection":"post_first_dose",
        "source_grade":"ORIGINAL_PUBLIC_DONOR_DATA_REPLAYED",
        "has_direct_patient_drug_PK":True,
        "has_full_PGx_panel":False,
        "original_reference":"Rotroff_2016_PMC4906013",
    },
    "ROTROFF_PLASMA_GCTOF_C":{
        "original_sample":"plasma",
        "original_organ":"blood","omic":"metabolomics",
        "method":"GC_TOF_untargeted",
        "collection":"post_second_dose",
        "source_grade":"ORIGINAL_PUBLIC_DONOR_DATA_REPLAYED",
        "has_direct_patient_drug_PK":True,
        "has_full_PGx_panel":False,
        "original_reference":"Rotroff_2016_PMC4906013",
    },
    "EATRIS_WHITE_BLOOD_CELLS_RNA":{
        "original_sample":"white_blood_cells",
        "original_organ":"blood","omic":"transcriptomics",
        "method":"WBC_mRNAseq",
        "collection":"baseline_healthy_no_drug_PK",
        "source_grade":"PUBLIC_PREPROCESSED_MULTIVIEW_NO_PUBLIC_FULL_GERM_LINE_OR_PK",
        "has_direct_patient_drug_PK":False,
        "has_full_PGx_panel":False,
        "original_reference":"EATRIS_Plus_v4_Zenodo_19136013",
    },
    "EATRIS_PLASMA_PROTEIN":{
        "original_sample":"plasma",
        "original_organ":"blood","omic":"proteomics",
        "method":"plasma_proteomic_panel",
        "collection":"baseline_healthy_no_drug_PK",
        "source_grade":"PUBLIC_PREPROCESSED_MULTIVIEW_NO_PUBLIC_FULL_GERM_LINE_OR_PK",
        "has_direct_patient_drug_PK":False,
        "has_full_PGx_panel":False,
        "original_reference":"EATRIS_Plus_v4_Zenodo_19136013",
    },
    "EATRIS_PLASMA_METABOLITES":{
        "original_sample":"plasma",
        "original_organ":"blood","omic":"metabolomics",
        "method":"targeted_plasma_metabolites",
        "collection":"baseline_healthy_no_drug_PK",
        "source_grade":"PUBLIC_PREPROCESSED_MULTIVIEW_NO_PUBLIC_FULL_GERM_LINE_OR_PK",
        "has_direct_patient_drug_PK":False,
        "has_full_PGx_panel":False,
        "original_reference":"EATRIS_Plus_v4_Zenodo_19136013",
    },
    "GTEX_V8_BULK_LIVER_RNA_MEDIAN":{
        "original_sample":"bulk_liver_tissue",
        "original_organ":"liver","omic":"transcriptomics",
        "method":"bulk_RNAseq_median_TPM",
        "collection":"reference_atlas_aggregate",
        "source_grade":"ORIGINAL_PUBLIC_TISSUE_MEDIAN_REPLAYED_NO_DONOR_VALUES",
        "has_direct_patient_drug_PK":False,
        "has_full_PGx_panel":False,
        "original_reference":"GTEx_V8_2017_DCC",
    },
    "GTEX_V8_WHOLE_BLOOD_RNA_MEDIAN":{
        "original_sample":"whole_blood",
        "original_organ":"blood","omic":"transcriptomics",
        "method":"bulk_RNAseq_median_TPM",
        "collection":"reference_atlas_aggregate",
        "source_grade":"ORIGINAL_PUBLIC_TISSUE_MEDIAN_REPLAYED_NO_DONOR_VALUES",
        "has_direct_patient_drug_PK":False,
        "has_full_PGx_panel":False,
        "original_reference":"GTEx_V8_2017_DCC",
    },
    "NCI60_RNA":{
        "original_sample":"cancer_cell_line",
        "original_organ":"cancer_cells","omic":"transcriptomics",
        "method":"cell_line_gene_expression_different_study",
        "collection":"cross_study_cultures",
        "source_grade":"ORIGINAL_PUBLIC_CELL_LINE_JOIN_NOT_SAME_SPECIMEN",
        "has_direct_patient_drug_PK":False,
        "has_full_PGx_panel":False,
        "original_reference":"IntLIM_NCI60_Mathelab",
    },
    "NCI60_SWATH_PROTEOMICS":{
        "original_sample":"cancer_cell_line",
        "original_organ":"cancer_cells","omic":"proteomics",
        "method":"Guo_2019_original_SWATH_MS",
        "collection":"cross_study_cultures",
        "source_grade":"ORIGINAL_PUBLIC_CELL_LINE_JOIN_NOT_SAME_SPECIMEN",
        "has_direct_patient_drug_PK":False,
        "has_full_PGx_panel":False,
        "original_reference":"Guo_2019_PMC6889472",
    },
    "NCI60_METABOLITES":{
        "original_sample":"cancer_cell_line",
        "original_organ":"cancer_cells","omic":"metabolomics",
        "method":"IntLIM_original_named_metabolites",
        "collection":"cross_study_cultures",
        "source_grade":"ORIGINAL_PUBLIC_CELL_LINE_JOIN_NOT_SAME_SPECIMEN",
        "has_direct_patient_drug_PK":False,
        "has_full_PGx_panel":False,
        "original_reference":"IntLIM_NCI60_Mathelab",
    },
    "CYP3A5_GEO_PHH_RNA":{
        "original_sample":"primary_hepatocytes",
        "original_organ":"liver","omic":"transcriptomics",
        "method":"GSE124074_original_RNAseq",
        "collection":"ex_vivo_baseline",
        "source_grade":"ORIGINAL_PUBLIC_SAMPLE_METADATA_CROSSWALK_ONLY",
        "has_direct_patient_drug_PK":False,
        "has_full_PGx_panel":False,
        "original_reference":"Zhong_2020_GSE124074_GSE147628",
    },
    "CYP2C9_DMS_VAMPSEQ":{
        "original_sample":"engineered_human_cell_culture",
        "original_organ":"engineered_cells","omic":"proteomics",
        "method":"VAMPseq_variant_abundance",
        "collection":"ex_vivo_construct_screen",
        "source_grade":"ORIGINAL_PUBLIC_ENGINEERED_VARIANT_REPLAYED",
        "has_direct_patient_drug_PK":False,
        "has_full_PGx_panel":False,
        "original_reference":"Amorosi_2021_AJHG",
    },
    "CYP2C9_DMS_CLICKSEQ":{
        "original_sample":"engineered_yeast_culture",
        "original_organ":"engineered_cells","omic":"functional_activity",
        "method":"Clickseq_variant_activity",
        "collection":"ex_vivo_construct_screen",
        "source_grade":"ORIGINAL_PUBLIC_ENGINEERED_VARIANT_REPLAYED",
        "has_direct_patient_drug_PK":False,
        "has_full_PGx_panel":False,
        "original_reference":"Amorosi_2021_AJHG",
    },
    "NUDT15_DMS_HEK293T_PROTEIN":{
        "original_sample":"engineered_human_cell_culture",
        "original_organ":"engineered_cells","omic":"proteomics",
        "method":"VAMPseq_variant_abundance",
        "collection":"ex_vivo_construct_screen",
        "source_grade":"ORIGINAL_PUBLIC_ENGINEERED_VARIANT_REPLAYED",
        "has_direct_patient_drug_PK":False,
        "has_full_PGx_panel":False,
        "original_reference":"Suiter_2020_PNAS",
    },
    "NUDT15_DMS_HEK293T_TOLERANCE":{
        "original_sample":"engineered_human_cell_culture",
        "original_organ":"engineered_cells","omic":"cellular_PD",
        "method":"thiopurine_cell_viability",
        "collection":"postdrug_response",
        "source_grade":"ORIGINAL_PUBLIC_ENGINEERED_VARIANT_REPLAYED",
        "has_direct_patient_drug_PK":False,
        "has_full_PGx_panel":False,
        "original_reference":"Suiter_2020_PNAS",
    },
}
TISSUE_LEVEL={"bulk_liver_tissue","liver_microsomes","primary_hepatocytes",
              "small_intestinal_tissue","kidney_tissue"}
FLUID_LEVEL={"plasma","serum","urine"}
BLOOD_CELL={"PBMC","white_blood_cells","whole_blood"}
POST={"post_first_dose","post_second_dose","postdrug_response"}
NONPREDICTIVE={"reference_atlas_aggregate","cross_study_cultures","postdrug_response"}
RESEARCH={"ex_vivo_baseline","ex_vivo_construct_screen","cross_study_cultures"}


def qualify_public_specimen_for_PK(*,source_id,drug,process,
                                   intended_prediction="baseline_predose",
                                   require_real_human_drug_PK=False):
    """Return physical-origin relevance and separate source evidence limits.

    Never infer individual clinical PK from an atlas, in-vitro probe,
    a published group mean, a cell line, or genetically-engineered clone.
    """
    if source_id not in SOURCES:
        raise PhysiopharmacoSpecimenError("Unknown original independently qualified specimen source")
    if drug not in PROCESS or process not in PROCESS[drug]:
        raise PhysiopharmacoSpecimenError("PK process not predeclared for original drug")
    if intended_prediction not in {"baseline_predose","ex_vivo_function",
                                   "on_treatment_future_landmark"}:
        raise PhysiopharmacoSpecimenError("Original prediction landmark not qualified")
    if not isinstance(require_real_human_drug_PK,bool):
        raise PhysiopharmacoSpecimenError("Explicit real PK evidence requirement flag")
    s=SOURCES[source_id]
    original=s["original_organ"]
    expected=PROCESS[drug][process]
    if original in expected:
        if s["original_sample"] in TISSUE_LEVEL:
            relevance="SOURCE_IN_TARGET_ORGAN_BUT_NOT_PROOF_IN_VIVO_FUNCTION"
        elif s["original_sample"] in BLOOD_CELL:
            relevance="SOURCE_BLOOD_COMPARTMENT_SAME_ORGAN_NOT_SAME_DRUG_PARTITION"
        else:
            relevance="ORIGINAL_ORGAN_MATCH_WITH_UNQUALIFIED_CELL_FRACTION"
    elif s["original_sample"] in FLUID_LEVEL:
        relevance="SYSTEMIC_BIOFLUID_PROXY_NOT_TISSUE_SPECIFIC_ENZYME_FUNCTION"
    elif s["original_sample"] in BLOOD_CELL:
        relevance="BLOOD_CELL_RNA_NOT_HEPATIC_OR_RENAL_TARGET_TISSUE"
    elif s["original_sample"].startswith("engineered"):
        relevance="ENGINEERED_VARIANT_CONSTRUCT_NOT_REAL_PATIENT_ORGAN"
    elif s["original_sample"]=="cancer_cell_line":
        relevance="CANCER_CELL_LINE_ONLY_NOT_PATIENT_ORGAN"
    else:
        relevance="DIFFERENT_ORGAN_NOT_DIRECT_TISSUE_BIOLOGY"
    if s["collection"] in POST and intended_prediction=="baseline_predose":
        temporal="POSTDOSE_NOT_ELIGIBLE_AS_BASELINE_PREDICTOR"
    elif s["collection"]=="reference_atlas_aggregate":
        temporal="REFERENCE_ATLAS_MEDIAN_NOT_INDIVIDUAL_BASELINE"
    elif s["collection"]=="cross_study_cultures":
        temporal="SOURCE_CELL_LINE_CULTURE_BATCH_NOT_PATIENT_TIME"
    elif s["collection"] in RESEARCH:
        temporal="EX_VIVO_ASSAY_NOT_A_SAMPLED_PATIENT_PREDOSE_MARKER"
    elif s["collection"]=="untreated_predose":
        temporal="MEASURED_PRETREATMENT_SYSTEMIC_SPECIMEN"
    else:
        temporal="COLLECTION_TIME_OR_LANDMARK_CONTEXT_NOT_INDIVIDUAL_PK_VALIDATED"
    if source_id=="ROTROFF_PLASMA_GCTOF_B" and intended_prediction=="on_treatment_future_landmark":
        temporal="POST_FIRST_DOSE_ONLY_POTENTIALLY_ELIGIBLE_FOR_LATER_FUTURE_OUTCOME_WITH_TIMING_QC"
    if source_id=="ROTROFF_PLASMA_GCTOF_C" and intended_prediction=="on_treatment_future_landmark":
        temporal="POST_SECOND_DOSE_ONLY_POTENTIALLY_ELIGIBLE_FOR_LATER_FUTURE_OUTCOME_WITH_TIMING_QC"
    is_source_real=bool(s["has_direct_patient_drug_PK"])
    can_fit=(
        temporal in {"MEASURED_PRETREATMENT_SYSTEMIC_SPECIMEN"}
        and is_source_real
        and s["source_grade"]=="ORIGINAL_PUBLIC_DONOR_DATA_REPLAYED"
    )
    if intended_prediction=="on_treatment_future_landmark":
        can_fit=False  # individual sample time vs later endpoint has not been qualified
    if require_real_human_drug_PK and not is_source_real:
        can_fit=False
    # Actual source has no complete clinical+PGx reference on any
    # listed study, including original Rotroff's true clinical PK.
    return {
        "status":"SOURCE_SPECIMEN_SITE_QUALIFIED_NOT_PREDICTIVE_BENEFIT",
        "source":source_id,"original_reference":s["original_reference"],
        "measured_source_omic":s["omic"],
        "actual_collected_fraction":s["original_sample"],
        "actual_anatomical_collection_organ":original,
        "assay_method":s["method"],
        "collection_landmark":s["collection"],
        "intended_PK_process":process,
        "candidate_target_organ_anatomy":sorted(expected),
        "mechanistic_relation_of_sample_to_process":relevance,
        "sample_timing_interpretation":temporal,
        "public_original_source_qualification":s["source_grade"],
        "original_paired_human_drug_PK_source_present":is_source_real,
        "full_PGx_only_model_benchmark_possible_now":bool(s["has_full_PGx_panel"]),
        "eligible_as_baseline_individual_PK_omic_proxy_under_source_contract":can_fit,
        "source_is_direct_measurement_of_expected_target_organ_function":
            relevance=="SOURCE_IN_TARGET_ORGAN_BUT_NOT_PROOF_IN_VIVO_FUNCTION"
            and s["omic"]=="functional_activity",
        "molecule_in_plasma_does_not_establish_unique_organ_or_enzyme":True,
        "blood_CYP3A_gene_RNA_not_liver_CYP3A_capacity":True,
        "cross_source_subject_PK_and_omics_join_independently_established":False,
        "patient_PK_improvement_vs_complete_genotype_proven":False,
        "original_subject_rows_and_matrices_exported":False,
    }


def qualify_public_omic_model_arms(*,arm_to_original_sources,
                                   drug,process,endpoint_class):
    """Guard against joining molecular arms from unrelated public cohorts.

    This is study/source metadata qualification ONLY. Even if two
    modalities originate from the same article, the original individual
    patient ID match, assay precision and collection time must still be
    checked by the primary-source predictor. No feature matrix generated.
    """
    if endpoint_class not in {
            "CLINICAL_PK","IN_VITRO_DRUG_METABOLISM_NOT_CLINICAL_PK",
            "BIOLOGICAL_FUNCTION_NOT_CLINICAL_PK","IN_VITRO_CELL_PD_NOT_PATIENT"}:
        raise PhysiopharmacoSpecimenError("Original independent phenotype class unqualified")
    if not isinstance(arm_to_original_sources,Mapping) or (
            "clinical_plus_genotype" not in arm_to_original_sources):
        raise PhysiopharmacoSpecimenError("Original clinical+measured genotype baseline mandatory")
    if arm_to_original_sources["clinical_plus_genotype"] not in ([],()):
        raise PhysiopharmacoSpecimenError("Clinical+genotype model arm cannot silently contain added omics")
    if len(arm_to_original_sources)<2 or len(arm_to_original_sources)>20:
        raise PhysiopharmacoSpecimenError("Original source omic comparison not specified")
    used=set();descriptions={}
    for arm,names in sorted(arm_to_original_sources.items()):
        if not isinstance(arm,str) or not arm.strip():
            raise PhysiopharmacoSpecimenError("Invalid original model arm identity")
        if not isinstance(names,(list,tuple)) or (arm!="clinical_plus_genotype" and not names):
            raise PhysiopharmacoSpecimenError("Every added omic model needs an actual measured source")
        if len(names)!=len(set(names)):
            raise PhysiopharmacoSpecimenError("Duplicate original assay in one omic arm")
        assays=[]
        for name in names:
            q=qualify_public_specimen_for_PK(
                source_id=name,drug=drug,process=process,
                intended_prediction="ex_vivo_function" if endpoint_class!= "CLINICAL_PK"
                else "baseline_predose",
                require_real_human_drug_PK=(endpoint_class=="CLINICAL_PK"),
            )
            if q["sample_timing_interpretation"] in {
                    "POSTDOSE_NOT_ELIGIBLE_AS_BASELINE_PREDICTOR",
                    "REFERENCE_ATLAS_MEDIAN_NOT_INDIVIDUAL_BASELINE",
                    "SOURCE_CELL_LINE_CULTURE_BATCH_NOT_PATIENT_TIME"}:
                raise PhysiopharmacoSpecimenError(
                    "Cannot insert post-dose, atlas median or cross-study culture profiles as individual predose predictors")
            if endpoint_class=="CLINICAL_PK":
                if not q["eligible_as_baseline_individual_PK_omic_proxy_under_source_contract"]:
                    raise PhysiopharmacoSpecimenError(
                        "No real same-original-patient predose omic and clinical PK join qualified by this public source")
            elif q["public_original_source_qualification"] not in {
                    "ORIGINAL_PUBLIC_DONOR_DATA_REPLAYED",
                    "OWNER_AUTHORIZED_REPLAY_NOT_OPEN_RAW",
            }:
                raise PhysiopharmacoSpecimenError(
                    "In-vitro function assay requires actual original individual human donor data")
            assays.append(q)
            used.add(name)
        descriptions[arm]=assays
    unique_refs={SOURCES[k]["original_reference"] for k in used}
    if len(unique_refs)!=1:
        raise PhysiopharmacoSpecimenError(
            "Original omic arms from different studies cannot be matched as one biological cohort")
    return {
        "status":"QUALIFIED_ASSAY_TISSUE_SITE_AND_ORIGINAL_STUDY_FOR_PRESPECIFIED_ARMS",
        "drug":drug,"biological_process":process,
        "endpoint_class":endpoint_class,
        "source_public_original_reference":next(iter(unique_refs)),
        "n_distinct_original_measured_assays":len(used),
        "arm_original_collection_sites":descriptions,
        "source_site_does_not_prove_actual_original_subject_time_join":True,
        "no_added_omic_accuracy_or_causal_effect_inferred":True,
        "no_clinical_superiority_over_complete_genotype_claim":True,
    }
