"""Synthetic contracts for provenance; source sample sites are curated facts."""
import unittest
from omics2pk.specimen_physiopharmaco_gate import (
    SOURCES,qualify_public_specimen_for_PK,qualify_public_omic_model_arms,
    qualify_pharmacogenetic_sampling_origin,PhysiopharmacoSpecimenError,
)

def q(s,drug,process,**kw):
    return qualify_public_specimen_for_PK(source_id=s,drug=drug,process=process,**kw)

class SourceTissuePhysioGate(unittest.TestCase):
    def test_curated_public_source_all_explicitly_organ_method_landmark(self):
        self.assertGreaterEqual(len(SOURCES),15)
        for source,v in SOURCES.items():
            with self.subTest(source=source):
                self.assertTrue(v["original_sample"])
                self.assertTrue(v["original_organ"])
                self.assertTrue(v["method"])
                self.assertTrue(v["collection"])
                self.assertFalse(v["has_full_PGx_panel"])

    def test_GTex_actual_original_v8_liver_blood_not_a_donor(self):
        liver=q("GTEX_V8_BULK_LIVER_RNA_MEDIAN","tacrolimus","hepatic_CYP3A_metabolism")
        blood=q("GTEX_V8_WHOLE_BLOOD_RNA_MEDIAN","tacrolimus","hepatic_CYP3A_metabolism")
        self.assertEqual(liver["mechanistic_relation_of_sample_to_process"],
                         "SOURCE_IN_TARGET_ORGAN_BUT_NOT_PROOF_IN_VIVO_FUNCTION")
        self.assertEqual(blood["mechanistic_relation_of_sample_to_process"],
                         "BLOOD_CELL_RNA_NOT_HEPATIC_OR_RENAL_TARGET_TISSUE")
        self.assertEqual(liver["sample_timing_interpretation"],
                         "REFERENCE_ATLAS_MEDIAN_NOT_INDIVIDUAL_BASELINE")
        self.assertFalse(liver["eligible_as_baseline_individual_PK_omic_proxy_under_source_contract"])
        self.assertFalse(blood["eligible_as_baseline_individual_PK_omic_proxy_under_source_contract"])

    def test_rotroff_real_plasma_baseline_vs_postdose_ABC(self):
        a=q("ROTROFF_PLASMA_GCTOF_A","metformin","renal_transporter_secretion",
            require_real_human_drug_PK=True)
        b=q("ROTROFF_PLASMA_GCTOF_B","metformin","renal_transporter_secretion")
        c=q("ROTROFF_PLASMA_GCTOF_C","metformin","renal_transporter_secretion")
        self.assertEqual(a["actual_collected_fraction"],"plasma")
        self.assertEqual(a["actual_anatomical_collection_organ"],"blood")
        self.assertTrue(a["eligible_as_baseline_individual_PK_omic_proxy_under_source_contract"])
        self.assertFalse(a["full_PGx_only_model_benchmark_possible_now"])
        for x in (b,c):
            self.assertFalse(x["eligible_as_baseline_individual_PK_omic_proxy_under_source_contract"])
            self.assertEqual(x["sample_timing_interpretation"],
                             "POSTDOSE_NOT_ELIGIBLE_AS_BASELINE_PREDICTOR")
        self.assertFalse(q("ROTROFF_PLASMA_GCTOF_B","metformin",
                           "renal_transporter_secretion",
                           intended_prediction="on_treatment_future_landmark")
                         ["eligible_as_baseline_individual_PK_omic_proxy_under_source_contract"])
        self.assertEqual(a["mechanistic_relation_of_sample_to_process"],
                         "SYSTEMIC_BIOFLUID_PROXY_NOT_TISSUE_SPECIFIC_ENZYME_FUNCTION")
        self.assertIn("BLOOD_OR_URINE_SAMPLE",a["potential_clinical_collection_feasibility"])
        self.assertFalse(a["assay_clinically_validated_for_dosing"])

    def test_Latham_original_pHh_protein_not_function_or_clinical_PK(self):
        prot=q("LATHAM_PHH_CYP_PROTEIN","imatinib","hepatic_CYP_metabolism")
        act=q("LATHAM_PHH_SHARED_PROBE","imatinib","hepatic_CYP_metabolism")
        self.assertEqual(prot["actual_collected_fraction"],"primary_hepatocytes")
        self.assertEqual(act["actual_collected_fraction"],"primary_hepatocytes")
        self.assertFalse(prot["source_is_direct_measurement_of_expected_target_organ_function"])
        self.assertTrue(act["source_is_direct_measurement_of_expected_target_organ_function"])
        self.assertFalse(act["original_paired_human_drug_PK_source_present"])
        self.assertFalse(act["eligible_as_baseline_individual_PK_omic_proxy_under_source_contract"])
        self.assertIn("INVASIVE_BIOPSY_OR_EX_VIVO",
                      act["potential_clinical_collection_feasibility"])

    def test_EATRIS_blood_cell_transcript_not_plasma_or_liver(self):
        rna=q("EATRIS_WHITE_BLOOD_CELLS_RNA","tacrolimus","hepatic_CYP3A_metabolism")
        p=q("EATRIS_PLASMA_PROTEIN","tacrolimus","hepatic_CYP3A_metabolism")
        m=q("EATRIS_PLASMA_METABOLITES","tacrolimus","hepatic_CYP3A_metabolism")
        self.assertEqual(rna["mechanistic_relation_of_sample_to_process"],
                         "BLOOD_CELL_RNA_NOT_HEPATIC_OR_RENAL_TARGET_TISSUE")
        self.assertEqual(p["mechanistic_relation_of_sample_to_process"],
                         "SYSTEMIC_BIOFLUID_PROXY_NOT_TISSUE_SPECIFIC_ENZYME_FUNCTION")
        self.assertEqual(m["actual_collected_fraction"],"plasma")
        self.assertFalse(any(x["original_paired_human_drug_PK_source_present"]
                             for x in (rna,p,m)))

    def test_tacrolimus_immune_PD_blood_cells_can_be_proxy_not_tissue_CYP(self):
        wbc=q("EATRIS_WHITE_BLOOD_CELLS_RNA","tacrolimus","immune_T_cell_PD")
        self.assertEqual(wbc["mechanistic_relation_of_sample_to_process"],
                         "IMMUNE_CELL_MIXTURE_POTENTIAL_PD_PROXY_NOT_TCELL_SPECIFIC_EFFECT")
        self.assertFalse(wbc["original_paired_human_drug_PK_source_present"])
        self.assertFalse(wbc["eligible_as_baseline_individual_PK_omic_proxy_under_source_contract"])
        liver=q("EATRIS_WHITE_BLOOD_CELLS_RNA","tacrolimus",
                "hepatic_CYP3A_metabolism")
        self.assertEqual(liver["mechanistic_relation_of_sample_to_process"],
                         "BLOOD_CELL_RNA_NOT_HEPATIC_OR_RENAL_TARGET_TISSUE")

    def test_HLC_liver_native_vs_liver_microsomes(self):
        a=q("HLC_LIVER_RNA","tacrolimus","hepatic_CYP3A_metabolism")
        b=q("HLC_LIVER_MICROSOME_CYP_PROBE","tacrolimus","hepatic_CYP3A_metabolism")
        self.assertEqual(a["actual_collected_fraction"],"bulk_liver_tissue")
        self.assertEqual(b["actual_collected_fraction"],"liver_microsomes")
        self.assertFalse(a["source_is_direct_measurement_of_expected_target_organ_function"])
        self.assertTrue(b["source_is_direct_measurement_of_expected_target_organ_function"])
        self.assertFalse(a["full_PGx_only_model_benchmark_possible_now"])

    def test_NCI60_cross_lab_lines_not_human_tacrolimus_biopsy(self):
        for assay in ("NCI60_RNA","NCI60_SWATH_PROTEOMICS","NCI60_METABOLITES"):
            x=q(assay,"doxorubicin","cancer_cellular_PD")
            self.assertEqual(x["actual_collected_fraction"],"cancer_cell_line")
            self.assertFalse(x["original_paired_human_drug_PK_source_present"])
            self.assertFalse(x["eligible_as_baseline_individual_PK_omic_proxy_under_source_contract"])
        x=q("NCI60_RNA","tacrolimus","hepatic_CYP3A_metabolism")
        self.assertEqual(x["mechanistic_relation_of_sample_to_process"],
                         "CANCER_CELL_LINE_ONLY_NOT_PATIENT_ORGAN")

    def test_engineered_cyp2c9_yeast_not_liver_patient(self):
        for assay in ("CYP2C9_DMS_VAMPSEQ","CYP2C9_DMS_CLICKSEQ",
                      "NUDT15_DMS_HEK293T_PROTEIN","NUDT15_DMS_HEK293T_TOLERANCE"):
            x=q(assay,"tacrolimus","hepatic_CYP3A_metabolism")
            self.assertEqual(x["mechanistic_relation_of_sample_to_process"],
                             "ENGINEERED_VARIANT_CONSTRUCT_NOT_REAL_PATIENT_ORGAN")
            self.assertFalse(x["original_paired_human_drug_PK_source_present"])
            self.assertFalse(x["patient_PK_improvement_vs_complete_genotype_proven"])

    def test_same_original_human_phh_protein_probe_arms_and_no_cross_cohort_join(self):
        model={
            "clinical_plus_genotype":[],
            "add_proteomics":["LATHAM_PHH_CYP_PROTEIN"],
            "add_probe_activity":["LATHAM_PHH_SHARED_PROBE"],
            "add_both":["LATHAM_PHH_CYP_PROTEIN","LATHAM_PHH_SHARED_PROBE"],
        }
        evidence=qualify_public_omic_model_arms(
            arm_to_original_sources=model,drug="imatinib",
            process="hepatic_CYP_metabolism",
            endpoint_class="IN_VITRO_DRUG_METABOLISM_NOT_CLINICAL_PK")
        self.assertEqual(evidence["n_distinct_original_measured_assays"],2)
        self.assertEqual(evidence["source_public_original_reference"],
                         "Latham_2025_PMC12990092_Table4")
        self.assertEqual(
            evidence["arm_original_collection_sites"]["add_both"][0]["actual_collected_fraction"],
            "primary_hepatocytes")
        bad={**model,"add_both":["LATHAM_PHH_CYP_PROTEIN","HLC_LIVER_RNA"]}
        with self.assertRaisesRegex(PhysiopharmacoSpecimenError,"different studies"):
            qualify_public_omic_model_arms(
                arm_to_original_sources=bad,drug="imatinib",
                process="hepatic_CYP_metabolism",
                endpoint_class="IN_VITRO_DRUG_METABOLISM_NOT_CLINICAL_PK")

    def test_metformin_predose_plasma_clinical_PK_allowed_only_as_proxy_no_PGx_claim(self):
        result=qualify_public_omic_model_arms(
            arm_to_original_sources={
                "clinical_plus_genotype":[],
                "add_metabolomics":["ROTROFF_PLASMA_GCTOF_A"],
            },drug="metformin",process="renal_transporter_secretion",
            endpoint_class="CLINICAL_PK")
        self.assertFalse(result["no_clinical_superiority_over_complete_genotype_claim"] is False)
        # The source does not have original genuine complete PGx: this
        # is specimen-method eligibility only, NOT benchmark readiness.
        self.assertFalse(
            result["arm_original_collection_sites"]["add_metabolomics"][0]
                  ["full_PGx_only_model_benchmark_possible_now"])
        for post in ("ROTROFF_PLASMA_GCTOF_B","ROTROFF_PLASMA_GCTOF_C"):
            with self.assertRaisesRegex(PhysiopharmacoSpecimenError,"post-dose"):
                qualify_public_omic_model_arms(
                    arm_to_original_sources={
                        "clinical_plus_genotype":[],"add_metabolomics":[post],
                    },drug="metformin",process="renal_transporter_secretion",
                    endpoint_class="CLINICAL_PK")
        with self.assertRaisesRegex(PhysiopharmacoSpecimenError,"No real same-original"):
            qualify_public_omic_model_arms(
                arm_to_original_sources={
                    "clinical_plus_genotype":[],
                    "add_RNA":["HLC_LIVER_RNA"],
                },drug="tacrolimus",process="hepatic_CYP3A_metabolism",
                endpoint_class="CLINICAL_PK")

    def test_organ_donor_versus_recipient_sample_genotype_origin(self):
        def geno(context,organ,sample,identity=True):
            return qualify_pharmacogenetic_sampling_origin(
                organ_transplant_status=context,target_function_organ=organ,
                original_DNA_sample=sample,original_alleles_covered_verified=True,
                donor_recipient_identity_verified=identity,
                hematopoietic_chimerism_qualified=True)
        # CYP hepatic metabolism after liver transplant: donor genotype
        # matches transplanted liver, recipient peripheral blood does not.
        donor=geno("liver_transplant","liver","donor_source_genotype")
        self.assertEqual(donor["status"],
                         "ORIGINAL_GENOTYPE_ORGAN_LINEAGE_ELIGIBLE_NOT_FUNCTION_ASSAY")
        blood=geno("liver_transplant","liver","recipient_predonation_blood")
        self.assertEqual(blood["status"],"ABSTAIN_TARGET_ORGAN_GERMLINE_UNQUALIFIED")
        # Intestinal first pass remains recipient-derived in liver recipient.
        intestine=geno("liver_transplant","intestine","recipient_predonation_blood")
        self.assertEqual(intestine["status"],
                         "ORIGINAL_GENOTYPE_ORGAN_LINEAGE_ELIGIBLE_NOT_FUNCTION_ASSAY")
        # A liver graft sample is not the genotype of kidney or intestine.
        incorrect=geno("liver_transplant","intestine","transplant_graft_tissue")
        self.assertEqual(incorrect["status"],"ABSTAIN_TARGET_ORGAN_GERMLINE_UNQUALIFIED")
        kidney=geno("kidney_transplant","kidney","transplant_graft_tissue")
        self.assertEqual(kidney["status"],
                         "ORIGINAL_GENOTYPE_ORGAN_LINEAGE_ELIGIBLE_NOT_FUNCTION_ASSAY")
        no_alleles=qualify_pharmacogenetic_sampling_origin(
            organ_transplant_status="liver_transplant",target_function_organ="liver",
            original_DNA_sample="donor_source_genotype",
            original_alleles_covered_verified=False,
            donor_recipient_identity_verified=True,
            hematopoietic_chimerism_qualified=True)
        self.assertEqual(no_alleles["status"],
                         "ABSTAIN_TARGET_ORGAN_GERMLINE_UNQUALIFIED")

    def test_hematopoietic_donor_blood_not_recipient_germline_after_allo_HSCT(self):
        a=qualify_pharmacogenetic_sampling_origin(
            organ_transplant_status="allogeneic_hsct",target_function_organ="liver",
            original_DNA_sample="current_peripheral_blood_leukocytes",
            original_alleles_covered_verified=True,
            donor_recipient_identity_verified=True,
            hematopoietic_chimerism_qualified=False)
        self.assertEqual(a["status"],"ABSTAIN_TARGET_ORGAN_GERMLINE_UNQUALIFIED")
        self.assertIn("MAY_BE_DONOR",a["source_origin_inference_limit"])
        b=qualify_pharmacogenetic_sampling_origin(
            organ_transplant_status="allogeneic_hsct",target_function_organ="liver",
            original_DNA_sample="recipient_native_nonhematopoietic_tissue",
            original_alleles_covered_verified=True,
            donor_recipient_identity_verified=True,
            hematopoietic_chimerism_qualified=True)
        self.assertEqual(b["status"],
                         "ORIGINAL_GENOTYPE_ORGAN_LINEAGE_ELIGIBLE_NOT_FUNCTION_ASSAY")
    def test_stop_unknown_noncatalogue_and_mismatched_process(self):
        with self.assertRaises(PhysiopharmacoSpecimenError):
            q("Unverified_RNA_CYP3A5_from_plasma","tacrolimus","hepatic_CYP3A_metabolism")
        with self.assertRaises(PhysiopharmacoSpecimenError):
            q("ROTROFF_PLASMA_GCTOF_A","metformin","hepatic_CYP3A_metabolism")
        with self.assertRaises(PhysiopharmacoSpecimenError):
            q("ROTROFF_PLASMA_GCTOF_A","metformin","renal_transporter_secretion",
              intended_prediction="predose_or_whatever")

if __name__=="__main__":
    unittest.main()
