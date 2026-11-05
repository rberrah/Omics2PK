"""No clinical data: adversarial, transcript-aware assay triage contracts."""
import unittest

from omics2pk.assay_aware_physio_triage import (
    METHODS, OmicTriageError, _classical_nmd_rule, plan_physio_omic_acquisition,
)

SOURCE="https://pmc.ncbi.nlm.nih.gov/articles/PMC10565868/"
NMD="https://pmc.ncbi.nlm.nih.gov/articles/PMC10045457/"


def a(name, method, gene="CYP3A5", tissue="liver", **kwargs):
    return {"assay_id":name,"method":method,"gene":gene,
            "tissue":tissue,"status":"proposed","collection_phase":"predose",**kwargs}


def variant(consequence="splice_disrupting", **kwargs):
    return {"gene":"CYP3A5","call_status":"complete",
            "source_variant_id":"rs776746",
            "source_url":SOURCE,"consequence":consequence,
            "transcript_id":"NM_000777.5",**kwargs}


def plan(assays, variants=(), drug="tacrolimus", endpoint="AUC"):
    return plan_physio_omic_acquisition(
        drug=drug, endpoint=endpoint, assays=assays, variants=variants)


def byname(result, name):
    return next(v for v in result["assay_candidates_by_interpretability_not_expected_model_performance"]
                if v["assay_id"]==name)


class PhysiopharmMethodTriageTests(unittest.TestCase):
    def test_early_stop_ejc_rule_is_plausibility_not_certainty(self):
        v=variant("premature_stop",premature_stop_nt=200,last_exon_junction_nt=350)
        result=plan([a("bulk","RNA_GENE_BULK"),
                     a("junction","RNA_SPLICE_JUNCTION",transcript_id="NM_000777.5"),
                     a("long","RNA_LONG_READ_ISOFORM",transcript_id="NM_000777.5")],
                    variants=[v])
        self.assertEqual(result["variant_source_consequence_interpretations"][0]
                         ["canonical_NMD_rule_not_function_proof"],
                         "CLASSICAL_EJC_NMD_PLAUSIBLE_NOT_CONFIRMED")
        self.assertEqual(byname(result,"bulk")["candidate_explanatory_role"],
                         "total_rna_proxy")
        self.assertIn("TOTAL_RNA_CANNOT_DISTINGUISH_FUNCTIONAL_FROM_ABERRANT_TRANSCRIPT",
                      byname(result,"bulk")["limitations"])
        self.assertEqual(byname(result,"junction")
                         ["candidate_explanatory_role"],"variant_or_isoform_resolution")
        self.assertLess(byname(result,"long")["heuristic_acquisition_priority_tier_not_predictive_rank"],
                        byname(result,"bulk")["heuristic_acquisition_priority_tier_not_predictive_rank"])
        self.assertTrue(result["no_quantitative_activity_or_effect_sizes_invented"])
        self.assertTrue(result["original_outcome_not_used_for_assay_prioritization"])
        self.assertFalse(result["no_PopPK_CYP_CL_or_dose_inferred"] is False)

    def test_last_exon_stop_and_near_boundary_not_called_nmd(self):
        for stop,last,want in (
            (1050,1000,"CLASSICAL_EJC_NMD_NOT_PREDICTED_LAST_EXON_OR_AFTER_FINAL_JUNCTION"),
            (950,1000,"CLASSICAL_EJC_NMD_NOT_PREDICTED_WITHIN_LAST_50NT"),
            (948,1000,"NMD_BOUNDARY_51_TO_54NT_INDETERMINATE"),
            (945,1000,"CLASSICAL_EJC_NMD_PLAUSIBLE_NOT_CONFIRMED"),
        ):
            with self.subTest(stop=stop,last=last):
                v=variant("premature_stop",premature_stop_nt=stop,last_exon_junction_nt=last)
                r=plan([a("rna","RNA_3PRIME_TAG")],[v])
                self.assertEqual(r["variant_source_consequence_interpretations"][0]
                                 ["canonical_NMD_rule_not_function_proof"],want)
                self.assertIn("3PRIME_TAG_MAY_PERSIST_DESPITE_UPSTREAM_STOP",
                              byname(r,"rna")["limitations"])

    def test_missing_nmd_geometry_is_unknown_not_always_zero_RNA(self):
        v=variant("frameshift",star_allele="*7",source_variant_id="rs41303343")
        r=plan([a("total","RNA_GENE_BULK")],[v])
        self.assertEqual(r["variant_source_consequence_interpretations"][0]
                         ["canonical_NMD_rule_not_function_proof"],
                         "UNKNOWN_NO_TRANSCRIPT_RESOLVED_STOP_AND_LAST_JUNCTION")
        self.assertIn("PREMATURE_TERMINATION_NMD_GEOMETRY_UNRESOLVED",
                      byname(r,"total")["limitations"])
        for vbad in (
            variant("frameshift",premature_stop_nt=200),
            variant("frameshift",premature_stop_nt="200",last_exon_junction_nt=400),
            variant("frameshift",premature_stop_nt=True,last_exon_junction_nt=400),
        ):
            with self.subTest(vbad=vbad), self.assertRaises(OmicTriageError):
                plan([a("rna","RNA_GENE_BULK")],[vbad])

    def test_cyp3a5_star3_star6_are_splice_defects_not_simple_stop_genotypes(self):
        for star,consequence in (("*3","splice_disrupting"),("*6","splice_disrupting"),
                                 ("*7","frameshift")):
            with self.subTest(star=star):
                r=plan([a("iso","RNA_LONG_READ_ISOFORM",
                          transcript_id="NM_000777.5")],
                       [variant(consequence,star_allele=star)])
                self.assertEqual(r["variant_source_consequence_interpretations"][0]
                                 ["canonical_NMD_rule_not_function_proof"],
                                 "UNKNOWN_NO_TRANSCRIPT_RESOLVED_STOP_AND_LAST_JUNCTION")
                self.assertFalse(r["no_false_CYP3A5_star_allele_from_gene_level_stop_call"] is False)
        with self.assertRaises(OmicTriageError):
            plan([a("rna","RNA_GENE_BULK")],
                 [variant("premature_stop",star_allele="*3")])
        with self.assertRaises(OmicTriageError):
            plan([a("rna","RNA_GENE_BULK")],
                 [variant("unknown",star_allele="*1")])

    def test_stop_allele_dosage_is_not_inferred_and_wildtype_can_contribute(self):
        r=plan([a("total","RNA_GENE_BULK")],
               [variant("premature_stop",allele_state="heterozygous")])
        self.assertEqual(r["variant_source_consequence_interpretations"][0]["allele_state"],
                         "heterozygous")
        self.assertIn("WILDTYPE_ALLELE_MAY_STILL_PRODUCE_FUNCTIONAL_TRANSCRIPT",
                      byname(r,"total")["limitations"])
        with self.assertRaises(OmicTriageError):
            plan([a("total","RNA_GENE_BULK")],
                 [variant("premature_stop",allele_state="triple_mutant")])

    def test_proteins_distinguish_truncated_fragment_from_intact_function(self):
        variants=[variant("premature_stop",premature_stop_nt=100,
                          last_exon_junction_nt=400)]
        r=plan([
            a("peptideN","PROTEIN_TARGETED_PEPTIDE",peptide_region="n_terminal"),
            a("peptideC","PROTEIN_TARGETED_PEPTIDE",peptide_region="c_terminal"),
            a("bottom","PROTEIN_BOTTOM_UP_SHARED_PEPTIDE"),
            a("top_unverified","PROTEIN_TOP_DOWN_INTACT",intact_identity_verified=False),
            a("top_verified","PROTEIN_TOP_DOWN_INTACT",intact_identity_verified=True,
              original_spectrum_url="https://proteomecentral.proteomexchange.org/"),
        ],variants)
        self.assertIn("N_TERMINAL_PEPTIDE_CAN_PERSIST_IN_TRUNCATED_INACTIVE_PROTEIN",
                      byname(r,"peptideN")["limitations"])
        self.assertIn("C_TERMINAL_PEPTIDE_MAY_HELP_DISTINGUISH_TRUNCATION_IF_COVERAGE_VERIFIED",
                      byname(r,"peptideC")["evidence_observations"])
        self.assertIn("SHARED_PEPTIDES_MAY_INCLUDE_TRUNCATED_OR_PARALOGOUS_PROTEIN",
                      byname(r,"bottom")["limitations"])
        self.assertIn("INTACT_PROTEOFORM_ABUNDANCE_NOT_ENZYME_FUNCTION",
                      byname(r,"top_verified")["limitations"])
        self.assertEqual(byname(r,"top_verified")["candidate_explanatory_role"],
                         "variant_or_isoform_resolution")
        self.assertEqual(byname(r,"top_unverified")["candidate_explanatory_role"],
                         "protein_or_isoform_abundance")
        with self.assertRaises(OmicTriageError):
            plan([a("false_top","PROTEIN_TOP_DOWN_INTACT",intact_identity_verified=True)])

    def test_selective_vs_shared_functional_activity_and_laboratory_doxorubicin_warning(self):
        r=plan([
            a("specific","ACTIVITY_CYP3A5_SELECTIVE_PROBE",assay_validated=True,
              selectivity_evidence_url="https://doi.org/10.1124/dmd.113.054726",
              probe_assay_code="T5_N_OXIDATION"),
            a("shared","ACTIVITY_CYP3A_SHARED_PROBE",assay_validated=True,
              probe_assay_code="MIDAZOLAM_1_HYDROXYLATION",
              probe_evidence_url="https://doi.org/10.1124/dmd.119.088823"),
            a("unqualified","ACTIVITY_CYP3A5_SELECTIVE_PROBE",assay_validated=False,
              selectivity_evidence_url="https://doi.org/10.1124/dmd.113.054726",
              probe_assay_code="T5_N_OXIDATION"),
        ])
        self.assertEqual(byname(r,"specific")["candidate_explanatory_role"],
                         "direct_measured_function")
        self.assertIn("SELECTIVE_CYP3A5_ACTIVITY_PROXY_NOT_ACTUAL_TACROLIMUS_CLEARANCE",
                      byname(r,"specific")["evidence_observations"])
        self.assertIn("SHARED_CYP3A4_CYP3A5_ACTIVITY_NOT_ISOFORM_SPECIFIC",
                      byname(r,"shared")["limitations"])
        self.assertEqual(byname(r,"unqualified")["candidate_explanatory_role"],
                         "nonspecific_signal")
        self.assertEqual(r["prospective_fold_safe_evaluation_plan"]
                         ["required_clinical_covariates"][4],"hematocrit")

    def test_metabolite_is_not_isoform_selective_or_identity_confirmed_from_ms1(self):
        r=plan([
            {"assay_id":"ms","method":"METABOLITE_MS1_UNTARGETED",
             "tissue":"plasma","candidate_genes":["CYP3A4","CYP3A5"],
             "chemical_identity_grade":"source_name_only",
             "status":"proposed","collection_phase":"predose"},
            {"assay_id":"met","method":"METABOLITE_MSMS_TARGETED",
             "tissue":"plasma","candidate_genes":["CYP3A5"],
             "chemical_identity_grade":"confirmed_standard",
             "status":"proposed","collection_phase":"predose"},
        ])
        self.assertIn("METABOLITE_IS_NOT_SPECIFIC_TO_ONE_CYP_OR_TRANSPORTER",
                      byname(r,"ms")["limitations"])
        self.assertIn("SMALL_MOLECULE_IDENTITY_UNVERIFIED",
                      byname(r,"ms")["limitations"])
        self.assertIn("CIRCULATING_PROXY_TISSUE_NOT_DIRECT_ORGAN_ACTIVITY",
                      byname(r,"met")["limitations"])
        self.assertEqual(byname(r,"ms")["candidate_explanatory_role"],"nonspecific_signal")
        with self.assertRaises(OmicTriageError):
            plan([{"assay_id":"false","method":"METABOLITE_MS1_UNTARGETED",
                   "tissue":"plasma","candidate_genes":["CYP3A5"],
                   "chemical_identity_grade":"confirmed_standard",
                   "status":"proposed","collection_phase":"predose"}])

    def test_renal_cl_mechanism_avoids_tacrolimus_Pgp_and_hepatic_oct1(self):
        good=[a("oct2","PROTEIN_TARGETED_PEPTIDE",gene="SLC22A2",
                tissue="kidney",peptide_region="c_terminal"),
              a("oct1","RNA_GENE_BULK",gene="SLC22A1",tissue="liver")]
        r=plan(good,drug="metformin",endpoint="CLr")
        self.assertEqual(byname(r,"oct2")["status"],
                         "PROSPECTIVE_MEASUREMENT_NOT_YET_ELIGIBLE")
        self.assertEqual(byname(r,"oct1")["status"],
                         "ABSTAIN_WRONG_DRUG_PROCESS_FOR_ENDPOINT")
        self.assertFalse(byname(r,"oct1")["eligible_for_outcome_ablation"])
        self.assertEqual(r["drug"],"metformin")
        self.assertIn("original_renal_function",r["prospective_fold_safe_evaluation_plan"]
                      ["required_clinical_covariates"])
        with self.assertRaises(OmicTriageError):
            plan([a("pgp","RNA_GENE_BULK",gene="ABCB1")],
                 drug="metformin",endpoint="CLr",variants=[
                     {"gene":"ABCB1","consequence":"unknown","call_status":"missing"}])

    def test_original_measured_not_falsely_source_authored_and_postdose_forbidden(self):
        r=plan([
            a("real","RNA_SPLICE_JUNCTION",transcript_id="NM_000777.5",
              status="original_measured",source_url=SOURCE,
              original_feature_id="junction_1"),
            a("postdose","RNA_GENE_BULK",status="original_measured",
              source_url=SOURCE,original_feature_id="gene_count",
              collection_phase="postdose"),
        ])
        self.assertTrue(byname(r,"real")["eligible_for_outcome_ablation"])
        self.assertFalse(byname(r,"postdose")["eligible_for_outcome_ablation"])
        self.assertEqual(byname(r,"postdose")["status"],
                         "ABSTAIN_PREDOSE_LEAKAGE_OR_TIMING_UNRESOLVED")
        self.assertFalse(r["prospective_fold_safe_evaluation_plan"]
                         ["actual_same_subject_modalities_intersection_verified"])
        self.assertEqual(r["prospective_fold_safe_evaluation_plan"]
                         ["candidate_arms_not_fitted_models"][0]["modalities"],[])
        with self.assertRaises(OmicTriageError):
            plan([a("not_actually_measured","RNA_GENE_BULK",
                    status="original_measured")])

    def test_gene_source_timing_and_duplicate_inputs_fail_closed(self):
        with self.assertRaises(OmicTriageError):
            plan([a("rna","RNA_GENE_BULK")], drug="unknown")
        with self.assertRaises(OmicTriageError):
            plan([a("rna","RNA_GENE_BULK")], endpoint="GI50")
        with self.assertRaises(OmicTriageError):
            plan([a("rna","RNA_GENE_BULK"),a("rna","RNA_GENE_BULK")])
        with self.assertRaises(OmicTriageError):
            plan([a("wrong","RNA_GENE_BULK",gene="CYP3A5",tissue="kidney")],
                 [variant("unknown",source_variant_id="rs776746",source_url="http://bad")])
        with self.assertRaises(OmicTriageError):
            plan([a("known","RNA_GENE_BULK")],[variant("splice_disrupting",
                                                    transcript_id="BAD.1")])
        with self.assertRaises(OmicTriageError):
            plan([a("none","RNA_GENE_BULK")],
                 [{"gene":"CYP3A5","consequence":"premature_stop",
                   "call_status":"missing","source_variant_id":"rs776746"}])

    def test_assay_transcript_mismatch_is_explicit_not_splice_evidence(self):
        r=plan([a("rna","RNA_LONG_READ_ISOFORM",transcript_id="NM_000777.2")],
               [variant("splice_disrupting",star_allele="*3")])
        self.assertIn("ASSAY_TRANSCRIPT_MAY_NOT_MATCH_VARIANT_TRANSCRIPT",
                      byname(r,"rna")["limitations"])
        self.assertFalse(byname(r,"rna")["mechanism_confirmed_for_subject"])
        self.assertFalse(r["prospective_fold_safe_evaluation_plan"]
                         ["external_validation_required"] is False)


if __name__=="__main__":
    unittest.main()
