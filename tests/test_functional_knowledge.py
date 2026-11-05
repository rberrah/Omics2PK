"""Synthetic fail-closed biological evidence tests; not a PK validation study."""
import unittest

from omics2pk.functional_knowledge import (
    CYP3A5_ALLELES, FunctionalEvidenceError, assess_cyp3a5_evidence,
)


class CYP3A5FunctionalEvidenceTests(unittest.TestCase):
    def test_cpic_function_and_allele_mechanisms(self):
        self.assertEqual(CYP3A5_ALLELES["*1"]["function"], "normal")
        for allele in ("*3", "*6", "*7"):
            self.assertEqual(CYP3A5_ALLELES[allele]["function"], "no_function")
        self.assertIn("insertion", CYP3A5_ALLELES["*7"]["mechanism"])
        self.assertIn("synonymous", CYP3A5_ALLELES["*6"]["mechanism"])

    def test_expected_diplotype_classes_without_percent_loss(self):
        for diplotype, expected in (("*1/*1", "expresser"),
                                   ("*1/*3", "intermediate_expresser"),
                                   ("*3/*3", "nonexpresser"),
                                   ("*6/*7", "nonexpresser")):
            with self.subTest(diplotype=diplotype):
                result = assess_cyp3a5_evidence(diplotype=diplotype, call_status="complete")
                self.assertEqual(result["diplotype_classification"], expected)
                self.assertIsNone(result["percent_activity_loss"])
                self.assertIsNone(result["estimated_enzyme_activity"])
                self.assertFalse(result["in_vivo_pk_gain_demonstrated"])

    def test_total_transcript_detected_in_nonexpresser_is_not_a_contradiction(self):
        out = assess_cyp3a5_evidence(
            diplotype="*3/*3", call_status="complete",
            observations=[{"omic": "transcriptomics", "target": "total_gene",
                           "signal": "detected", "tissue": "liver"}])
        self.assertEqual(out["status"], "INTERPRETATION_ONLY")
        self.assertIn("TOTAL_RNA_NOT_FUNCTIONAL_TRANSCRIPT", out["evidence_warnings"])
        self.assertNotIn("GENOTYPE_FUNCTIONAL_ACTIVITY_DISCORDANCE_REVIEW", out["evidence_warnings"])

    def test_full_length_transcript_is_still_not_function(self):
        out = assess_cyp3a5_evidence(
            diplotype="*3/*3", call_status="complete",
            observations=[{"omic": "transcriptomics", "target": "canonical_splice_junction",
                           "signal": "detected", "tissue": "liver"}])
        self.assertEqual(out["status"], "INTERPRETATION_ONLY")
        self.assertIn("SPLICE_JUNCTION_NOT_ENZYME_ACTIVITY", out["evidence_warnings"])

    def test_pan_cyp3a_protein_is_not_specific(self):
        out = assess_cyp3a5_evidence(
            diplotype="*3/*3", call_status="complete",
            observations=[{"omic": "proteomics", "target": "pan_cyp3a",
                           "signal": "detected", "tissue": "liver"}])
        self.assertIn("PROTEIN_ASSAY_NOT_FULL_LENGTH_SPECIFIC", out["evidence_warnings"])
        self.assertEqual(out["status"], "INTERPRETATION_ONLY")

    def test_metabolite_concentration_with_two_routes_not_specific(self):
        out = assess_cyp3a5_evidence(
            diplotype="*1/*3", call_status="complete",
            observations=[{"omic": "metabolomics", "target": "metabolite_abundance",
                           "signal": "detected", "candidate_enzymes": ["CYP3A4", "CYP3A5"],
                           "tissue": "liver"}])
        self.assertIn("METABOLITE_PATHWAY_NOT_UNIQUE", out["evidence_warnings"])
        self.assertIn("METABOLITE_ABUNDANCE_NOT_FLUX", out["evidence_warnings"])
        self.assertIsNone(out["estimated_enzyme_activity"])

    def test_tracer_flux_still_not_proof_of_enzyme(self):
        out = assess_cyp3a5_evidence(
            diplotype="*3/*3", call_status="complete",
            observations=[{"omic": "metabolomics", "target": "tracer_flux",
                           "signal": "detected", "candidate_enzymes": ["CYP3A5"],
                           "tissue": "liver"}])
        self.assertIn("TRACER_FLUX_NOT_AUTOMATICALLY_ENZYME_SPECIFIC", out["evidence_warnings"])

    def test_validated_high_specific_activity_is_review_not_override(self):
        out = assess_cyp3a5_evidence(
            diplotype="*3/*3", call_status="complete",
            observations=[{"omic": "activity", "target": "cyp3a5_specific",
                           "signal": "above_validated_reference", "assay_validated": True,
                           "tissue": "liver"}])
        self.assertEqual(out["status"], "REVIEW_FUNCTIONAL_DISCORDANCE")
        self.assertEqual(out["diplotype_classification"], "nonexpresser")
        self.assertIsNone(out["estimated_enzyme_activity"])

    def test_mere_residual_activity_is_not_contradiction(self):
        out = assess_cyp3a5_evidence(
            diplotype="*3/*3", call_status="complete",
            observations=[{"omic": "activity", "target": "cyp3a5_specific",
                           "signal": "detected", "assay_validated": True,
                           "tissue": "liver"}])
        self.assertEqual(out["status"], "INTERPRETATION_ONLY")

    def test_unvalidated_activity_not_overridden(self):
        out = assess_cyp3a5_evidence(
            diplotype="*3/*3", call_status="complete",
            observations=[{"omic": "activity", "target": "cyp3a5_specific",
                           "signal": "above_validated_reference", "assay_validated": False,
                           "tissue": "liver"}])
        self.assertEqual(out["status"], "INTERPRETATION_ONLY")
        self.assertIn("ACTIVITY_ASSAY_UNVALIDATED", out["evidence_warnings"])

    def test_allele_coverage_fail_closed(self):
        for diplotype, status, warning in [
            ("*1/*1", "partial", "PARTIAL_ALLELE_COVERAGE"),
            ("*1/*8", "complete", "UNRESOLVED_ALLELE_FUNCTION"),
            ("*3/*3", "partial", "PARTIAL_ALLELE_COVERAGE")]:
            with self.subTest(diplotype=diplotype, status=status):
                out = assess_cyp3a5_evidence(diplotype=diplotype, call_status=status)
                self.assertEqual(out["status"], "ABSTAIN_GENOTYPE_UNQUALIFIED")
                self.assertIsNone(out["diplotype_classification"])
                self.assertIn(warning, out["evidence_warnings"])
        self.assertEqual(assess_cyp3a5_evidence()["status"], "ABSTAIN_GENOTYPE_UNQUALIFIED")

    def test_unknown_or_invalid_evidence_fails_closed(self):
        bad = [
            {"omic": "transcriptomics", "target": "protein", "signal": "detected"},
            {"omic": "metabolomics", "target": "metabolite_abundance", "signal": "detected"},
            {"omic": "activity", "target": "cyp3a5_specific", "signal": "detected"},
            {"omic": "transcriptomics", "target": "total_gene", "signal": "above_validated_reference"},
        ]
        for obs in bad:
            with self.subTest(obs=obs), self.assertRaises(FunctionalEvidenceError):
                assess_cyp3a5_evidence(diplotype="*3/*3", call_status="complete", observations=[obs])
        with self.assertRaises(FunctionalEvidenceError):
            assess_cyp3a5_evidence(diplotype="*1/*1", call_status="missing")
        with self.assertRaises(FunctionalEvidenceError):
            assess_cyp3a5_evidence(diplotype="*1", call_status="complete")

    def test_nonhepatic_sample_flagged(self):
        out = assess_cyp3a5_evidence(
            diplotype="*1/*1", call_status="complete",
            observations=[{"omic": "transcriptomics", "target": "total_gene",
                           "signal": "detected", "tissue": "whole_blood"}])
        self.assertIn("NON_TARGET_TISSUE_EXTRAPOLATION", out["evidence_warnings"])


if __name__ == "__main__":
    unittest.main()
