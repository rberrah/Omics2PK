"""Source-verified assay taxonomy, synthetic software tests only."""
import unittest

from omics2pk.cyp3a_activity_phenotype import (
    ASSAY_ROLES, CYP3AAssayError, classify_cyp3a_assay,
)
from omics2pk.functional_activity_benchmark import (
    FunctionalBenchmarkError, benchmark_cyp3a5_activity,
)
from test_functional_activity_benchmark import fixture


class CYP3AProbeQualificationTests(unittest.TestCase):
    def test_T5_is_selective_but_not_enzyme_exclusive(self):
        out=classify_cyp3a_assay(
            assay_code="T5_N_OXIDATION",assay_validated=True,
            source_doi="10.1124/dmd.113.054726")
        self.assertTrue(out["CYP3A5_selective_activity_endpoint_eligible"])
        self.assertIn("not fully",out["interpretation_note"].lower())
        self.assertFalse(out["independent_baseline_molecular_measurement_possible"])

    def test_midazolam_probe_and_testosterone_are_shared_cyp3a(self):
        for key in ("MIDAZOLAM_1_HYDROXYLATION","TESTOSTERONE_6B_HYDROXYLATION"):
            with self.subTest(key=key):
                out=classify_cyp3a_assay(
                    assay_code=key, assay_validated=True,
                    source_doi=ASSAY_ROLES[key]["reference"])
                self.assertEqual(out["assay_category"],"SHARED_CYP3A_PROBE")
                self.assertFalse(out["CYP3A5_selective_activity_endpoint_eligible"])

    def test_endogenous_4beta_HC_not_specific_and_conditionally_independent(self):
        out=classify_cyp3a_assay(
            assay_code="FOUR_BETA_HYDROXYCHOLESTEROL_ABUNDANCE",
            assay_validated=True,source_doi="10.1002/jcph.1391")
        self.assertTrue(out["independent_baseline_molecular_measurement_possible"])
        self.assertFalse(out["CYP3A5_selective_activity_endpoint_eligible"])

    def test_post_drug_metabolite_outcomes_not_predose_metabolomics(self):
        for key in ("LAPATINIB_M1_FORMATION","VINORELBINE_PARENT_DEPLETION"):
            with self.subTest(key=key):
                out=classify_cyp3a_assay(
                    assay_code=key,assay_validated=True,
                    source_doi=ASSAY_ROLES[key]["reference"])
                self.assertFalse(out["independent_baseline_molecular_measurement_possible"])
                self.assertFalse(out["CYP3A5_selective_activity_endpoint_eligible"])

    def test_unknown_unvalidated_or_wrong_reference_rejected(self):
        bad=[("unknown",True,"10.1124/dmd.113.054726"),
             ("T5_N_OXIDATION",False,"10.1124/dmd.113.054726"),
             ("T5_N_OXIDATION",True,"10.1124/dmd.119.088823")]
        for code,validated,doi in bad:
            with self.subTest(code=code,doi=doi),self.assertRaises(CYP3AAssayError):
                classify_cyp3a_assay(
                    assay_code=code,assay_validated=validated,source_doi=doi)

    def test_original_activity_requires_source_attested_assay(self):
        d=fixture()
        d["source_kind"]="original_source_declared"
        with self.assertRaisesRegex(FunctionalBenchmarkError,"assay evidence invalid"):
            benchmark_cyp3a5_activity(**d)
        d["activity_assay_code"]="MIDAZOLAM_1_HYDROXYLATION"
        d["activity_assay_evidence_doi"]="10.1124/dmd.119.088823"
        with self.assertRaisesRegex(FunctionalBenchmarkError,"Shared CYP3A"):
            benchmark_cyp3a5_activity(**d)

    def test_t5_typed_synthetic_test_still_not_claimed_clinical(self):
        d=fixture()
        d["activity_assay_code"]="T5_N_OXIDATION"
        d["activity_assay_evidence_doi"]="10.1124/dmd.113.054726"
        out=benchmark_cyp3a5_activity(**d)
        self.assertEqual(out["source_assay_identity"],"T5_N_OXIDATION")
        self.assertIn("SYNTHETIC",out["status"])
        self.assertFalse(out["in_vivo_PK_prediction_demonstrated"])


if __name__=="__main__":
    unittest.main()
