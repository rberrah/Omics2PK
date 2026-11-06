"""Contract for mechanism-first assays; pure stdlib and no protected source data."""
from __future__ import annotations
import json
from pathlib import Path
import unittest

PRIOR=Path(__file__).resolve().parents[1]/"research"/"cyp3a_mechanism_first_assay_priors_20261106.json"


class MechanismFirstScientificSourceContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=json.loads(PRIOR.read_text(encoding="utf-8"))
        cls.hypotheses={v["id"]:v for v in cls.data["hypotheses"]}

    def test_no_false_clinical_PK_mechanism_claim(self):
        self.assertEqual(
            self.data["status"],"RESEARCH_ONLY_MECHANISM_FIRST_PRIORS_NOT_TRAINED_MODEL")
        self.assertEqual(len(self.hypotheses),len(self.data["hypotheses"]))
        for item in self.hypotheses.values():
            self.assertIs(item["clinical_pk_direct"],False)
            self.assertTrue(item["preferred_specimens"])
            self.assertTrue(item["required_validation"])
            self.assertTrue(item["source_urls"])
            self.assertTrue(all(url.startswith("https://") for url in item["source_urls"]))

    def test_cyp3a5_star3_splice_never_missense_or_blood_RNA_proxy(self):
        x=self.hypotheses["cyp3a5_star3_splice"]
        self.assertEqual(x["gene"],"CYP3A5")
        self.assertIn("rs776746",x["variants"])
        self.assertIn("splice",x["variant_effect"])
        self.assertNotIn("missense",x["variant_effect"])
        self.assertNotIn("plasma",x["preferred_specimens"])
        self.assertIn("whole_blood_bulk_CYP3A5_RNA",x["not_proved_by"])

    def test_source_2026_GHCA_is_cyp3a4_not_claimed_specific_cyp3a5(self):
        x=self.hypotheses["cyp3a4_GCDCA_to_GHCA"]
        self.assertEqual(x["gene"],"CYP3A4")
        self.assertIn("glycohyocholate_GHCA",x["downstream"])
        self.assertIn("plasma",x["preferred_specimens"])
        self.assertIn("CYP3A5_star3_genotype_solely",x["not_proved_by"])

    def test_gnps_unidentified_ions_are_not_clinical_biomarkers(self):
        x=self.hypotheses["Jarmusch_unknown_selected_ions"]
        self.assertEqual(x["type"],"hypothesis_generation_only")
        self.assertIsNone(x["gene"])
        self.assertEqual(
            set(x["source_ids"]),{"369","3018","4321","6171","936"})
        self.assertIn("signal over batch-matched extraction/analytical blanks",
                      x["required_validation"])
        self.assertTrue(
            self.data["qualification_policy"]["source_ms1_scan_id_is_not_compound_identity"])


if __name__=="__main__":
    unittest.main()
