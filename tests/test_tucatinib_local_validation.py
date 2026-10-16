"""Synthetic-only HLM protein measured-clearance benchmark tests."""
import unittest
from omics2pk.tucatinib_local_validation import qualify,score_donors,TucatinibEvidenceError,_spearman,verify_published_table_correlations

def fake_22():
    rows=[]
    for i in range(21):
        x=3+i*0.4
        rows.append(dict(donor_id=f"FAKE_{i:02}",row_kind="individual",
          source_doi="10.1016/j.dmd.2025.100061",
          tucatinib_cl_app_ul_min_mg=str(2*x+0.25*i),
          cyp2c8_protein_pmol_mg=str(x),
          cyp3a_protein_pmol_mg=str(2*x+0.5),
          cyp2c8_activity_pmol_min_mg=str(30+i),
          cyp3a_activity_pmol_min_mg=str(60+3*i)))
    rows.append(dict(donor_id="FAKE_POOLED",row_kind="pool",source_doi="10.1016/j.dmd.2025.100061",
      tucatinib_cl_app_ul_min_mg="500000",cyp2c8_protein_pmol_mg="1",cyp3a_protein_pmol_mg="1",
      cyp2c8_activity_pmol_min_mg="1",cyp3a_activity_pmol_min_mg="1"))
    return rows

class TucatinibSyntheticTests(unittest.TestCase):
    def test_strictly_21_individual_donors_and_no_pool_pseudoreplication(self):
        result=score_donors(qualify(fake_22()))
        self.assertEqual(result["n_unique_human_donors"],21)
        self.assertTrue(result["pooled_hlm_excluded"])
        self.assertLess(result["protein_mae_ul_min_mg"],result["baseline_mae_ul_min_mg"])
        self.assertNotIn("FAKE_",repr(result))

    def test_duplicate_donors_rejected(self):
        bad=fake_22()
        bad[1]["donor_id"]="FAKE_00"
        with self.assertRaisesRegex(TucatinibEvidenceError,"unique"):
            qualify(bad)

    def test_pooled_lot_required_not_22nd_donor(self):
        bad=fake_22()
        bad[-1]["row_kind"]="individual"
        with self.assertRaisesRegex(TucatinibEvidenceError,"pooled"):
            qualify(bad)

    def test_missing_units_and_doi_rejected(self):
        bad=fake_22()
        del bad[0]["cyp3a_protein_pmol_mg"]
        with self.assertRaisesRegex(TucatinibEvidenceError,"column"):
            qualify(bad)
        bad=fake_22()
        bad[0]["source_doi"]="wrong"
        with self.assertRaisesRegex(TucatinibEvidenceError,"DOI"):
            qualify(bad)

    def test_average_rank_spearman_with_ties(self):
        self.assertAlmostEqual(_spearman([1, 2, 2, 4], [2, 3, 3, 5]), 1.0)
        self.assertAlmostEqual(_spearman([1, 2, 3, 4], [4, 3, 2, 1]), -1.0)

    def test_source_crosscheck_rejects_synthetic_biomarker_relationship(self):
        measured = qualify(fake_22())
        with self.assertRaisesRegex(TucatinibEvidenceError, "published CYP-clearance correlations mismatch"):
            verify_published_table_correlations(measured)

    def test_negative_or_imputed_clearance_not_accepted(self):
        for value in ("","-2","nan","0"):
            bad=fake_22()
            bad[2]["tucatinib_cl_app_ul_min_mg"]=value
            with self.assertRaises(TucatinibEvidenceError):
                qualify(bad)

if __name__=="__main__":
    unittest.main()
