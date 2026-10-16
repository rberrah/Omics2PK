"""Synthetic source-gated two-drug same-HLM-donor validation tests."""
import unittest
from omics2pk.matched_tki_metabolism import (
    MatchedMetabolismError, validate_pair,
)


def fixtures():
    a, b = [], []
    missing = {1, 4, 7, 10}
    for i in range(21):
        protein2 = round(4 + 1.1*i, 2)
        protein3 = round(12 + 2.3*i, 2)
        activity2 = round(5 + 0.9*i, 2)
        activity3 = round(15 + 1.8*i, 2)
        donor = "S%02d" % i
        a.append({
            "donor_id":donor, "row_kind":"individual",
            "source_doi":"10.1016/j.dmd.2025.100061",
            "tucatinib_cl_app_ul_min_mg":str(2+i*0.4),
            "cyp2c8_protein_pmol_mg":str(protein2),
            "cyp3a_protein_pmol_mg":str(protein3),
            "cyp2c8_activity_pmol_min_mg":str(activity2),
            "cyp3a_activity_pmol_min_mg":str(activity3),
        })
        b.append({
            "donor_id":donor,"row_kind":"individual",
            "source_doi":"10.1016/j.dmd.2025.100196",
            "imatinib_cl_app_ul_min_mg":"" if i in missing else str(3.0+i*1.7),
            "ndesmethyl_imatinib_peak_area_ratio":str(0.05+0.014*i),
            "cyp2c8_protein_pmol_mg":str(protein2),
            "cyp3a_protein_pmol_mg":str(protein3),
            "cyp2c8_activity_pmol_min_mg":str(activity2),
            "cyp3a_activity_pmol_min_mg":str(activity3),
        })
    pool_t={**a[0],"donor_id":"pooled","row_kind":"pool"}
    pool_i={**b[0],"donor_id":"pooled","row_kind":"pool"}
    return a+[pool_t], b+[pool_i]


class SharedDonorBenchmarkTests(unittest.TestCase):
    def test_four_nonquantifiable_imatinib_clearances_are_not_zeroes(self):
        tuc, ima=fixtures()
        report=validate_pair(tuc, ima,verify_reference=False)
        self.assertEqual(report["original_donors"],21)
        self.assertEqual(report["imatinib_missing_clearance_not_imputed"],4)
        self.assertEqual(report["imatinib_observed_cl_ul_min_mg"]["n_independent_biological_donors"],17)
        self.assertEqual(report["ndesmethyl_imatinib_peak_area_ratio"]["n_independent_biological_donors"],21)
        self.assertTrue(report["no_cross_study_independence"])

    def test_source_cyp_mismatch_rejected(self):
        t,i=fixtures()
        i[1]["cyp3a_protein_pmol_mg"]="99999"
        with self.assertRaisesRegex(MatchedMetabolismError,"phenotypes disagree"):
            validate_pair(t,i,verify_reference=False)

    def test_fake_join_by_row_position_rejected(self):
        t,i=fixtures()
        i[2]["donor_id"]="EXTERNAL"
        with self.assertRaisesRegex(MatchedMetabolismError,"Donor IDs differ"):
            validate_pair(t,i,verify_reference=False)

    def test_imatinib_missing_cannot_be_imputed_with_zero(self):
        t,i=fixtures()
        i[1]["imatinib_cl_app_ul_min_mg"]="0"
        with self.assertRaisesRegex(MatchedMetabolismError,"must be positive"):
            validate_pair(t,i,verify_reference=False)

    def test_pool_cannot_be_counted_as_a_new_individual(self):
        t,i=fixtures()
        i[-1]["row_kind"]="individual"
        with self.assertRaises(MatchedMetabolismError):
            validate_pair(t,i,verify_reference=False)


if __name__=="__main__":
    unittest.main()
