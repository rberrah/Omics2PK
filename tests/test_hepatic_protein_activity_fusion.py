"""Synthetic-only original two-compound hepatic protein+activity fusion tests."""
from __future__ import annotations

import unittest
from omics2pk.hepatic_protein_activity_fusion import (
    run_fusion_on_verified_same_donors,from_original_source_rows
)
from omics2pk.matched_tki_metabolism import MatchedMetabolismError
from test_matched_tki_metabolism import fixtures


class SameDonorProteinFunctionalFusionTests(unittest.TestCase):
    def test_comparison_identical_donor_folds_all_arms_and_two_hepatic_drugs(self):
        a,b=fixtures()
        # Synthetic values do not reproduce original paper Spearman by design.
        from omics2pk.matched_tki_metabolism import qualify_pair
        joined,_=qualify_pair(a,b,verify_reference=False)
        r=run_fusion_on_verified_same_donors(joined)
        self.assertEqual(r["original_biological_donors_total"],21)
        self.assertFalse(r["can_test_RNA_only_vs_RNA_plus_proteins_on_same_donors"])
        self.assertTrue(r["not_a_multiomic_comparison"])
        sizes=sorted(x["n_original_donors_same_for_all_arms"] for x in r["outcomes"].values())
        self.assertEqual(sizes,[17,21,21])
        self.assertNotIn("S00",str(r))
        for out in r["outcomes"].values():
            self.assertEqual(set(out["MAE_by_fixed_modality_arm_original_units"]),
                             {"cohort_mean","CYP_proteins","CYP_probe_activities",
                              "CYP_proteins_plus_activities"})
            self.assertTrue(out["original_donor_level_losses_and_predictions_not_emitted"])

    def test_duplicate_source_human_donor_rejected(self):
        a,b=fixtures()
        from omics2pk.matched_tki_metabolism import qualify_pair
        paired,_=qualify_pair(a,b,verify_reference=False)
        paired[1]=(paired[0][0],*paired[1][1:])
        with self.assertRaisesRegex(MatchedMetabolismError,"21 original"):
            run_fusion_on_verified_same_donors(paired)

    def test_imatinib_unquantifiable_source_outcomes_are_never_zeroes(self):
        a,b=fixtures()
        from omics2pk.matched_tki_metabolism import qualify_pair
        paired,_=qualify_pair(a,b,verify_reference=False)
        bad=list(paired)
        bad[1]=(bad[1][0],bad[1][1],0.0,bad[1][3],bad[1][4],bad[1][5])
        with self.assertRaisesRegex(MatchedMetabolismError,"Four original"):
            run_fusion_on_verified_same_donors(bad)


if __name__=="__main__":
    unittest.main()
