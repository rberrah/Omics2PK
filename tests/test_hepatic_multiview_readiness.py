import unittest

from omics2pk.hepatic_multiview_readiness import (
    HepaticMultiviewReadinessError,
    assess_same_donor_hepatic_multiview,
    qualify_candidate_summary,
)


class HepaticMultiviewReadinessTests(unittest.TestCase):
    def test_ready_requires_identical_complete_case_donors(self):
        donors=["D1","D2","D3"]
        out=assess_same_donor_hepatic_multiview(
            source_id="AUTHENTIC_SOURCE",
            rna_donors=donors,
            protein_donors=donors,
            outcome_donors=donors,
            outcome_kind="drug_specific_metabolic_activity",
            source_rows_acquired=True,
        )
        self.assertTrue(out.fit_allowed)
        self.assertEqual(out.status,"READY_SAME_DONOR_RNA_PROTEIN_DRUG_OUTCOME")
        self.assertEqual(out.n_complete_case,3)

    def test_cross_study_or_partial_overlap_abstains(self):
        out=assess_same_donor_hepatic_multiview(
            source_id="DO_NOT_FUSE_STUDIES",
            rna_donors=["EPA1","EPA2"],
            protein_donors=["TKI1","TKI2"],
            outcome_donors=["TKI1","TKI2"],
            outcome_kind="drug_specific_clint",
            source_rows_acquired=True,
        )
        self.assertFalse(out.fit_allowed)
        self.assertEqual(out.status,"ABSTAIN_NO_SAME_DONOR_INTERSECTION")

    def test_missing_view_has_explicit_status(self):
        out=assess_same_donor_hepatic_multiview(
            source_id="TKI21",
            rna_donors=[],
            protein_donors=["D1","D2"],
            outcome_donors=["D1","D2"],
            outcome_kind="drug_specific_clint",
            source_rows_acquired=True,
        )
        self.assertEqual(out.status,"ABSTAIN_NO_MATCHED_RNA")

    def test_article_level_claims_do_not_unlock_fit(self):
        out=assess_same_donor_hepatic_multiview(
            source_id="JAMWAL74",
            rna_donors=None,
            protein_donors=None,
            outcome_donors=None,
            outcome_kind="drug_specific_metabolic_activity",
            source_rows_acquired=False,
        )
        self.assertEqual(out.status,"ABSTAIN_SOURCE_ROWS_NOT_ACQUIRED")
        self.assertFalse(out.fit_allowed)

    def test_duplicate_biological_ids_are_rejected(self):
        with self.assertRaisesRegex(HepaticMultiviewReadinessError,"Duplicated"):
            assess_same_donor_hepatic_multiview(
                source_id="X",
                rna_donors=["D1","D1"],
                protein_donors=["D1"],
                outcome_donors=["D1"],
                outcome_kind="drug_specific_clint",
                source_rows_acquired=True,
            )

    def test_known_candidate_summary_codes(self):
        self.assertEqual(
            qualify_candidate_summary(
                "EPA100",has_rna=True,has_protein=False,
                has_drug_outcome=True,individual_rows_acquired=True),
            "ABSTAIN_NO_MATCHED_PROTEIN")
        self.assertEqual(
            qualify_candidate_summary(
                "TKI21",has_rna=False,has_protein=True,
                has_drug_outcome=True,individual_rows_acquired=True),
            "ABSTAIN_NO_MATCHED_RNA")
        self.assertEqual(
            qualify_candidate_summary(
                "JAMWAL74",has_rna=True,has_protein=True,
                has_drug_outcome=True,individual_rows_acquired=False),
            "ABSTAIN_SOURCE_ROWS_NOT_ACQUIRED")


if __name__ == "__main__":
    unittest.main()
