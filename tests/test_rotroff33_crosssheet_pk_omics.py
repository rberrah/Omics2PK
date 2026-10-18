"""Synthetic-only public Rotroff supplement cross-sheet patient-match screening."""
import unittest

from scripts.audit_rotroff33_crosssheet_pk_omics import id_candidates, summarize_sheet, PK_TERMS


class CrossSheetSyntheticTests(unittest.TestCase):
    def test_source_qualifying_wide_omics_and_clinical_candidates(self):
        omics=[{0:"Sample ID",1:"timepoint"}]
        for i in range(33):
            for t in ("A","B","C"):
                omics.append({0:f"FAKE_ID_{i}",1:t,
                              **{j+2:str((i+j+1)/11) for j in range(105)}})
        pk=[{0:"Subject ID",1:"Metformin AUC",2:"Cmax"}]
        for i in range(33):
            pk.append({0:f"FAKE_ID_{i}",1:str(100+i),2:str(10+i)})
        a=summarize_sheet(omics)
        b=summarize_sheet(pk)
        self.assertTrue(a["possible_omics_33x3_or_transposed_matrix"])
        self.assertFalse(a["possible_33person_rows_and_10plus_quantified_covariates"])
        self.assertTrue(b["header_pk_role_tokens_found_first_10_rows"]["AUC_token"])
        self.assertTrue(b["header_pk_role_tokens_found_first_10_rows"]["Cmax_token"])
        self.assertEqual(len(id_candidates(omics)&id_candidates(pk)),33)
        self.assertNotIn("FAKE_ID",str(a))
        self.assertNotIn("FAKE_ID",str(b))

    def test_demographic_metadata_without_pk_not_acceptable(self):
        pk=[{0:"subject",1:"glucose",2:"body mass"}]
        for i in range(33):
            pk.append({0:f"SYN_{i}",1:str(11+i),2:str(60+i)})
        report=summarize_sheet(pk)
        self.assertFalse(report["header_pk_role_tokens_found_first_10_rows"]["AUC_token"])
        self.assertFalse(report["header_pk_role_tokens_found_first_10_rows"]["Cmax_token"])

    def test_do_not_treat_numeric_values_as_subject_ids(self):
        rows=[{0:str(i),1:"male",2:"female"} for i in range(40)]
        self.assertEqual(len(id_candidates(rows)),0)

    def test_unjoined_ids_have_zero_intersection(self):
        a=[{0:f"SYN_A_{i}"} for i in range(33)]
        b=[{0:f"SYN_B_{i}"} for i in range(33)]
        self.assertEqual(len(id_candidates(a)&id_candidates(b)),0)


if __name__=="__main__":unittest.main()
