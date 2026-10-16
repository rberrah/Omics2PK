"""Synthetic-only original-source NAT/SMX anti-leakage tests."""
import copy
import unittest

from omics2pk.smx_nat_local_validation import (
    NATSourceError, _ridge, evaluate, qualify_original,
)


def fake_rows():
    rows=[]
    categories = ["slow"]*9 + ["intermediate"]*6 + ["rapid"]*5 + [""]*6
    for i in range(26):
        have_probe = i < 21
        rows.append({
            "donor_id":f"D{i:03d}",
            "row_kind":"individual",
            "source_doi":"10.1016/j.dmd.2025.100195",
            "na_smx_pmol_min_million_cells":str(0.10 + 0.002*i),
            "na_paba_pmol_min_million_cells":str(7+i) if have_probe else "",
            "na_smz_pmol_min_million_cells":str(5+i*0.4) if have_probe else "",
            "nat2_peptide_eflnshllpk_pmol_mg":str(0.2+i*.02),
            "nat2_peptide_dntdlvefk_pmol_mg":str(0.3+i*.015),
            "nat2_genotype_phenotype":categories[i],
            "age_years":str(18+i),
            "sex":"Female" if i%2 else "Male",
            "bmi":str(22+i*.2),
        })
    return rows


class NATDonorEvidenceTests(unittest.TestCase):
    def test_genuine_26_donors_and_separate_probe_genotype_coverage(self):
        result=evaluate(fake_rows())
        self.assertEqual(result["original_biological_donors"],26)
        self.assertEqual(result["arms"]["nat_probes"]["n_original_complete_case_donors"],21)
        self.assertEqual(result["arms"]["nat2_genotype"]["n_original_complete_case_donors"],20)
        self.assertEqual(result["arms"]["nat2_proteins"]["n_original_complete_case_donors"],26)
        shared=result["strict_four_arm_shared_donor_comparison"]
        self.assertEqual(shared["n_identical_complete_cases"],20)
        self.assertEqual(shared["status"],"IDENTICAL_DONOR_LODO_EXPLORATORY")
        self.assertEqual(set(shared["MAE"]),{"baseline","clinical","nat_probes","nat2_proteins","nat2_genotype"})
        for arm in result["arms"].values():
            self.assertTrue(arm["status"] == "DONOR_LODO_EXPLORATORY")
        self.assertTrue(result["original_rows_redistributed"] is False)

    def test_duplicate_donor_or_pool_is_not_a_27th_record(self):
        rows=fake_rows()
        rows[2]["donor_id"]=rows[0]["donor_id"]
        with self.assertRaisesRegex(NATSourceError,"Duplicate"):
            qualify_original(rows)
        rows=fake_rows()
        rows[-1]["row_kind"]="pool"
        with self.assertRaisesRegex(NATSourceError,"Pooled"):
            qualify_original(rows)

    def test_genotyped_subset_must_not_be_described_as_all_26(self):
        rows=fake_rows()
        rows[-1]["nat2_genotype_phenotype"]="rapid"
        with self.assertRaisesRegex(NATSourceError,"20 genotyped"):
            qualify_original(rows)

    def test_probe_assay_subset_must_not_expand_to_26(self):
        rows=fake_rows()
        rows[25]["na_smz_pmol_min_million_cells"]="30.2"
        rows[25]["na_paba_pmol_min_million_cells"]="15"
        with self.assertRaisesRegex(NATSourceError,"21 original donors"):
            qualify_original(rows)

    def test_missing_peptides_are_abstained_per_arm_not_imputed(self):
        rows=fake_rows()
        for r in rows[:18]:
            r["nat2_peptide_dntdlvefk_pmol_mg"]=""
        result=evaluate(rows)
        self.assertEqual(result["arms"]["nat2_proteins"]["status"],
                         "ABSTAIN_INSUFFICIENT_COMPLETE_CASES")
        self.assertEqual(result["arms"]["nat_probes"]["n_original_complete_case_donors"],21)
        shared=result["strict_four_arm_shared_donor_comparison"]
        self.assertEqual(shared["status"],"ABSTAIN_INSUFFICIENT_IDENTICAL_COMPLETE_CASE_DONORS")
        self.assertLess(shared["n_identical_complete_cases"],10)

    def test_zero_or_nonquantifiable_y_cannot_be_imputed(self):
        rows=fake_rows()
        rows[0]["na_smx_pmol_min_million_cells"]="0"
        with self.assertRaisesRegex(NATSourceError,"positive"):
            qualify_original(rows)

    def test_train_only_scaling_no_heldout_y_dependency(self):
        x=[[i,i*.2+1] for i in range(1,12)]
        y=[0.2+i*.04 for i in range(1,12)]
        p1=_ridge(x,y,[20,5.0])
        p2=_ridge(x,y,[20,5.0])
        self.assertAlmostEqual(p1,p2,places=12)
        self.assertGreaterEqual(p1,0)


if __name__=="__main__":
    unittest.main()
