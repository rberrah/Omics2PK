"""Purely synthetic Topletz-style test fixture; NO actual original donor data."""
import unittest

from omics2pk.topletz_vinorelbine_source_gate import (
    DOI, TopletzSourceError, benchmark_topletz_table2,
)


def fixture():
    # Synthetic labels deliberately different from original Topletz HLM donors.
    dips = ["*1/*1"] + ["*1/*3"] * 6 + ["*1/*6"] + [
        "*3/*3"] * 5 + ["*1/*3"] + ["*3/*3"] * 2 + [
        "*1/*7"] + ["*3/*3"] * 3
    assert len(dips) == 20
    rows = []
    for i, dip in enumerate(dips):
        p5 = float(20 + i) if i < 8 else (0.0 if i == 9 else 0.5 + 0.1 * i)
        p4 = 20 + 9 * (i % 9) + (i / 10)
        rows.append({
            "source_doi":DOI,
            "donor_id":f"IUL-{i+101}",
            "cyp3a5_diplotype":dip,
            "original_measured_CYP3A5_expression_group":"high" if i < 8 else "low",
            "cyp3a4_protein_pmol_mg":p4,
            "cyp3a5_protein_pmol_mg":p5,
            "vinorelbine_rate_pmol_min_mg":1.0 + 0.28 * p4 + i * 0.03,
            # These two derived source values are intentionally included
            # to ensure they are NOT used to model the measured substrate rate.
            "derived_CLint":1000 + i,
            "derived_CLh":200 + i,
        })
    return rows


class TopletzSourceTests(unittest.TestCase):
    def test_all_arms_same_donors_and_no_derived_clearance_target(self):
        out=benchmark_topletz_table2(fixture())
        self.assertEqual(out["n_genuine_biological_donors"],20)
        self.assertEqual(out["source_protein_group_counts"],{"high":8,"low":12})
        self.assertEqual(len(out["heldout_MAE_by_arm"]),7)
        self.assertTrue(out["not_derived_Clint_or_CLh"])
        self.assertTrue(out["no_individual_rows_donor_ids_or_predictions_exported"])
        self.assertNotIn("IUL-101",str(out))

    def test_source_input_permutation_invariant(self):
        rows=fixture()
        a=benchmark_topletz_table2(rows)["heldout_MAE_by_arm"]
        b=benchmark_topletz_table2(list(reversed(rows)))["heldout_MAE_by_arm"]
        self.assertEqual(a,b)

    def test_target_derived_clearances_never_leak(self):
        rows=fixture()
        before=benchmark_topletz_table2(rows)["heldout_MAE_by_arm"]
        for row in rows:
            row["derived_CLint"] *= 999
            row["derived_CLh"] *= 999
        after=benchmark_topletz_table2(rows)["heldout_MAE_by_arm"]
        self.assertEqual(before,after)

    def test_invalid_or_duplicate_source_denied(self):
        rows=fixture()
        rows[0]["donor_id"]=rows[1]["donor_id"]
        with self.assertRaisesRegex(TopletzSourceError,"Duplicate"):
            benchmark_topletz_table2(rows)
        rows=fixture()[:-1]
        with self.assertRaisesRegex(TopletzSourceError,"exactly 20"):
            benchmark_topletz_table2(rows)

    def test_measured_zero_cyp3a5_protein_is_allowed_but_unobserved_missing_not(self):
        rows=fixture()
        self.assertEqual(rows[9]["cyp3a5_protein_pmol_mg"],0.0)
        benchmark_topletz_table2(rows)
        rows[9]["cyp3a5_protein_pmol_mg"]=None
        with self.assertRaisesRegex(TopletzSourceError,"measured CYP3A5 protein"):
            benchmark_topletz_table2(rows)

    def test_missing_measured_low_CYP3A5_in_star1_carriers_denied(self):
        rows=fixture()
        rows[13]["cyp3a5_diplotype"]="*3/*3"
        with self.assertRaisesRegex(TopletzSourceError,"low measured CYP3A5"):
            benchmark_topletz_table2(rows)

    def test_genotype_allele_coverage_or_doi_mismatch(self):
        rows=fixture()
        rows[0]["cyp3a5_diplotype"]="*1/*8"
        with self.assertRaisesRegex(TopletzSourceError,"Unresolved CYP3A5"):
            benchmark_topletz_table2(rows)
        rows=fixture()
        rows[0]["source_doi"]="10.1000/fiction"
        with self.assertRaisesRegex(TopletzSourceError,"DOI"):
            benchmark_topletz_table2(rows)

    def test_no_fixed_activity_percentages_implicit_in_regression(self):
        out=benchmark_topletz_table2(fixture())
        self.assertNotIn("activity_loss_pct",str(out))
        self.assertNotIn("genotype_percent_activity",str(out))
        self.assertFalse(out["no_automatic_source_byte_verification"] is False)


if __name__=="__main__":
    unittest.main()
