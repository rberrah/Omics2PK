"""Synthetic-only tests for the authenticated multi-CYP source gate."""
import unittest

from omics2pk.multicyp_rna_protein_activity_source_gate import (
    LABELS,
    MultiCYPSourceNotQualified,
    evaluate_original_source,
    source_contract,
)


def rows_for(doi, n, cyps=("CYP3A4",), protein_driven=True):
    rows=[]
    for cidx,cyp in enumerate(cyps):
        for i in range(n):
            # Deterministic positive values; no original donor data.
            rna=1.0 + ((i*7+cidx*3)%23)/7
            protein=2.0 + ((i*11+cidx*5)%29)/5
            activity=(4.0*protein + 0.2*rna + (i%5)*0.03
                      if protein_driven else 4.0*rna + 0.2*protein + (i%5)*0.03)
            rows.append({
                "donor_id":f"SYNTH_{i:03d}",
                "source_doi":doi,
                "row_kind":"individual",
                "cyp_isoform":cyp,
                "probe_substrate":"midazolam" if cyp=="CYP3A4" else "synthetic_probe",
                "rna_value":str(rna),
                "rna_unit":"relative_expression",
                "rna_status":"observed",
                "protein_value":str(protein),
                "protein_unit":"pmol/mg",
                "protein_status":"observed",
                "activity_value":str(activity),
                "activity_unit":"pmol/min/mg",
                "activity_status":"observed",
            })
    return rows


class MultiCYPSourceGateTests(unittest.TestCase):
    def test_zhang_n100_two_cyp_source_contract_and_benchmark(self):
        rows=rows_for("10.1124/jpet.116.233635",100,("CYP3A4","CYP2C9"))
        doi,profile,by_cyp,units=source_contract(rows)
        self.assertEqual(doi,"10.1124/jpet.116.233635")
        self.assertEqual(profile["expected_donors_per_included_cyp"],100)
        self.assertEqual(set(by_cyp),{"CYP3A4","CYP2C9"})
        self.assertEqual(len(by_cyp["CYP3A4"]),100)
        out=evaluate_original_source(rows)
        self.assertEqual(out["n_included_cyp_endpoints"],2)
        self.assertEqual(out["results_by_cyp"]["CYP3A4"]["status"],
                         "QUALIFIED_COMPLETE_CASE_LODO")
        self.assertTrue(out["not_clinical_plasma_PK"])
        self.assertTrue(out["no_cross_study_fusion"])
        self.assertFalse(out["original_individual_rows_or_predictions_exported"])
        self.assertNotIn("SYNTH_",str(out))
        # With this synthetic construction protein contains the stronger signal.
        self.assertLess(
            out["results_by_cyp"]["CYP3A4"]["mae_by_model"]["protein_only"],
            out["results_by_cyp"]["CYP3A4"]["mae_by_model"]["RNA_only"])

    def test_ohtsuki_requires_exact_17_donors_for_included_cyp(self):
        rows=rows_for("10.1124/dmd.111.042259",17)
        out=evaluate_original_source(rows)
        self.assertEqual(out["expected_source_donors_per_included_cyp"],17)
        self.assertEqual(out["minimum_complete_donors_required_per_included_cyp"],8)
        self.assertEqual(out["results_by_cyp"]["CYP3A4"]["n_held_donor_folds"],17)
        self.assertEqual(out["results_by_cyp"]["CYP3A4"]["n_source_donors"],17)
        with self.assertRaisesRegex(MultiCYPSourceNotQualified,"exactly 17"):
            evaluate_original_source(rows[:-1])


    def test_explicit_blq_keeps_all_source_donors_but_uses_same_complete_cases(self):
        rows=rows_for("10.1124/dmd.111.042259",17)
        rows[0]["protein_value"]=""
        rows[0]["protein_status"]="blq"
        rows[1]["rna_value"]=""
        rows[1]["rna_status"]="missing"
        out=evaluate_original_source(rows)
        cyp=out["results_by_cyp"]["CYP3A4"]
        self.assertEqual(cyp["status"],"QUALIFIED_COMPLETE_CASE_LODO")
        self.assertEqual(cyp["n_source_donors"],17)
        self.assertEqual(cyp["n_complete_case_donors"],15)
        self.assertEqual(cyp["n_held_donor_folds"],15)
        self.assertEqual(cyp["missingness"]["protein_nonobserved"],1)
        self.assertEqual(cyp["missingness"]["rna_nonobserved"],1)
        self.assertEqual(cyp["missingness"]["any_nonobserved"],2)
        self.assertEqual(cyp["complete_case_policy"],
                         "no imputation; same donors in every compared arm")

    def test_nonobserved_status_cannot_hide_numeric_value_or_coerce_zero(self):
        rows=rows_for("10.1124/dmd.111.042259",17)
        rows[0]["protein_status"]="blq"
        with self.assertRaisesRegex(MultiCYPSourceNotQualified,"must be blank"):
            source_contract(rows)
        rows=rows_for("10.1124/dmd.111.042259",17)
        rows[0]["protein_value"]="0"
        with self.assertRaisesRegex(MultiCYPSourceNotQualified,"strictly positive"):
            source_contract(rows)

    def test_insufficient_complete_donors_abstains_without_imputation(self):
        rows=rows_for("10.1124/dmd.111.042259",17)
        for row in rows[:10]:
            row["protein_value"]=""
            row["protein_status"]="blq"
        out=evaluate_original_source(rows)
        cyp=out["results_by_cyp"]["CYP3A4"]
        self.assertEqual(cyp["status"],"ABSTAIN_INSUFFICIENT_COMPLETE_DONORS")
        self.assertEqual(cyp["n_source_donors"],17)
        self.assertEqual(cyp["n_complete_case_donors"],7)
        self.assertNotIn("mae_by_model",cyp)

    def test_unknown_or_blank_measurement_status_is_rejected(self):
        rows=rows_for("10.1124/dmd.111.042259",17)
        rows[0]["rna_status"]="estimated"
        with self.assertRaisesRegex(MultiCYPSourceNotQualified,"Invalid measurement status"):
            source_contract(rows)
        rows=rows_for("10.1124/dmd.111.042259",17)
        rows[0]["activity_status"]=""
        with self.assertRaisesRegex(MultiCYPSourceNotQualified,"Invalid measurement status"):
            source_contract(rows)

    def test_zhang_partial_donor_subset_is_rejected(self):
        rows=rows_for("10.1124/jpet.116.233635",99)
        with self.assertRaisesRegex(MultiCYPSourceNotQualified,"exactly 100"):
            source_contract(rows)

    def test_included_cyps_must_retain_identical_source_donor_cohort(self):
        rows=rows_for("10.1124/dmd.111.042259",17,("CYP3A4","CYP2C9"))
        for row in rows:
            if row["cyp_isoform"]=="CYP2C9" and row["donor_id"]=="SYNTH_016":
                row["donor_id"]="SYNTH_REPLACEMENT"
                break
        with self.assertRaisesRegex(
                MultiCYPSourceNotQualified,"same source donor cohort"):
            source_contract(rows)

    def test_duplicate_donor_cyp_is_rejected(self):
        rows=rows_for("10.1124/dmd.111.042259",17)
        rows.append(dict(rows[0]))
        with self.assertRaisesRegex(MultiCYPSourceNotQualified,"Duplicate"):
            source_contract(rows)

    def test_mixed_units_within_cyp_are_rejected(self):
        rows=rows_for("10.1124/dmd.111.042259",17)
        rows[4]["protein_unit"]="fmol/ug"
        with self.assertRaisesRegex(MultiCYPSourceNotQualified,"Units vary"):
            source_contract(rows)

    def test_multiple_probe_assays_for_one_cyp_are_not_mixed(self):
        rows=rows_for("10.1124/dmd.111.042259",17)
        rows[3]["probe_substrate"]="nifedipine"
        with self.assertRaisesRegex(MultiCYPSourceNotQualified,"one prespecified probe"):
            source_contract(rows)

    def test_unknown_source_and_summary_rows_fail_closed(self):
        rows=rows_for("10.1124/dmd.111.042259",17)
        rows[0]["source_doi"]="10.0000/not-preregistered"
        with self.assertRaises(MultiCYPSourceNotQualified):
            source_contract(rows)
        rows=rows_for("10.1124/dmd.111.042259",17)
        rows[0]["row_kind"]="mean"
        with self.assertRaisesRegex(MultiCYPSourceNotQualified,"Pooled"):
            source_contract(rows)

    def test_schema_and_positive_measurements_are_strict(self):
        rows=rows_for("10.1124/dmd.111.042259",17)
        rows[0].pop("activity_unit")
        with self.assertRaisesRegex(MultiCYPSourceNotQualified,"exact"):
            source_contract(rows)
        rows=rows_for("10.1124/dmd.111.042259",17)
        rows[0]["rna_value"]="0"
        with self.assertRaisesRegex(MultiCYPSourceNotQualified,"strictly positive"):
            source_contract(rows)


if __name__=="__main__":
    unittest.main()
