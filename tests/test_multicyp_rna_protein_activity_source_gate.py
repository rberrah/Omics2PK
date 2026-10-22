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
                "protein_value":str(protein),
                "protein_unit":"pmol/mg",
                "activity_value":str(activity),
                "activity_unit":"pmol/min/mg",
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
        self.assertEqual(out["expected_complete_donors_per_included_cyp"],17)
        self.assertEqual(out["results_by_cyp"]["CYP3A4"]["n_held_donor_folds"],17)
        with self.assertRaisesRegex(MultiCYPSourceNotQualified,"exactly 17"):
            evaluate_original_source(rows[:-1])

    def test_zhang_partial_donor_subset_is_rejected(self):
        rows=rows_for("10.1124/jpet.116.233635",99)
        with self.assertRaisesRegex(MultiCYPSourceNotQualified,"exactly 100"):
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
