"""Synthetic measured-PHH source-table and leakage-free 6-donor CV tests."""
import math
import unittest

from omics2pk.phh_cyp_evidence import CYPS, EvidenceError, extract_mean, held_out_biomarker_test, parse_published_functional_table

D=("330","409","416","211","025","789")
def make_table(kind):
    head=["Enzyme", "Substrate (Concentration µM)" if kind=="substrate_clearance" else "Metabolite"]+[f"PHH{d}" for d in D]
    rows=[head]
    for i,cyp in enumerate(CYPS):
        rows.append([cyp,"standard"]+[f"{i+j+1}.1 ± 0.2" for j in range(6)])
    return rows

class FunctionalTests(unittest.TestCase):
    def test_mean_uncensored(self):
        self.assertEqual(extract_mean("12.3 ± 0.4"),12.3)
        self.assertEqual(extract_mean("0.0004"),0.0004)
        for x in ("ND","<0.1"," > 3","n.a.","","-5","mean 2.3"):
            self.assertIsNone(extract_mean(x))

    def test_six_donors_two_assay_types_with_original_order(self):
        for kind in ("substrate_clearance","metabolite_formation"):
            table=parse_published_functional_table(make_table(kind),kind=kind)
            self.assertEqual(len(table),6)
            self.assertEqual(set(table["CYP2B6"]),set(D))
            self.assertEqual(table["CYP3A4"]["330"],6.1)

    def test_wrong_donor_or_duplicate(self):
        table=make_table("substrate_clearance")
        table[0][2]="PHH999"
        with self.assertRaisesRegex(EvidenceError,"Unrecognized"):
            parse_published_functional_table(table,kind="substrate_clearance")
        table=make_table("substrate_clearance")
        table[0][3]="PHH330"
        with self.assertRaisesRegex(EvidenceError,"[Ss]ix donor"):
            parse_published_functional_table(table,kind="substrate_clearance")

    def test_real_donor_splits_and_train_only_scaling(self):
        x={d:float(i) for i,d in enumerate(D)}
        y={d:float(i*2+1) for i,d in enumerate(D)}
        report=held_out_biomarker_test(x,y)
        self.assertEqual(report["n_independent_donors"],6)
        self.assertEqual(report["n_donor_heldout_folds"],6)
        self.assertGreater(report["baseline_mae_original_assay_units"],report["rna_mae_original_assay_units"])

    def test_feature_is_not_dose_response_without_omics(self):
        y={d:float(i) for i,d in enumerate(D)}
        with self.assertRaisesRegex(EvidenceError,"Require 4"):
            held_out_biomarker_test({"330":1.0,"211":2.0},y)
        report=held_out_biomarker_test({d:0.0 for d in D},y)
        self.assertEqual(report["rna_mae_original_assay_units"],report["baseline_mae_original_assay_units"])
        self.assertIsNone(report["donor_expression_activity_pearson_r_exploratory"])

    def test_nan_and_posthoc_penalty_rejected(self):
        y={d:float(i) for i,d in enumerate(D)}
        with self.assertRaisesRegex(EvidenceError,"observed nonnegative"):
            held_out_biomarker_test({d:float("nan") for d in D},y)
        with self.assertRaisesRegex(EvidenceError,"post-hoc"):
            held_out_biomarker_test({d:float(i) for i,d in enumerate(D)},y,fixed_ridge_alpha=0)

if __name__=="__main__":
    unittest.main()
