"""Strict source-only tests; source values are intentionally synthetic."""
import unittest
from scripts.validate_public_tacrolimus_36_c0_baseline import (
    EXPECTED,TacrolimusClinicalSourceError,
    qualify_clinical_rows,baseline_loocv
)


def fake():
    rows=[list(EXPECTED)]
    for i in range(1,37):
        dose=3+0.5*(i%5)
        weight=65+2*(i%10)
        dose_kg=dose/weight
        concentration=5+0.2*i
        rows.append([
            i,45+i%25, "Male" if i%2 else "Female","SYNTH_STAGE",
            dose,170,weight,weight/2.89,
            dose_kg,concentration,concentration/dose_kg,
            "NONE","NONE","SYNTH","SYNTH","SYNTH"
        ])
    return rows


class SyntheticHumanTac2026Tests(unittest.TestCase):
    def test_36_original_units_source_formula_and_patient_heldout_fit(self):
        rows=qualify_clinical_rows(fake())
        self.assertEqual(len(rows),36)
        result=baseline_loocv(rows)
        self.assertEqual(result["n_patient_disjoint_holdout_folds"],36)
        self.assertTrue(result["no_omics_model_fitted"])
        self.assertNotIn("SYNTH_STAGE",repr(result))

    def test_mismatched_dose_per_weight_rejected(self):
        r=fake();r[10][8]=0.001
        with self.assertRaisesRegex(TacrolimusClinicalSourceError,"dose per kg"):
            qualify_clinical_rows(r)

    def test_source_clinical_covariate_not_normalized_by_future(self):
        r=fake();r[10][10]*=3
        with self.assertRaisesRegex(TacrolimusClinicalSourceError,"C/D"):
            qualify_clinical_rows(r)

    def test_count_or_schema_drift_rejected(self):
        with self.assertRaisesRegex(TacrolimusClinicalSourceError,"layout"):
            qualify_clinical_rows(fake()[:30])
        r=fake();r[0][0]="Not the original DOI source header"
        with self.assertRaisesRegex(TacrolimusClinicalSourceError,"layout"):
            qualify_clinical_rows(r)

    def test_source_patient_duplicate_rejected(self):
        r=fake();r[9][0]=r[3][0]
        with self.assertRaisesRegex(TacrolimusClinicalSourceError,"distinct"):
            qualify_clinical_rows(r)

if __name__=="__main__":unittest.main()
