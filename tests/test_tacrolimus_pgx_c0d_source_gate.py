"""Synthetic tests for tacrolimus clinical+PGx C0/D source qualification."""
import copy
import unittest

from omics2pk.tacrolimus_pgx_c0d_source_gate import (
    TacrolimusPGxBaselineError,
    qualify_tacrolimus_clinical_pgx_c0d,
)


def rows():
    out=[]
    for i in range(12):
        out.append({
            "subject_id":f"SYN_{i:02d}",
            "age_years":25+i,
            "sex":"male" if i%2 else "female",
            "weight_kg":55+i,
            "hematocrit_fraction":0.30+0.005*i,
            "c0_ng_ml":9.5+0.2*i,
            "daily_dose_mg":6.0+0.25*(i%4),
            "stable_maintenance":True,
            "rs776746_effect_allele":"T",
            "rs776746_effect_allele_dosage":i%3,
            "rs75125371_effect_allele":"C",
            "rs75125371_effect_allele_dosage":0 if i<7 else (1 if i<11 else 2),
        })
    return out


class TacrolimusPGxC0DGateTests(unittest.TestCase):
    def test_complete_same_subject_source_gate(self):
        out=qualify_tacrolimus_clinical_pgx_c0d(
            rows(),require_rs75125371=True)
        self.assertEqual(out["n_unique_subjects"],12)
        self.assertEqual(out["n_female"],6)
        self.assertEqual(out["n_male"],6)
        self.assertEqual(out["rs776746_effect_allele"],"T")
        self.assertEqual(out["rs776746_dosage_counts"],{"0":4,"1":4,"2":4})
        self.assertTrue(out["rs75125371_included"])
        self.assertEqual(out["rs75125371_effect_allele"],"C")
        self.assertTrue(out["stable_maintenance_required"])
        self.assertEqual(out["hematocrit_unit"],"fraction")
        self.assertFalse(out["individual_subject_ids_or_values_exported"])
        self.assertTrue(out["published_R2_is_in_sample_explanatory_not_heldout_prediction"])
        self.assertNotIn("SYN_00",str(out))

    def test_duplicate_subject_fails_closed(self):
        x=rows()
        x[1]["subject_id"]=x[0]["subject_id"]
        with self.assertRaisesRegex(TacrolimusPGxBaselineError,"Duplicate"):
            qualify_tacrolimus_clinical_pgx_c0d(x)

    def test_hct_percent_or_permille_is_not_silently_guessed(self):
        x=rows()
        x[0]["hematocrit_fraction"]=31.0
        with self.assertRaisesRegex(TacrolimusPGxBaselineError,"explicit fraction"):
            qualify_tacrolimus_clinical_pgx_c0d(x)
        x=rows()
        x[0]["hematocrit_fraction"]=310
        with self.assertRaisesRegex(TacrolimusPGxBaselineError,"explicit fraction"):
            qualify_tacrolimus_clinical_pgx_c0d(x)

    def test_stable_maintenance_is_required_by_default(self):
        x=rows()
        x[0]["stable_maintenance"]=False
        with self.assertRaisesRegex(TacrolimusPGxBaselineError,"stable maintenance"):
            qualify_tacrolimus_clinical_pgx_c0d(x)
        out=qualify_tacrolimus_clinical_pgx_c0d(
            x,require_stable_maintenance=False)
        self.assertFalse(out["stable_maintenance_required"])

    def test_effect_allele_orientation_must_be_constant(self):
        x=rows()
        x[0]["rs776746_effect_allele"]="C"
        with self.assertRaisesRegex(TacrolimusPGxBaselineError,"orientation"):
            qualify_tacrolimus_clinical_pgx_c0d(x)
        x=rows()
        x[0]["rs75125371_effect_allele"]="T"
        with self.assertRaisesRegex(TacrolimusPGxBaselineError,"rs75125371"):
            qualify_tacrolimus_clinical_pgx_c0d(
                x,require_rs75125371=True)

    def test_genotype_must_vary_for_predictive_comparison(self):
        x=rows()
        for row in x:
            row["rs776746_effect_allele_dosage"]=0
        with self.assertRaisesRegex(TacrolimusPGxBaselineError,"no genotype variation"):
            qualify_tacrolimus_clinical_pgx_c0d(x)
        x=rows()
        for row in x:
            row["rs75125371_effect_allele_dosage"]=0
        with self.assertRaisesRegex(TacrolimusPGxBaselineError,"rs75125371"):
            qualify_tacrolimus_clinical_pgx_c0d(
                x,require_rs75125371=True)

    def test_nonpositive_target_inputs_rejected(self):
        for field in ("weight_kg","c0_ng_ml","daily_dose_mg"):
            x=rows()
            x[0][field]=0
            with self.assertRaises(TacrolimusPGxBaselineError):
                qualify_tacrolimus_clinical_pgx_c0d(x)

    def test_candidate_variant_can_be_omitted_without_fabrication(self):
        x=rows()
        for row in x:
            row.pop("rs75125371_effect_allele")
            row.pop("rs75125371_effect_allele_dosage")
        out=qualify_tacrolimus_clinical_pgx_c0d(x)
        self.assertFalse(out["rs75125371_included"])
        self.assertIsNone(out["rs75125371_effect_allele"])

    def test_partial_candidate_variant_columns_fail_closed(self):
        x=rows()
        x[0].pop("rs75125371_effect_allele_dosage")
        with self.assertRaises(TacrolimusPGxBaselineError):
            qualify_tacrolimus_clinical_pgx_c0d(x)


if __name__=="__main__":
    unittest.main()
