"""Synthetic-only HLC original-source CV tests. No real donor observations."""
from __future__ import annotations

import unittest
import numpy as np

from omics2pk.hlc_real_cyp_validation import (
    evaluate_fold_matched_cyp,OriginalHLCQualificationError
)
from omics2pk.hlc_real_fullarray_validation import evaluate_original_array_PCA8


def fictional():
    rng=np.random.default_rng(204)
    n=65
    age=rng.integers(23,70,size=n)
    sex=rng.integers(0,2,size=n)
    clinical=np.column_stack((age,sex)).astype(float)
    target=rng.normal(0,1,size=(n,2))
    panel=np.column_stack((target,rng.normal(0,1,size=(n,16))))
    activity=np.exp(1+.6*target[:,0]+.02*age+.1*sex)
    whole=np.column_stack((panel,rng.normal(0,1,size=(n,25)))).astype("float32")
    return clinical,target,panel,activity,whole


class OriginalHumanLiverLeakageTests(unittest.TestCase):
    def test_fictional_source_RNA_target_nuisance_and_age_stratified_permutations(self):
        c,t,p,y,_=fictional()
        r=evaluate_fold_matched_cyp(c,t,p,y,seed=20261018,permutations=19)
        self.assertEqual(r["n_original_human_livers"],65)
        self.assertEqual(r["sex_and_age_quartile_stratified_RNA_shuffles"],19)
        self.assertIn("age_sex_plus_target_CYP_RNA",r["MAE_by_model_original_source_units"])
        self.assertLess(r["delta_MAE_target_RNA_minus_age_sex"],0)
        self.assertTrue(r["no_individual_original_predictions_or_source_values_exported"])

    def test_liver_source_must_have_original_measured_activity_on_every_eligible_donor(self):
        c,t,p,y,_=fictional()
        y[1]=np.nan
        with self.assertRaisesRegex(OriginalHLCQualificationError,"Missing"):
            evaluate_fold_matched_cyp(c,t,p,y,permutations=2)

    def test_original_source_same_person_row_count_required(self):
        c,t,p,y,_=fictional()
        with self.assertRaisesRegex(OriginalHLCQualificationError,"matched"):
            evaluate_fold_matched_cyp(c,t[:-1],p,y,permutations=2)

    def test_fullarray_pca_each_original_fold_training_only(self):
        c,t,p,y,w=fictional()
        r=evaluate_original_array_PCA8(c,w,np.array([0,1]),y)
        self.assertEqual(r["n_original_heldout_biological_liver_samples"],65)
        self.assertEqual(r["n_original_unannotated_microarray_probes"],43)
        self.assertTrue(all(np.isfinite(list(r["MAE_by_model_in_original_assay_units"].values()))))
        self.assertTrue(r["no_original_individual_samples_or_measured_values_redistributed"])

    def test_fullarray_refuses_non_quantitative_measured_CYP(self):
        c,t,p,y,w=fictional()
        y[2]=float("nan")
        with self.assertRaises(OriginalHLCQualificationError):
            evaluate_original_array_PCA8(c,w,np.array([0,1]),y)


if __name__=="__main__":
    unittest.main()
