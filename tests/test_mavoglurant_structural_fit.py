"""Synthetic-only numerical fitting smoke; original data is accessed in CI runtime."""
import importlib.util
import math
import unittest
from omics2pk.iv_pk_structural import iv_1c_ng_ml
from scripts.benchmark_mavoglurant_iv_structural import heldout_comparison


@unittest.skipUnless(importlib.util.find_spec("scipy") and
                     importlib.util.find_spec("numpy"), "pk-structural optional dependencies absent")
class SyntheticHeldoutTests(unittest.TestCase):
    def test_pooled_fixed_effect_iv_validation_is_patient_disjoint(self):
        profiles=[]
        for patient in range(10):
            for occ in range(2):
                points=tuple(
                    (t,iv_1c_ng_ml(t,10+occ,1/6,8,30))
                    for t in (0.2,0.5,1,3,8,16,24)
                )
                profiles.append({
                    "subject":str(patient),"occasion":str(occ),
                    "dose_mg":10+occ,"duration_h":1/6,"points":points,
                })
        outcome=heldout_comparison(profiles,folds=5)
        self.assertTrue(outcome["no_nlme_or_popPK_claim"])
        self.assertEqual(outcome["original_biological_humans"],10)
        self.assertEqual(outcome["original_qualified_occasions"],20)
        self.assertEqual(outcome["original_measured_concentrations"],140)
        self.assertEqual(outcome["1C"]["n_biological_test_subjects_all_folds"],10)
        self.assertLess(outcome["1C"]["root_mean_individual_mean_squared_log_error"],0.001)


if __name__=="__main__": unittest.main()
