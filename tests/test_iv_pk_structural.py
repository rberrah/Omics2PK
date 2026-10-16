"""Synthetic-only mathematical and anti-leakage IV benchmark tests."""
import math
import unittest

from omics2pk.iv_pk_structural import (
    StructuralPKError, iv_1c_ng_ml, iv_2c_ng_ml, subject_disjoint_folds
)


class IVAnalyticalTests(unittest.TestCase):
    def test_1c_infusion_matches_closed_form_and_units(self):
        dose, duration, cl, volume = 10.0, 1/6, 5.0, 25.0
        result = iv_1c_ng_ml(duration, dose, duration, cl, volume)
        ref = 1000 * (dose / duration) / cl * (1 - math.exp(-cl/volume*duration))
        self.assertAlmostEqual(result, ref, places=8)
        self.assertEqual(iv_1c_ng_ml(0, dose, duration, cl, volume), 0)
        post = iv_1c_ng_ml(duration+3, dose, duration, cl, volume)
        self.assertAlmostEqual(post/result, math.exp(-cl/volume*3), places=10)

    def test_2c_approaches_1c_when_q_is_small(self):
        t,dose,T,cl,volume = 1.5,10,1/6,8,30
        c1 = iv_1c_ng_ml(t,dose,T,cl,volume)
        c2 = iv_2c_ng_ml(t,dose,T,cl,volume,1e-5,100)
        self.assertAlmostEqual(c2,c1,delta=1e-3)
        self.assertEqual(iv_2c_ng_ml(0,dose,T,cl,volume,5,100),0)

    def test_concentration_is_nonnegative_and_decays(self):
        for model in (lambda t: iv_1c_ng_ml(t,10,1/6,8,30),
                      lambda t: iv_2c_ng_ml(t,10,1/6,8,30,12,100)):
            values = [model(t) for t in (0,0.1,1/6,1,4,24)]
            self.assertTrue(all(math.isfinite(x) and x>=0 for x in values))
            self.assertGreater(values[2],values[3])
            self.assertGreater(values[3],values[4])
            self.assertGreater(values[4],values[5])

    def test_rejects_invalid_doses_parameters_and_times(self):
        with self.assertRaises(StructuralPKError):
            iv_1c_ng_ml(-1,10,1/6,8,30)
        with self.assertRaises(StructuralPKError):
            iv_2c_ng_ml(1,10,1/6,8,30,-2,100)
        with self.assertRaises(StructuralPKError):
            iv_1c_ng_ml(1,10,0,8,30)
        with self.assertRaises(StructuralPKError):
            iv_1c_ng_ml(float("nan"),10,1/6,8,30)

    def test_individual_folds_preserve_all_occasions_without_leakage(self):
        subjects=[str(i) for i in range(120)]
        folds = subject_disjoint_folds(subjects)
        self.assertEqual(folds,subject_disjoint_folds(reversed(subjects)))
        self.assertEqual(sorted(len(test) for _,test in folds),[24]*5)
        seen=[]
        occasions = {(s, o) for s in subjects for o in (1,2)}
        for train,test in folds:
            self.assertFalse(train & test)
            self.assertEqual({s for s,o in occasions if s in test},test)
            seen.extend(test)
        self.assertEqual(set(seen),set(subjects))
        self.assertEqual(len(seen),len(subjects))

    def test_invalid_subject_groups_are_not_silently_used(self):
        for subjects in (["1","1","2"], ["1",""], ["1"]):
            with self.assertRaises(StructuralPKError):
                subject_disjoint_folds(subjects,2)


if __name__=="__main__": unittest.main()
