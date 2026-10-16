"""Synthetic profile QA; never fit drug elimination from a fabricated AUCinf."""
import unittest
from omics2pk.single_dose_nca import single_dose_descriptive_nca, LongitudinalNCAError

class SingleDoseNCATests(unittest.TestCase):
    def test_linear_trapezoid_and_nonzero_baseline_are_preserved(self):
        result=single_dose_descriptive_nca([(4.,2.),(0.,0.5),(2.,2.)])
        self.assertAlmostEqual(result["AUC0_last_original_mg_h_l_linear_trapezoid"],6.5)
        self.assertEqual(result["Cmax_original_mg_l"],2.)
        self.assertEqual(result["Tmax_h"],2.)
        self.assertEqual(result["last_original_sample_h"],4.)
        self.assertTrue(result["does_not_impute_terminal_AUC"])

    def test_zero_or_negative_OR_BLQ_imputed_unacceptable(self):
        for observations in (
            [(0,0),(2,1),(2,3)],
            [(0,0),(1,-0.001),(2,1)],
            [(0,0),(1,float("nan")),(2,1)],
            [(1,0),(2,1),(3,1)],
        ):
            with self.assertRaises(LongitudinalNCAError):
                single_dose_descriptive_nca(observations)

    def test_too_sparse_cannot_imply_pk_curve(self):
        with self.assertRaisesRegex(LongitudinalNCAError,"At least three"):
            single_dose_descriptive_nca([(0,1),(2,2)])

if __name__=="__main__":
    unittest.main()
