"""Synthetic-only multi-analyte PBPK source profile mathematical invariants."""
import unittest
from scripts.validate_apap_three_analyte_literature_nca import compute_area


class OriginalThreeAnalyteTests(unittest.TestCase):
    def test_only_original_interpolation_between_observed_timestamps(self):
        result=compute_area([(0,1),(1,2),(3,1)])
        self.assertEqual(result["n_source_time_points"],3)
        self.assertAlmostEqual(result["observed_AUC_first_to_last_mg_h_l"],4.5)
        self.assertEqual(result["tmax_h"],1.0)
        self.assertEqual(result["cmax_mg_l"],2.0)

    def test_duplicates_and_sparse_traces_must_abstain(self):
        for trace in ([(0,0),(1,1)],[(0,1),(1,2),(1,3)],
                      [(0,1),(1,-2),(2,3)]):
            self.assertIsNone(compute_area(trace))

    def test_distinguish_real_first_sample_from_invented_t0(self):
        result=compute_area([(2,2),(3,4),(5,1)])
        self.assertEqual(result["first_time_h"],2)
        self.assertEqual(result["last_time_h"],5)
        self.assertAlmostEqual(result["observed_AUC_first_to_last_mg_h_l"],8.5)


if __name__=="__main__":
    unittest.main()
