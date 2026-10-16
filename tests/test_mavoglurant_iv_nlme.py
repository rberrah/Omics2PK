"""Synthetic-only, subject-level NLME probability and anti-leakage tests."""
import importlib.util
import math
import unittest

from omics2pk.iv_pk_structural import iv_2c_ng_ml
from scripts.benchmark_mavoglurant_iv_nlme import (
    _data, _grid_prediction, _marginal_loglike, fit_iiv_cl,
    score_population_heldout,
)


@unittest.skipUnless(importlib.util.find_spec("numpy") and
                     importlib.util.find_spec("scipy"), "Pinned PK numerical extras absent")
class NLMEProofTests(unittest.TestCase):
    @staticmethod
    def example(n_people=12):
        import numpy as np
        rng = np.random.default_rng(25)
        profiles = []
        for i in range(n_people):
            eta = rng.normal(0, 0.32)
            for occ in (0, 1):
                dose = 10 + 4 * occ
                points = tuple(
                    (t, float(iv_2c_ng_ml(t, dose, 0.2, 9 * math.exp(eta), 25, 12, 75)
                              * math.exp(rng.normal(0, 0.09))))
                    for t in (0.3, 0.6, 1, 3, 8, 16, 24)
                )
                profiles.append({"subject": str(i), "occasion": str(occ),
                                 "dose_mg": dose, "duration_h": 0.2, "points": points})
        return profiles

    def test_joint_subject_not_independent_occasions(self):
        import numpy as np
        profiles = self.example(4)
        data = _data(profiles)
        pars = (9, 25, 12, 75, 0.32, 0.09)
        grid, weights = _grid_prediction(pars, data[0], 9)
        self.assertEqual(grid.shape, (9, 56))
        self.assertAlmostEqual(float(np.exp(weights).sum()), 1, places=12)
        self.assertEqual(len(data[1]), 4)
        shared = _marginal_loglike(pars, data)
        # Deliberately wrong independent-ETA-per-occasion comparator.
        wrong = [{**p, "subject": p["subject"] + "-" + p["occasion"]}
                 for p in profiles]
        separate = _marginal_loglike(pars, _data(wrong))
        self.assertGreater(abs(shared - separate), 0.01)

    def test_heldout_prediction_is_not_empirical_bayes(self):
        profiles = self.example(5)
        pars = (9, 25, 12, 75, 0.3, 0.15)
        original = score_population_heldout(pars, profiles)
        altered = [{**p, "points": tuple((t, dv * 2) for t, dv in p["points"])}
                   for p in profiles]
        changed = score_population_heldout(pars, altered)
        self.assertEqual(original["n_independent_heldout_humans"], 5)
        self.assertEqual(original["n_heldout_occasions"], 10)
        self.assertNotAlmostEqual(original["person_mse_log"], changed["person_mse_log"], places=4)
        self.assertAlmostEqual(changed["person_mean_bias_log"] - original["person_mean_bias_log"],
                               -math.log(2), places=10)

    def test_mdv_and_blq_cannot_enter_log_likelihood_as_zero(self):
        profiles = self.example(4)
        first = profiles[0]
        profiles[0] = {**first, "points": ((0.3, 0.0),) + first["points"][1:]}
        with self.assertRaisesRegex(ValueError, "positive measured"):
            _data(profiles)

    def test_fitted_iiv_is_positive_and_finite(self):
        profiles = self.example(12)
        pars = fit_iiv_cl(profiles, [9, 25, 12, 75], gh_nodes=7, maxiter=150)
        self.assertEqual(len(pars), 6)
        self.assertTrue(all(math.isfinite(x) and x > 0 for x in pars))
        self.assertLess(pars[4], 1.5)


if __name__ == "__main__":
    unittest.main()
