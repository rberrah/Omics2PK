"""Synthetic-only checks: per-human posterior, exact IV physics and convergence."""
import math
import unittest
import importlib.util

PK_NUMERICS_AVAILABLE = bool(importlib.util.find_spec("numpy") and
                             importlib.util.find_spec("scipy"))
if PK_NUMERICS_AVAILABLE:
    import numpy as np
    from omics2pk.adaptive_iv_nlme import (
        MarginalLikelihoodError, adaptive_subject_loglik,
        adaptive_total_loglik, log_concentration_per_record,
    )
    from omics2pk.iv_pk_structural import iv_2c_ng_ml
    from scripts.benchmark_mavoglurant_iv_nlme import _data, _marginal_loglike
    from scripts.benchmark_mavoglurant_iv_adaptive_nlme import fit_adaptive


def profiles(n=5):
    rng = np.random.default_rng(2026)
    records = []
    for human in range(n):
        eta = rng.normal(0, 0.25)
        for occasion in (0, 1):
            dose = 12.0 + occasion * 3.0
            points = tuple(
                (time, iv_2c_ng_ml(time, dose, 0.2, 10*math.exp(eta), 28, 15, 70)
                 * math.exp(rng.normal(0, 0.18)))
                for time in (0.3, 0.6, 1.1, 3, 9, 24)
            )
            records.append({"subject": str(human), "occasion": str(occasion),
                            "dose_mg": dose, "duration_h": 0.2, "points": points})
    return records


@unittest.skipUnless(PK_NUMERICS_AVAILABLE, 'Optional pinned PK numerical extras not installed')
class AdaptiveMarginalTests(unittest.TestCase):
    def test_analytic_iv_2c_matches_zero_eta(self):
        data = _data(profiles())
        pars = (10., 28., 15., 70., 0.25, 0.18)
        logs = log_concentration_per_record(pars, data[0], data[2], np.zeros(len(data[1])))
        reference = np.array([
            math.log(iv_2c_ng_ml(t, d, duration, 10., 28., 15., 70.))
            for t, d, duration in zip(data[0][0], data[0][2], data[0][3])
        ])
        self.assertTrue(np.allclose(logs, reference, atol=1e-10))

    def test_adaptive_density_is_finite_and_stable_by_order(self):
        data = _data(profiles())
        pars = (10., 28., 15., 70., 0.32, 0.18)
        seven = adaptive_total_loglik(pars, data, 7)
        eleven = adaptive_total_loglik(pars, data, 11)
        twenty_one = adaptive_total_loglik(pars, data, 21)
        self.assertTrue(all(math.isfinite(x) for x in (seven, eleven, twenty_one)))
        self.assertLess(abs(eleven-twenty_one)/len(data[1]), 0.15)
        self.assertEqual(len(adaptive_subject_loglik(pars, data, 11)), 5)

    def test_both_occasions_share_same_random_clearance(self):
        data = _data(profiles(4))
        pars = (10., 28., 15., 70., 0.32, 0.18)
        original = adaptive_total_loglik(pars, data, 13)
        split = [{**p, "subject": p["subject"] + "-" + p["occasion"]}
                 for p in profiles(4)]
        fake = adaptive_total_loglik(pars, _data(split), 13)
        self.assertGreater(abs(original-fake), 0.05)

    def test_invalid_and_nonobserved_concentrations_fail_closed(self):
        rows = profiles(4)
        rows[0] = {**rows[0], "points": ((0.3, 0.0),) + rows[0]["points"][1:]}
        with self.assertRaises(ValueError):
            _data(rows)
        data = _data(profiles(4))
        with self.assertRaises(MarginalLikelihoodError):
            adaptive_total_loglik((10.,28.,15.,70.,-0.5,0.2), data)

    def test_training_fit_on_synthetic_people(self):
        records = profiles(12)
        params = fit_adaptive(records, pooled=[10.,28.,15.,70.], nodes=9,maxiter=120)
        self.assertTrue(np.all(np.isfinite(params)))
        self.assertEqual(len(params), 6)
        self.assertGreater(params[4],0)


if __name__ == "__main__":
    unittest.main()
