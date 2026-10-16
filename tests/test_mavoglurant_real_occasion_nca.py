"""Synthetic-only original mavoglurant IV single-occasion analytical invariants."""
import unittest
from scripts.validate_mavoglurant_120_real_occasion_nca import (
    OriginalPKCurveError,observed_curve
)

def fake():
    dose={"EVID":1,"MDV":1,"TIME":0,"AMT":20,"RATE":120,"DV":None}
    obs=[
        {"EVID":0,"MDV":0,"TIME":0.5,"DV":10},
        {"EVID":0,"MDV":0,"TIME":1.0,"DV":8},
        {"EVID":0,"MDV":0,"TIME":2.0,"DV":4},
    ]
    return [dose,*obs]

class ObservedMavoglurantTests(unittest.TestCase):
    def test_observed_profile_nca_is_bounded_by_real_samples(self):
        q=observed_curve(fake())
        self.assertEqual(q["n_original_samples"],3)
        self.assertAlmostEqual(q["source_infusion_duration_h"],1/6)
        self.assertAlmostEqual(q["AUC_first_to_last_original_DV_h"],10.5)
        self.assertEqual(q["first_original_observation_time_h"],0.5)
        self.assertEqual(q["Tmax_observed_h_from_dose"],0.5)
        self.assertFalse(q["AUC_0_infinity_calculated"])

    def test_multiple_doses_in_same_occasion_abstain(self):
        events=fake()
        events.append(dict(events[0],TIME=12))
        with self.assertRaisesRegex(OriginalPKCurveError,"exactly one"):
            observed_curve(events)

    def test_zero_rate_or_nonphysiologic_duration_abstain(self):
        events=fake();events[0]["RATE"]=0
        with self.assertRaisesRegex(OriginalPKCurveError,"positive"):
            observed_curve(events)

    def test_duplicate_times_and_missing_assay_abstain(self):
        events=fake()
        events[3]["TIME"]=1.0
        with self.assertRaisesRegex(OriginalPKCurveError,"Replicate"):
            observed_curve(events)
        events=fake()
        events[1]["MDV"]=1
        with self.assertRaisesRegex(OriginalPKCurveError,"MDV"):
            observed_curve(events)

if __name__=="__main__":unittest.main()
