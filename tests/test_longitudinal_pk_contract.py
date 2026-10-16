"""Synthetic PK_LONG_V1 regression safety: whole dose history, repeated molecular state."""
import unittest
from omics2pk.longitudinal_pk_contract import (
    LongitudinalPKError, qualify_longitudinal_pk, align_molecular_states_to_landmarks
)

def study():
    return [
        dict(subject_id="SYN-A",drug="tacrolimus",time_h=0,event_type="dose",data_origin="observed",dose_id="D1",
             amount_mg=2,route="oral",formulation="IR"),
        dict(subject_id="SYN-A",drug="tacrolimus",time_h=0,event_type="observation",data_origin="observed",
             sample_id="C0",concentration_mg_l=0.005,analyte="tacrolimus",specimen="whole_blood",sampling_phase="predose"),
        dict(subject_id="SYN-A",drug="tacrolimus",time_h=11.8,event_type="observation",data_origin="observed",
             sample_id="C1",concentration_mg_l=0.008,analyte="tacrolimus",specimen="whole_blood"),
        dict(subject_id="SYN-A",drug="tacrolimus",time_h=12,event_type="dose",data_origin="observed",dose_id="D2",
             amount_mg=2,route="oral",formulation="IR"),
        dict(subject_id="SYN-A",drug="tacrolimus",time_h=24,event_type="observation",data_origin="observed",
             sample_id="C2",concentration_mg_l=None,is_blq=True,lloq_mg_l=0.001,analyte="tacrolimus",
             specimen="whole_blood",sampling_phase="predose"),
        dict(subject_id="SYN-A",drug="tacrolimus",time_h=24,event_type="dose",data_origin="simulated",dose_id="D3",
             amount_mg=2,route="oral",formulation="IR"),
        dict(subject_id="SYN-A",drug="tacrolimus",time_h=36,event_type="observation",data_origin="simulated",
             sample_id="C3",concentration_mg_l=0.006,analyte="tacrolimus",specimen="whole_blood"),
    ]


class PKLongitudinalTests(unittest.TestCase):
    def test_repeated_doses_and_measured_vs_simulated_never_mixed_in_counts(self):
        report=qualify_longitudinal_pk(study(),drug="tacrolimus",
          time_origin="study_elapsed_h",dataset_kind="multiple_dose",dosing_history_status="complete")
        self.assertEqual(report["n_independent_subjects"],1)
        self.assertEqual(report["n_dose_events"],3)
        self.assertEqual(report["n_observed_concentrations"],3)
        self.assertEqual(report["n_simulated_concentrations"],1)
        self.assertEqual(report["n_below_quantification_limit"],1)
        self.assertEqual(report["n_predose_observations"],2)
        self.assertNotIn("SYN-A",repr(report))

    def test_same_time_concentration_needs_predose_order(self):
        rows=study()
        rows[1]["sampling_phase"]="postdose"
        with self.assertRaisesRegex(LongitudinalPKError,"Time-zero"):
            qualify_longitudinal_pk(rows,drug="tacrolimus",
               time_origin="t0",dataset_kind="multiple_dose",dosing_history_status="complete")

    def test_blq_is_not_zero(self):
        rows=study()
        rows[4]["concentration_mg_l"]=0.0
        with self.assertRaisesRegex(LongitudinalPKError,"BLQ must not"):
            qualify_longitudinal_pk(rows,drug="tacrolimus",
               time_origin="t0",dataset_kind="multiple_dose",dosing_history_status="complete")

    def test_dose_event_must_not_contain_concentration(self):
        rows=study()
        rows[0]["concentration_mg_l"]=0.0
        with self.assertRaisesRegex(LongitudinalPKError,"may not contain"):
            qualify_longitudinal_pk(rows,drug="tacrolimus",
               time_origin="t0",dataset_kind="multiple_dose",dosing_history_status="complete")

    def test_unknown_dosing_history_never_autocompletes(self):
        report=qualify_longitudinal_pk(study(),drug="tacrolimus",
          time_origin="t0",dataset_kind="multiple_dose",dosing_history_status="unknown")
        self.assertFalse(report["dosing_history_suitable_for_mechanistic_fit"])

    def test_single_dose_rejects_repeated(self):
        with self.assertRaisesRegex(LongitudinalPKError,"Single-dose"):
            qualify_longitudinal_pk(study(),drug="tacrolimus",
               time_origin="t0",dataset_kind="single_dose",dosing_history_status="complete")

    def test_missing_observed_dose_rejected(self):
        rows=study()[1:]
        with self.assertRaises(LongitudinalPKError):
            qualify_longitudinal_pk(rows,drug="tacrolimus",
               time_origin="t0",dataset_kind="multiple_dose",dosing_history_status="complete")

    def test_repeated_exposure_and_omics_timeline_allow_postdose_but_before_future_landmark(self):
        assays=[
            dict(subject_id="SYN-A",sample_id="O-BASE",modality="RNA",time_h=0),
            dict(subject_id="SYN-A",sample_id="O-FOLLOWUP",modality="RNA",time_h=16),
            dict(subject_id="SYN-A",sample_id="O-FUTURE",modality="RNA",time_h=60),
        ]
        landmarks=[dict(subject_id="SYN-A",time_h=t) for t in (12,24,72)]
        qc=align_molecular_states_to_landmarks(assays,landmarks,modality="RNA",max_lookback_h=24)
        self.assertEqual(qc["n_prediction_landmarks"],3)
        self.assertEqual(qc["n_landmarks_with_preceding_molecular_measurement"],3)
        self.assertGreater(qc["n_future_omic_samples_excluded_by_landmark"],0)
        self.assertTrue(qc["no_future_molecular_state_used"])
        self.assertNotIn("SYN-A",repr(qc))

    def test_simultaneous_molecular_sample_is_ambiguous(self):
        assays=[dict(subject_id="SYN-A",sample_id="RNA0",modality="RNA",time_h=12),
                dict(subject_id="SYN-A",sample_id="RNA1",modality="RNA",time_h=10)]
        qc=align_molecular_states_to_landmarks(assays,[dict(subject_id="SYN-A",time_h=12)],
          modality="RNA",max_lookback_h=24)
        self.assertEqual(qc["n_landmarks_with_ambiguous_simultaneous_assay"],1)

    def test_no_overclaim_without_real_measured_samples(self):
        rows=study()
        for row in rows:row["data_origin"]="simulated"
        with self.assertRaisesRegex(LongitudinalPKError,"genuinely measured"):
            qualify_longitudinal_pk(rows,drug="tacrolimus",
             time_origin="t0",dataset_kind="multiple_dose",dosing_history_status="complete")

if __name__=="__main__":
    unittest.main()
