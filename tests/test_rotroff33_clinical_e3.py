"""Synthetic-only E3 numerical model and source leakage guards: NO human data."""
from __future__ import annotations
import importlib.util
import unittest

from omics2pk.rotroff33_clinical_e3 import (
    E3SourceError, _report_one_target, evaluate_e3_original
)

def synthetic():
    rows=[]
    for i in range(33):
        age=22.+i%24
        weight=55.+(i*7)%27
        male=float(i%2)
        omic=[float(((i+4)*(j+7))%77)/11. for j in range(360)]
        omic[10]=None if i%3==0 else omic[10]
        auc=1000.+5*age+3*weight+(50 if male else 0)+4*omic[4]
        clr=250.+age*0.75+2*omic[19]
        if i==32:auc=clr=None
        rows.append({"internal_patient_id":f"FAKE_N32_{i:02d}",
                     "clinical":[age,male,weight],
                     "predose_omics":omic,
                     "targets":{"metformin_AUClast_0_24_original_units":auc,
                                "metformin_CLr_renal_clearance_original_units":clr}})
    return rows

@unittest.skipUnless(importlib.util.find_spec("numpy") and
                     importlib.util.find_spec("sklearn"),"optional scientific packages not installed")
class RealPatientSourceProtocolSyntheticTests(unittest.TestCase):
    def test_exploratory_heldout_pca_foldwise_no_output_subject_ids(self):
        records=synthetic()
        r=evaluate_e3_original(records,{"synthetic_only":True})
        self.assertEqual(r["n_original_related_measured_pk_endpoints"],2)
        for result in r["endpoint_results"].values():
            self.assertEqual(result["n_original_PK_measured_subjects"],32)
            self.assertEqual(result["n_original_subject_heldout_folds"],32)
            self.assertEqual(result["n_original_source_omics_features_considered"],360)
            self.assertTrue(result["no_subject_IDs_numeric_or_predictions_exported"])
            self.assertIn("clinical_plus_3_predose_omics_PCs",result["performance"])
        self.assertNotIn("FAKE_N32",str(r))
        self.assertFalse(r["original_subject_IDs_and_medical_measurements_exported"])

    def test_unauthorized_target_count_stops_computation(self):
        records=synthetic()
        records[31]["targets"]["metformin_AUClast_0_24_original_units"]=None
        with self.assertRaisesRegex(E3SourceError,"32"):
            _report_one_target(records,"metformin_AUClast_0_24_original_units")


if __name__=="__main__":unittest.main()
