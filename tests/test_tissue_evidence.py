"""Fail-closed source tissular and clinical landmark audit: no patient records."""
import unittest
from omics2pk.tissue_evidence import audit_tissue_source_provenance,TissueEvidenceError
class TissueEvidenceTests(unittest.TestCase):
 def test_blood_to_liver_is_candidate_not_direct_function(self):
  rows=[
   {"tissue":"blood","specimen_type":"PBMC","cell_type_or_bulk":"bulk",
    "cell_composition_status":"unknown","pharmacology_target_organ":"liver","relative_to_dose":"pre_dose"},
   {"tissue":"blood","specimen_type":"plasma","cell_type_or_bulk":"cell_free",
    "pharmacology_target_organ":"liver","relative_to_dose":"post_dose"},
   {"tissue":"liver","specimen_type":"liver_biopsy","cell_type_or_bulk":"bulk",
    "cell_composition_status":"estimated","pharmacology_target_organ":"liver","relative_to_dose":"not_applicable"},
  ]
  r=audit_tissue_source_provenance(rows)
  self.assertEqual(r["n_original_tissue_declared"],3)
  self.assertEqual(r["n_blood_vs_hepatic_or_intestinal_target_candidates_unvalidated"],2)
  self.assertEqual(r["n_original_predose_declared"],1)
  self.assertEqual(r["n_original_postdose_declared"],1)
  self.assertEqual(r["n_bulk_assays_without_measured_or_estimated_cell_composition"],1)
  self.assertFalse(r["whole_blood_omics_can_be_qualified_direct_hepatic_CYP3A_activity"])
  self.assertTrue(r["statistical_use_of_blood_as_indirect_Pk_predictor_is_not_prohibited"])
  self.assertTrue(r["source_omics_measurements_and_patient_ids_not_returned"])
 def test_missing_source_is_not_assumed_hepatic(self):
  r=audit_tissue_source_provenance([{"assay_id":"synthetic","subject_id":"never_echo"}])
  self.assertEqual(r["n_original_tissue_declared"],0)
  self.assertEqual(r["n_original_predose_declared"],0)
  self.assertNotIn("subject_id",str(r))
  self.assertEqual(r["status"],"SOURCE_TISSUE_AND_SPECIMEN_NOT_FULLY_DECLARED")
 def test_T0_does_not_mean_predose(self):
  with self.assertRaises(TissueEvidenceError):
   audit_tissue_source_provenance([{"relative_to_dose":"T0"}])
 def test_unsupported_cell_mode_requires_protocol_mapping(self):
  with self.assertRaises(TissueEvidenceError):
   audit_tissue_source_provenance([{"cell_type_or_bulk":"hepatocyte"}])
if __name__=="__main__":unittest.main()
