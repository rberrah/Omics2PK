"""Synthetic public BinDiscover schema tests, never a fake molecular match."""
import json
import unittest
from scripts.probe_bindiscover_rotroff_shap_bins import summarize_source_reply,BINS
class BinSourceProbeTests(unittest.TestCase):
 def test_documented_public_reply_is_not_original_cohort_identity_proof(self):
  mock={"bin_id":222169,"retentionIndex":244189.0,
        "spectrum":"85:0.1 86:0.2 87:0.3 100:0.4 110:0.5 120:0.6 130:0.7"}
  out=summarize_source_reply(json.dumps(mock).encode(),"222169")
  self.assertEqual(out["n_spectral_ions_in_public_response"],7)
  self.assertTrue(out["documented_bin_id_explicitly_echoed_and_matches"])
  self.assertFalse(out["original_Rotroff_bin_namespace_correspondence_proved"])
  self.assertFalse(out["chemical_structure_identified"])
  self.assertFalse(out["spectrum_raw_values_exported"])
  self.assertNotIn("spectrum",out)
 def test_absent_spectrum_does_not_invent_molecular_mass(self):
  out=summarize_source_reply(json.dumps({"retentionIndex":244189.0}).encode(),"222169")
  self.assertIsNone(out["n_spectral_ions_in_public_response"])
  self.assertIsNone(out["documented_bin_id_explicitly_echoed_and_matches"])
  self.assertNotIn("nominal_mass",out)
 def test_unsupported_payload_fails_closed(self):
  with self.assertRaises(ValueError):
   summarize_source_reply(b"<html>login</html>","222169")
  with self.assertRaises(ValueError):
   summarize_source_reply(json.dumps([{"x":[1]},{"x":[2]}]).encode(),"222169")
 def test_top_ids_are_original_selected(self):
  self.assertEqual(BINS,("222169","221574","218765"))
if __name__=="__main__":unittest.main()
