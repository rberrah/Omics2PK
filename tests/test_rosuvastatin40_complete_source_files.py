import json
import unittest
from scripts.audit_rosuvastatin40_complete_original_filelist import scan_original_public_list

class OriginalMetaboLightsFullManifestTests(unittest.TestCase):
    def test_detects_only_metadata_and_raw(self):
        src={"study":[
            {"file":"a_MTBLS3128_positive.txt","type":"assay"},
            {"file":"FMTB_ROS_PRIVATE_DONOR_001_PRE.mzML","type":"raw"},
            {"file":"ESI_positive_processed_quant.csv","type":"derived"},
            {"file":"human_quant_PK_auc.csv","type":"derived"},
        ]}
        report=scan_original_public_list(json.dumps(src).encode())
        self.assertEqual(report["n_raw_MS_files"],1)
        self.assertEqual(report["n_possible_processed_structured_files"],2)
        self.assertTrue(report["possible_processed_structured_files"][1]["filename_has_pk_keyword"])
        self.assertNotIn("PRIVATE_DONOR_001",repr(report))
        self.assertNotIn("ESI_positive_processed_quant.csv",repr(report))
        self.assertFalse(report["source_sample_IDs_or_PK_values_exported"])

    def test_full_schema_identity_required(self):
        with self.assertRaisesRegex(ValueError,"schema missing"):
            scan_original_public_list(b'{"data":[{"file":"foo"}]}')

    def test_does_not_follow_unsafe_sample_path(self):
        with self.assertRaisesRegex(ValueError,"Unsafe"):
            scan_original_public_list(b'{"study":[{"file":"../../original_people.csv"}]}')


if __name__=="__main__":
    unittest.main()
