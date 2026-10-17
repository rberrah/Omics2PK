import unittest
from scripts.audit_rosuvastatin40_public_ftp_census import summarize_source_names

class OriginalFTPManifestPrivacyTests(unittest.TestCase):
    def test_only_counts_and_roles_no_original_sample_filenames(self):
        items=["secret_subject_01/PRE_original_0001.mzML",
               "SECRET_SUBJECT_2/features_processed_aligned.csv",
               "OTHER/PK_original_auc_cmax.xlsx",
               "secret_subject_03/PRE_original_0003.mzML"]
        result=summarize_source_names(items)
        self.assertEqual(result["original_listed_data_file_count"],4)
        self.assertEqual(result["n_raw_MS_format_files"],2)
        self.assertEqual(result["structured_candidate_counts"]["structured_file"],2)
        self.assertEqual(result["structured_candidate_counts"]["structured_name_suggests_PK_target"],1)
        self.assertNotIn("secret_subject",repr(result).lower())
        self.assertNotIn("PK_original_auc_cmax.xlsx",repr(result))
        self.assertFalse(result["raw_spectra_downloaded"])

    def test_raw_only_means_no_candidate_processed_matrix(self):
        result=summarize_source_names(["PERSON1/a.raw","PERSON2/b.mzML"])
        self.assertEqual(result["n_raw_MS_format_files"],2)
        self.assertEqual(result["structured_candidate_counts"],{})

if __name__=="__main__":
    unittest.main()
