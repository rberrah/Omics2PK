import unittest
from scripts.audit_rosuvastatin40_isa_patient_crosswalk import find_metadata_filenames,summarize_isa


class RosuvastatinISATest(unittest.TestCase):
    def test_exact_original_study_listing(self):
        p=b'{"study":[{"file":"s_MTBLS3128.txt","type":"metadata_sample"},{"file":"a_MTBLS3128_lcms.txt","type":"metadata_assay"},{"file":"../../pseudonyms.txt"}]}'
        with self.assertRaisesRegex(ValueError,"path"):
            find_metadata_filenames(p)

    def test_safe_metadata_and_no_volunteer_rows(self):
        s="Sample Name\tSource Name\tParameter Value[Time]\nPERSONPRIVATE_A\tPERSONPRIVATE_A\tpre-dose\nPERSONPRIVATE_B\tPERSONPRIVATE_B\tpre-dose\n"
        r,ids=summarize_isa(s.encode(),"metadata_sample")
        self.assertEqual(r["n_original_metadata_rows"],2)
        self.assertEqual(r["n_unique_source_names"],2)
        self.assertTrue(r["has_predose_factor"])
        self.assertNotIn("PERSONPRIVATE_A",str(r))
        self.assertTrue(r["sample_values_exported"] is False)
        self.assertEqual(len(ids),2)

    def test_pk_fields_are_not_inferred_without_real_columns(self):
        r,_=summarize_isa(b"Sample Name\tSome Clinical Descriptor\nPRIVATE_A\tPK baseline\n","sample")
        self.assertFalse(r["has_pk_auc_column"])
        self.assertFalse(r["has_pk_cmax_column"])
        self.assertFalse(r["has_drug_concentration_column"])

if __name__=="__main__":unittest.main()
