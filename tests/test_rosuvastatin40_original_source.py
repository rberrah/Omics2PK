"""Original rosuvastatin MTBLS3128: identity, source provenance, no volunteer rows."""
import io
import json
import unittest
from zipfile import ZipFile

from scripts.audit_rosuvastatin40_original_mtbls3128 import json_structure,zip_structure


class OriginalRosuvastatin40Test(unittest.TestCase):
    def test_real_accession_and_no_sample_identifiers_exported(self):
        src={"accession":"MTBLS3128","files":[{"fileName":"ORIGINAL_VOLUNTEER_A_PRE.mzML"}],
             "samples":[{"sampleId":"PERSON_PRIVATE","subject_id":"PERSON_PRIVATE"}]}
        got=json_structure(json.dumps(src).encode())
        self.assertTrue(got["study_ID_checked"])
        self.assertEqual(got["file_extensions_or_metadata_collection_counts"]["key_samples_entries"],1)
        self.assertNotIn("ORIGINAL_VOLUNTEER_A",str(got))
        self.assertNotIn("PERSON_PRIVATE",str(got))

    def test_other_dataset_rejected(self):
        with self.assertRaisesRegex(ValueError,"identity mismatch"):
            json_structure(b'{"accession":"MTBLS9999","files":[]}')

    def test_unsafe_zip_rejected(self):
        with io.BytesIO() as bio:
            with ZipFile(bio,"w") as z:
                z.writestr("../secrets.csv","protected")
            with self.assertRaisesRegex(ValueError,"Unsafe"):
                zip_structure(bio.getvalue())


if __name__=="__main__":
    unittest.main()
