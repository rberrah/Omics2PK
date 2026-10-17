import json
import unittest
from unittest.mock import patch
from pathlib import Path

from scripts.audit_rosuvastatin_author_gdown_csv import (
    author_code_drive_ids, classify_csv, public_gdown_inspect,
)

class RosuvastatinOriginalAuthorGDownTests(unittest.TestCase):
    def test_source_only_cells_extract_real_gdrive_id(self):
        source={"cells":[
          {"cell_type":"markdown","source":["https://drive.google.com/uc?id=PRIVATE_MARKDOWN_ID"]},
          {"cell_type":"code","source":["df = pd.read_csv('https://drive.google.com/uc?export=download&id=ORIGINAL_MOCK_ID_12345')"],"outputs":[{"text":"VERY_PRIVATE_PERSON_001"}]}
        ]}
        self.assertEqual(author_code_drive_ids(json.dumps(source).encode()),["ORIGINAL_MOCK_ID_12345"])

    def test_original_csv_has_40_donors_but_no_values_logged(self):
        lines=["patient_id,AUC,Cmax,feature_1,feature_2"]
        lines += [f"PERSON_SENSITIVE_{i},{10+i},{2+i},{100+i},{200+i}" for i in range(40)]
        data=("\n".join(lines)+"\n").encode()
        result=classify_csv(data)
        self.assertEqual(result["n_original_data_rows"],40)
        self.assertTrue(result["has_auc_related_column"])
        self.assertTrue(result["has_cmax_related_column"])
        self.assertTrue(result["has_subject_label_column"])
        self.assertNotIn("PERSON_SENSITIVE",repr(result))

    def test_original_gdown_source_never_leaks_cloud_id_or_observations(self):
        orig_id="SOURCE_PRIVATE_TEST_TOKEN"
        cell={"cells":[{"cell_type":"code","source":[
            "data=pd.read_csv('https://drive.google.com/uc?export=download&id="+orig_id+"')"]}]}
        def fake_download(*,id,output,quiet):
            self.assertEqual(id,orig_id)
            Path(output).write_text("patient_id,AUC,Cmax\nPERSON_PRIVATE_X,60.1,10.1\nPERSON_PRIVATE_Y,70.2,20.2\n")
            return output
        with patch("scripts.audit_rosuvastatin_author_gdown_csv.gdown.download",side_effect=fake_download):
            report=public_gdown_inspect(json.dumps(cell).encode())
        self.assertEqual(report["results"][0]["status"],"ORIGINAL_AUTHOR_CSV_RETRIEVED")
        self.assertNotIn(orig_id,str(report))
        self.assertNotIn("PERSON_PRIVATE",str(report))

    def test_google_drive_interstitial_never_misread_as_csv(self):
        with self.assertRaisesRegex(ValueError,"HTML"):
            classify_csv(b'<!doctype html><h1>Sign in to download</h1>')

if __name__=="__main__":unittest.main()
