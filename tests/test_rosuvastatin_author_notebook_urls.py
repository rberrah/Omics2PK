import json
import unittest
from scripts.audit_rosuvastatin_author_notebook_data_urls import classify_notebook_payload

class AuthorNotebookLinkTests(unittest.TestCase):
    def test_output_data_not_exposed_and_opaque_cloud_links_masked(self):
        doc={"cells":[
            {"cell_type":"code","source":[
                "import pandas as pd\n",
                "data=pd.read_csv('https://raw.githubusercontent.com/TEST/test/main/processed.csv')\n",
                "private='https://drive.google.com/file/d/PSEUDOPRIVATE123/view'\n"
            ],"outputs":[{"text":"PRIVATE_VOLUNTEER_ABC;AUC=33.11"}]},
            {"cell_type":"markdown","source":["https://example.org/hello"]}]}
        out=classify_notebook_payload(json.dumps(doc).encode())
        self.assertEqual(out["n_masked_cloud_drive_references"],1)
        self.assertEqual(out["n_cells_with_uninspected_embedded_outputs"],1)
        self.assertEqual(out["public_source_data_url_paths"],
                         ["https://raw.githubusercontent.com/TEST/test/main/processed.csv"])
        self.assertNotIn("PSEUDOPRIVATE123",str(out))
        self.assertNotIn("PRIVATE_VOLUNTEER_ABC",str(out))
        self.assertNotIn("33.11",str(out))

    def test_wrong_nonnotebook_source_rejected(self):
        with self.assertRaisesRegex(ValueError,"Not a Jupyter"):
            classify_notebook_payload(b'{"cells":"evil"}')

    def test_no_url_does_not_create_omics_dataset(self):
        result=classify_notebook_payload(b'{"cells":[{"cell_type":"code","source":["import numpy as np"]}]}')
        self.assertEqual(result["public_source_data_url_paths"],[])
        self.assertTrue(result["no_original_notebook_outputs_or_patient_data_exported"])

if __name__=="__main__":unittest.main()
