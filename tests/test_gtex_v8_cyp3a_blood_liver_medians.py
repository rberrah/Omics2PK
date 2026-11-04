import unittest,gzip,io,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_gtex_v8_cyp3a_blood_liver_medians import extract_original_gct,GTExSourceError
class GtexTissueSourceTests(unittest.TestCase):
    def test_not_accept_fake_tiny_source(self):
        b=gzip.compress(("#1.2\n50000\t54\nName\tDescription\tLiver\tWhole_Blood\tSmall_Intestine_Terminal_Ileum\n"
                         "ENSGX\tCYP3A4\t100\t1\t10\n").encode())
        with self.assertRaises(GTExSourceError):
            extract_original_gct(io.BytesIO(b),"synthetic")
    def test_missing_true_tissue_columns_fails_closed(self):
        b=gzip.compress(("#1.2\n50000\t54\nName\tDescription\tLiver\tPBMC\n"
                         "ENSGX\tCYP3A4\t100\t1\n").encode())
        with self.assertRaises(GTExSourceError):
            extract_original_gct(io.BytesIO(b),"synthetic")
if __name__=="__main__":unittest.main()
