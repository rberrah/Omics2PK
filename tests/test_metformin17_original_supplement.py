"""Fail-closed MDPI/PMC original supplemental ZIP inspection: no human row leakage."""
import io
from zipfile import ZipFile
import unittest
from scripts.audit_metformin17_441k_original_supplement import members_no_rows


class Metformin17SupplementTests(unittest.TestCase):
    def test_original_zip_member_dimensions_no_donor_records(self):
        with io.BytesIO() as bio:
            with ZipFile(bio,"w") as arch:
                arch.writestr("PK_Table_S5.csv","id,drug_concentration\nPRIVATE_SUBJECT_A,999.1\n")
            results=members_no_rows(bio.getvalue())
        self.assertEqual(results[0]["filename"],"PK_Table_S5.csv")
        self.assertEqual(results[0]["csv_line_count_approx"],2)
        self.assertNotIn("PRIVATE_SUBJECT_A",repr(results))
        self.assertNotIn("999.1",repr(results))

    def test_wrong_content_fails_as_source(self):
        with self.assertRaisesRegex(ValueError,"Not one"):
            members_no_rows(b"<html>publisher interstitial</html>")

    def test_path_traversal_archive_is_rejected(self):
        with io.BytesIO() as bio:
            with ZipFile(bio,"w") as arch:
                arch.writestr("../outside.txt","not permitted")
            with self.assertRaisesRegex(ValueError,"Unsafe"):
                members_no_rows(bio.getvalue())


if __name__=="__main__":
    unittest.main()
