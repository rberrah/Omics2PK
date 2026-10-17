import io
import unittest
from zipfile import ZipFile

from scripts.audit_jarmusch14_europepmc_supplement_zip import inspect_archive


class EuropeanOriginalSourceNoLeakTests(unittest.TestCase):
    def test_safe_nested_original_archive(self):
        with io.BytesIO() as inner:
            with ZipFile(inner,"w") as z:
                z.writestr("study.pdf",b"pdf aggregate only")
            with io.BytesIO() as outer:
                with ZipFile(outer,"w") as z:
                    z.writestr("original-study-supp.zip",inner.getvalue())
                o=inspect_archive(outer.getvalue())
        self.assertEqual(o[0]["nested_items"][0]["file"],"study.pdf")
        self.assertNotIn("subject_id",repr(o))

    def test_archive_path_traversal_rejected(self):
        with io.BytesIO() as src:
            with ZipFile(src,"w") as z:
                z.writestr("../secret.csv",b"private patient records")
            with self.assertRaisesRegex(ValueError,"Unsafe"):
                inspect_archive(src.getvalue())

    def test_not_original_zip_rejected(self):
        with self.assertRaisesRegex(ValueError,"Invalid"):
            inspect_archive(b"<html>browser robot</html>")


if __name__=="__main__":unittest.main()
