"""Synthetic-only closed-document scanner and SHA/ZIP path tests.

The original author-provided Wiley ZIP is never used, uploaded or copied in
GitHub Actions: this CI verifies parsing guards only.
"""
from __future__ import annotations

from io import BytesIO
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from zipfile import ZipFile,ZIP_DEFLATED

from scripts.audit_local_degraeve_supplement import (
    SupplementSourceError,_docx_structure,_powerpoint_structure,
    _ensure_safe_zip,audit_original_authorized_supplement
)

W='http://schemas.openxmlformats.org/wordprocessingml/2006/main'

def zip_bytes(entries):
    b=BytesIO()
    with ZipFile(b,'w',ZIP_DEFLATED) as z:
        for name,contents in entries.items():z.writestr(name,contents)
    return b.getvalue()


class SupplementAuditTests(unittest.TestCase):
    def test_source_zip_has_exact_3_members_not_user_defined_traversal(self):
        original=zip_bytes({"2023-0506-s01.docx":b"fake","2023-0506-s02.pptx":b"fake",
                            "2023-0506-s03.docx":b"fake"})
        with _ensure_safe_zip(original,{"2023-0506-s01.docx","2023-0506-s02.pptx","2023-0506-s03.docx"}) as z:
            self.assertEqual(len(z.namelist()),3)
        with self.assertRaisesRegex(SupplementSourceError,"traversal"):
            _ensure_safe_zip(zip_bytes({"../clinical.csv":"oops"}))

    def test_reads_only_synthetic_word_tables_without_content_export(self):
        xml=f'''<w:document xmlns:w="{W}"><w:body>
          <w:p><w:r><w:t>SYNTHETIC aggregate methods</w:t></w:r></w:p>
          <w:tbl><w:tr><w:tc><w:p><w:r><w:t>Age</w:t></w:r></w:p></w:tc>
          <w:tc><w:p><w:r><w:t>38.1%</w:t></w:r></w:p></w:tc></w:tr></w:tbl>
          </w:body></w:document>'''
        b=zip_bytes({"word/document.xml":xml})
        text,tables=_docx_structure(b)
        self.assertIn("SYNTHETIC",text)
        self.assertEqual(tables,[[["Age","38.1%"]]])

    def test_counts_only_public_synthetic_graphpad_ole_objects(self):
        deck=zip_bytes({
            "ppt/slides/slide1.xml":"<p:sld/>",
            "ppt/slides/slide2.xml":"<p:sld/>",
            "ppt/embeddings/oleObject1.bin":b"Prism9.Document arbitrary synthetic bytes",
        })
        result=_powerpoint_structure(deck)
        self.assertEqual(result["n_supplementary_figure_slides"],2)
        self.assertEqual(result["n_ole_objects_recognized_as_GraphPad_Prism9"],1)
        self.assertTrue(result["graphpad_original_data_not_exported_or_decoded"])

    def test_original_source_is_sha_locked_and_never_released_from_tests(self):
        with TemporaryDirectory() as td:
            path=Path(td)/"unlicensed_source.zip"
            path.write_bytes(zip_bytes({"not_original.txt":"SYNTHETIC"}))
            with self.assertRaisesRegex(SupplementSourceError,"SHA256"):
                audit_original_authorized_supplement(path)


if __name__=="__main__":
    unittest.main()
