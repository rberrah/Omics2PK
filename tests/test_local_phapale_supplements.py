"""No copyrighted Wiley source bytes in public unit tests."""
import io
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from zipfile import ZipFile,ZIP_DEFLATED
from scripts.audit_local_phapale_supplements import (
    OriginalPhapaleSourceError,EXPECTED_FILES,OLE_HEADER,
    inspect_archive_structure,inspect_original_local
)

def _source(entries):
    out=io.BytesIO()
    with ZipFile(out,"w",ZIP_DEFLATED) as z:
        z.writestr("paphale/",b"")
        for k,v in entries.items():z.writestr(k,v)
    return out.getvalue()


def synthetic_13():
    return {f"paphale/cptclpt2009296-sup-{n}.doc":OLE_HEADER+b"THIS_IS_SYNTHETIC_ONLY"
            for n in EXPECTED_FILES}


class OriginalPhapaleSupplementTests(unittest.TestCase):
    def test_exact_original_supplement_number_without_releasing_source_values(self):
        qc=inspect_archive_structure(_source(synthetic_13()),require_original_sha=False)
        self.assertEqual(qc["n_genuine_Wiley_legacy_Word_supplements"],13)
        self.assertTrue(qc["all_original_legacy_Word_ole2_headers_verified"])
        self.assertFalse(qc["patient_PK_Ct_data_released"])

    def test_no_faking_patient_matrix_by_named_figure(self):
        entries=synthetic_13()
        entries.pop("paphale/cptclpt2009296-sup-0013.doc")
        with self.assertRaisesRegex(OriginalPhapaleSourceError,"13 original"):
            inspect_archive_structure(_source(entries),require_original_sha=False)

    def test_legacy_source_member_magic_required(self):
        entries=synthetic_13()
        entries["paphale/cptclpt2009296-sup-0011.doc"]=b"not a word document"
        with self.assertRaisesRegex(OriginalPhapaleSourceError,"OLE2"):
            inspect_archive_structure(_source(entries),require_original_sha=False)

    def test_zip_traversal_is_refused(self):
        entries=synthetic_13()
        entries["../PRIVATE.xlsx"]=b"never"
        with self.assertRaisesRegex(OriginalPhapaleSourceError,"Unsafe"):
            inspect_archive_structure(_source(entries),require_original_sha=False)

    def test_uploaded_original_sha_must_match_before_any_local_source_parsing(self):
        with TemporaryDirectory() as folder:
            p=Path(folder)/"fake.zip"
            p.write_bytes(_source(synthetic_13()))
            with self.assertRaisesRegex(OriginalPhapaleSourceError,"SHA-256 mismatch"):
                inspect_original_local(p)


if __name__=="__main__":
    unittest.main()
