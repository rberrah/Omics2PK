"""All synthetic tests; no original EBI spectra or patient data bundled."""
from __future__ import annotations

from io import BytesIO
from unittest.mock import patch
from zipfile import ZipFile, ZIP_DEFLATED
import unittest

from scripts.audit_rosuvastatin40_zip_member_formats import (
    TAIL_BYTES, SourceZipError, archive_member_extensions,
    zip_end_metadata, probe_zip_central
)


def synthetic_zip():
    buf=BytesIO()
    with ZipFile(buf,"w",ZIP_DEFLATED) as archive:
        archive.writestr("private_label/raw-spectrum.raw","NOT_REAL")
        archive.writestr("private_label/summary.csv","NO PATIENTS")
        archive.writestr("private_label/metadata.xml","FAKE")
    return buf.getvalue()


class EBIOriginalZipFormatTests(unittest.TestCase):
    def test_zip_footer_and_central_directory_extract_formats_without_names(self):
        data=synthetic_zip()
        offset,length,count=zip_end_metadata(data[-TAIL_BYTES:],file_size=len(data))
        report=archive_member_extensions(data[offset:offset+length],expected_count=count)
        self.assertEqual(report["n_original_archive_members"],3)
        self.assertEqual(report["n_candidate_tabular_data_members"],1)
        self.assertEqual(report["original_member_extension_counts"][".csv"],1)
        self.assertNotIn("private_label",str(report))

    def test_probes_only_bounded_ranges_and_leaks_no_original_member_names(self):
        source=synthetic_zip()
        calls=[]
        def fake_range(url,spec,*,maximum):
            calls.append(spec)
            if spec.startswith("-"):
                return source[-int(spec[1:]):]
            begin,end=map(int,spec.split("-"))
            return source[begin:end+1]
        with patch("scripts.audit_rosuvastatin40_zip_member_formats._range",side_effect=fake_range):
            report=probe_zip_central("https://ftp.ebi.ac.uk/FAKE.zip",len(source))
        self.assertEqual(len(calls),2)
        self.assertNotIn("FAKE.zip",str(report))
        self.assertNotIn("private_label",str(report))
        self.assertEqual(report["n_candidate_tabular_data_members"],1)

    def test_reject_nonzip_and_truncated_central_directory(self):
        with self.assertRaises(SourceZipError):
            zip_end_metadata(b"INVALID",file_size=1234)
        with self.assertRaises(SourceZipError):
            archive_member_extensions(b"PK\x01\x02",expected_count=1)

    def test_reject_implausible_file_sizes(self):
        with self.assertRaises(SourceZipError):
            probe_zip_central("https://example.org/empty.zip",1)


if __name__=="__main__":
    unittest.main()
