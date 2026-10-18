"""Synthetic-only tests for direct anonymous EBI public MTBLS3128 retrieval."""
from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import json
import tempfile
import unittest

from scripts.download_mtbls3128_original_public import (
    ROOT, DownloadGateError, _safe_relative, download_study, select_files,
    source_inventory,
)


class FakeFTP:
    def __init__(self):
        self.files = {
            f"{ROOT}/i_Investigation.txt": b"INVESTIGATION\tSYNTHETIC\n",
            f"{ROOT}/s_source.txt": b"Sample Name\tSource Name\nSYNTH_A\tSYNTH_A\n",
            f"{ROOT}/a_positive.txt": b"Sample Name\tSource Name\nSYNTH_A\tSYNTH_A\n",
            f"{ROOT}/a_negative.txt": b"Sample Name\tSource Name\nSYNTH_A\tSYNTH_A\n",
            f"{ROOT}/m_positive.tsv": b"Metabolite Identification\tSYNTH_A\n",
            f"{ROOT}/m_negative.tsv": b"Metabolite Identification\tSYNTH_A\n",
            f"{ROOT}/FILES/ESI_POS/SYNTH_A.mzML": b"<mzML>synthetic spectrum</mzML>",
            f"{ROOT}/FILES/ESI_NEG/SYNTH_B.mzML": b"<mzML>synthetic spectrum</mzML>",
            f"{ROOT}/FILES/RAW/SYNTH_A.zip": b"PK\x03\x04SYNTHETIC",
        }
        self.directories = [
            ROOT, ROOT + "/FILES",
            ROOT + "/FILES/ESI_POS",
            ROOT + "/FILES/ESI_NEG",
            ROOT + "/FILES/RAW",
        ]

    def mlsd(self, directory):
        base = directory.rstrip("/") + "/"
        result = []
        for subdir in self.directories:
            if subdir != directory and subdir.startswith(base):
                rest = subdir[len(base):]
                if "/" not in rest:
                    result.append((rest, {"type": "dir"}))
        for path in self.files:
            if path.startswith(base):
                rest = path[len(base):]
                if "/" not in rest:
                    result.append((rest, {"type": "file"}))
        return result

    def retrbinary(self, cmd, callback, blocksize=131072):
        assert cmd.startswith("RETR ")
        data = self.files[cmd[5:]]
        for start in range(0, len(data), 8):
            callback(data[start:start + 8])


class LocalSourceTests(unittest.TestCase):
    def test_reject_outside_study_paths(self):
        for key in (ROOT + "/../MISMATCH.txt", ROOT + "/FILES/..\\evil",
                    "/tmp/other_study.txt"):
            with self.assertRaises(DownloadGateError):
                _safe_relative(key)

    def test_public_metadata_and_raw_list_separate(self):
        ftp = FakeFTP()
        meta, mzml, zipped = source_inventory(ftp, include_raw=False)
        self.assertEqual(len(meta), 6)
        self.assertEqual(mzml, [])
        self.assertEqual(zipped, [])
        meta, mzml, zipped = source_inventory(ftp, include_raw=True)
        self.assertEqual((len(meta), len(mzml), len(zipped)), (6, 2, 1))
        self.assertEqual(len(select_files(meta, mzml, zipped, mzml_max_files=1, zip_max_files=0)), 7)

    def test_real_bytes_copy_with_sha_local_manifest_and_no_stdout_ids(self):
        fake = FakeFTP()
        with tempfile.TemporaryDirectory() as td:
            private = Path(td) / "private" / "MTBLS3128"
            outcome = download_study(
                fake, output_dir=private, mzml_max_files=2,
                max_total_bytes=1024 * 1024, max_file_bytes=1024 * 1024,
            )
            self.assertEqual(outcome["downloaded_file_count"], 8)
            self.assertFalse(outcome["clinical_E3_results_qualified"] if "clinical_E3_results_qualified" in outcome else False)
            self.assertNotIn("SYNTH_A", repr(outcome))
            self.assertNotIn("SYNTH_B", repr(outcome))
            m = json.loads((private / "private_source_manifest.json").read_text())
            self.assertEqual(len(m["original_source_manifest"]), 8)
            for entry in m["original_source_manifest"]:
                copied = private / entry["relative_original_filename"]
                self.assertEqual(copied.read_bytes(), fake.files[ROOT + "/" + entry["relative_original_filename"]])
                self.assertEqual(entry["original_download_sha256"], sha256(copied.read_bytes()).hexdigest())

    def test_metadata_only_default_and_preview_no_materialization(self):
        with tempfile.TemporaryDirectory() as td:
            private = Path(td) / "private" / "MTBLS3128"
            preview = download_study(FakeFTP(), output_dir=private, list_only=True)
            self.assertEqual(preview["selected_mzml_files"], 0)
            self.assertFalse(private.exists())
            got = download_study(FakeFTP(), output_dir=private, max_total_bytes=1024)
            self.assertEqual(got["downloaded_file_count"], 6)
            self.assertFalse((private / "FILES").exists())

    def test_existing_source_cannot_be_silently_overwritten(self):
        with tempfile.TemporaryDirectory() as td:
            private = Path(td) / "private" / "MTBLS3128"
            download_study(FakeFTP(), output_dir=private, max_total_bytes=1024)
            with self.assertRaisesRegex(DownloadGateError, "exists"):
                download_study(FakeFTP(), output_dir=private, max_total_bytes=1024)

    def test_limits_cannot_be_exceeded_or_leave_partial_file(self):
        with tempfile.TemporaryDirectory() as td:
            private = Path(td) / "private" / "MTBLS3128"
            with self.assertRaisesRegex(DownloadGateError, "maximum bytes"):
                download_study(FakeFTP(), output_dir=private, max_total_bytes=1024, max_file_bytes=5)
            self.assertFalse(any(private.rglob(".partial-*")))
            self.assertFalse((private / "private_source_manifest.json").exists())

    def test_refuse_more_raw_files_than_original_inventory(self):
        with self.assertRaisesRegex(DownloadGateError, "more original files"):
            download_study(FakeFTP(), output_dir=Path("/tmp/NO_WRITE_NEEDED"),
                           mzml_max_files=99, list_only=True)


if __name__ == "__main__":
    unittest.main()
