"""Offline source audit tests: NEVER require web or public biological data."""
import io
import json
import unittest
import zipfile
from scripts.audit_public_pk_sources import (
    SOURCES, audit, validate_payload, workbook_structure, zenodo_metadata
)


def minimal_xlsx():
    workbook = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<sheets><sheet name="Proteins" sheetId="1"/>'
        '<sheet name="Quality control" sheetId="2"/></sheets></workbook>'
    )
    destination = io.BytesIO()
    with zipfile.ZipFile(destination, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("xl/workbook.xml", workbook)
        archive.writestr("[Content_Types].xml", "<Types/>")
    return destination.getvalue()


class PrimarySourceGateTests(unittest.TestCase):
    def test_xlsx_reports_worksheet_structure_not_assay_values(self):
        status = workbook_structure(minimal_xlsx())
        self.assertEqual(status["n_sheets"], 2)
        self.assertEqual(status["sheet_names"], ["Proteins", "Quality control"])

    def test_reject_unexpected_html_block_page_instead_of_pdf_or_xlsx(self):
        for fmt in ("xlsx", "pdf"):
            with self.subTest(fmt=fmt):
                with self.assertRaises(ValueError):
                    validate_payload(fmt, b"<html>Forbidden</html>")

    def test_record_identity_fixed_not_confused_with_other_zenodo_releases(self):
        valid = {
            "id": 15723817,
            "metadata": {"title": "Test release title"},
            "files": [{"key": "data.csv"}, {"key": "reproduction.R"}],
        }
        out = zenodo_metadata(json.dumps(valid).encode())
        self.assertEqual(out["record_id"], 15723817)
        self.assertEqual(out["n_public_file_entries"], 2)
        valid["id"] = 1111111
        with self.assertRaisesRegex(ValueError, "not the cited"):
            zenodo_metadata(json.dumps(valid).encode())

    def test_blocking_conditions_preserved_even_when_all_sources_exist(self):
        fake_content = {
            "wegler_si_methods_pdf": b"%PDF-1.7\n%example",
            "wegler_s1_proteomics": minimal_xlsx(),
            "wegler_s2_bioinformatics": minimal_xlsx(),
            "yin_battail_metadata": json.dumps({
                "id": 15723817, "metadata": {"title": "Test release"}, "files": []
            }).encode(),
        }
        lookup = {u: fake_content[k] for k, (u, _) in SOURCES.items()}
        report = audit(fetch=lambda u: lookup[u])
        self.assertEqual(report["format_verified_count"], 4)
        self.assertIn("NO_PAIRED_CLINT_TARGET_CONFIRMED", report["research_status"])
        self.assertGreaterEqual(len(report["blocking_conditions"]), 4)
        self.assertNotIn("measured_clint", json.dumps(report))

    def test_unavailable_sources_are_explicit_not_fake_pass(self):
        report = audit(fetch=lambda _: b"<html>bad gateway</html>")
        self.assertEqual(report["format_verified_count"], 0)
        self.assertTrue(all(
            item["state"] == "UNAVAILABLE_OR_INVALID"
            for item in report["sources"].values()
        ))


if __name__ == "__main__":
    unittest.main()
