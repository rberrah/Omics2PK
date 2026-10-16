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
        valid["id"] = 15723818
        valid["doi"] = "10.5281/zenodo.15723818"
        valid["conceptdoi"] = "10.5281/zenodo.15723817"
        version = zenodo_metadata(json.dumps(valid).encode())
        self.assertEqual(version["record_id"], 15723818)
        self.assertEqual(version["cited_concept_record"], 15723817)
        valid["id"] = 1111111
        valid["conceptdoi"] = "10.5281/zenodo.OTHER"
        with self.assertRaisesRegex(ValueError, "Zenodo record mismatch"):
            zenodo_metadata(json.dumps(valid).encode())

    def test_blocking_conditions_preserved_even_when_all_sources_exist(self):
        fake_content = {
            "wegler_si_methods_pdf": b"%PDF-1.7\n%example",
            "wegler_s1_proteomics": minimal_xlsx(),
            "wegler_s2_bioinformatics": minimal_xlsx(),
            "wegler_cloud_metadata": json.dumps({"pmcid":"PMC8041379",
                "doi":"10.1021/acs.molpharmaceut.1c00053", "version":1,
                "license_code":"CC BY 4.0",
                "media_urls":["https://pmc-oa-opendata.s3.amazonaws.com/PMC8041379.1/mp1c00053_si_002.xlsx"]}).encode(),
            "yin_battail_metadata": json.dumps({
                "id": 15723817, "metadata": {"title": "Test release"}, "files": []
            }).encode(),
        }
        lookup = {u: fake_content[k] for k, (u, _) in SOURCES.items()}
        lookup["https://pmc-oa-opendata.s3.amazonaws.com/PMC8041379.1/mp1c00053_si_002.xlsx"] = minimal_xlsx()
        report = audit(fetch=lambda u: lookup[u])
        self.assertEqual(report["format_verified_count"], 5)
        self.assertIn("NO_PAIRED_CLINT_TARGET_CONFIRMED", report["research_status"])
        self.assertIn("mp1c00053_si_002.xlsx", report["sources"]["wegler_cloud_metadata"]["metadata_only"]["listed_target_supplements"])
        self.assertGreaterEqual(len(report["blocking_conditions"]), 4)
        self.assertNotIn("measured_clint", json.dumps(report))

    def test_pmc_2026_official_media_fallback_after_legacy_404(self):
        fixed = ("https://pmc-oa-opendata.s3.amazonaws.com/PMC8041379.1/")
        available = {fixed+"mp1c00053_si_001.pdf": b"%PDF-1.7\nexample",
                     fixed+"mp1c00053_si_002.xlsx": minimal_xlsx(),
                     fixed+"mp1c00053_si_003.xlsx": minimal_xlsx()}
        meta_url = SOURCES["wegler_cloud_metadata"][0]
        public_metadata = json.dumps({
            "pmcid": "PMC8041379",
            "doi": "10.1021/acs.molpharmaceut.1c00053",
            "version": 1, "license_code": "CC BY",
            "media_urls": list(available),
        }).encode()
        def fetch(source_url):
            if source_url == meta_url:
                return public_metadata
            if source_url in available:
                return available[source_url]
            raise OSError("old PMC URL not accessible")
        result = audit(fetch=fetch)
        self.assertEqual(result["format_verified_count"], 4)
        self.assertEqual(result["sources"]["wegler_s1_proteomics"]["state"], "FORMAT_VERIFIED")
        self.assertIn("2026 AWS", result["sources"]["wegler_s1_proteomics"]["source_distribution"])
        self.assertEqual(result["sources"]["yin_battail_metadata"]["state"], "UNAVAILABLE_OR_INVALID")

    def test_unavailable_sources_are_explicit_not_fake_pass(self):
        report = audit(fetch=lambda _: b"<html>bad gateway</html>")
        self.assertEqual(report["format_verified_count"], 0)
        self.assertTrue(all(
            item["state"] == "UNAVAILABLE_OR_INVALID"
            for item in report["sources"].values()
        ))


if __name__ == "__main__":
    unittest.main()
