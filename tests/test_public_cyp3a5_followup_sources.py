"""Synthetic structural fixtures; no original donor data in tests."""
import unittest

from scripts.audit_public_cyp3a5_followup_sources import (
    PublicSourceAuditError, _xml_metadata, _html_metadata,
    SOURCES, audit_primary_source,
)

DOI = "10.1124/dmd.113.051094"


def xml_source(doi=DOI):
    return f"""<article><front><article-meta><article-id pub-id-type="doi">{doi}</article-id></article-meta></front><body>
<table-wrap><label>TABLE 2.</label><table><tbody>
<tr><td>D1</td><td>*3/*3</td><td>0.8</td></tr>
<tr><td>D2</td><td>*1/*6</td><td>22.1</td></tr>
</tbody></table></table-wrap>
</body><back><supplementary-material mimetype="application" mime-subtype="pdf">
<label>Supplement</label></supplementary-material></back></article>""".encode()


class SourceAuditTests(unittest.TestCase):
    def test_doi_locked_source_table_counts_no_individual_rows(self):
        r = _xml_metadata(xml_source(), DOI)
        self.assertEqual(r["tables"][0]["rows_total"], 2)
        self.assertEqual(r["tables"][0]["data_cell_column_count_distribution"], [(3, 2)])
        self.assertEqual(r["n_original_donor_numeric_matrix_verified"], 0)
        self.assertNotIn("D1", str(r))
        self.assertNotIn("*1/*6", str(r))

    def test_wrong_doi_fails(self):
        with self.assertRaisesRegex(PublicSourceAuditError,"DOI mismatch"):
            _xml_metadata(xml_source(doi="10.99/wrong"),DOI)

    def test_supplements_are_inventory_only(self):
        r = _xml_metadata(xml_source(), DOI)
        self.assertEqual(r["supplementary_materials"][0]["extension"], "application/pdf")
        self.assertTrue(r["no_public_patient_values_emitted"])

    def test_unsupported_source_label_fails_closed(self):
        with self.assertRaisesRegex(PublicSourceAuditError,"allowlisted"):
            audit_primary_source("other")

    def test_only_three_pinned_original_sources(self):
        self.assertEqual({k:v[0] for k,v in SOURCES.items()}, {
            "grangeon":"PMC12735975", "bissada":"PMC6815940",
            "topletz":"PMC3876810"})

    def test_html_fallback_metadata_synthetic_if_parser_available(self):
        try:
            from bs4 import BeautifulSoup  # noqa:F401
        except ImportError:
            self.skipTest("beautifulsoup4 available only in optional external-source job")
        html = ('<html><body>'+DOI+
                '<table><caption>TABLE 2</caption><tr><th>HLM</th></tr>'
                '<tr><td>D1</td><td>12</td><td>3</td></tr></table>'
                '<a href="article-s001.pdf">Supplement S1</a></body></html>').encode()
        r=_html_metadata(html, DOI)
        self.assertEqual(r["tables"][0]["n_rows_with_3plus_data_cells"],1)
        self.assertEqual(r["n_original_donor_numeric_matrix_verified"],0)
        self.assertNotIn("D1",str(r))
        self.assertEqual(r["supplementary_materials"][0]["extension"],"pdf")


if __name__=="__main__":
    unittest.main()
