"""Synthetic tests for Ohtsuki 2012 official supplemental-route audit."""
import unittest

from scripts.audit_ohtsuki2012_official_supplement import (
    CURRENT_ARTICLE_URL,
    LEGACY_SUPPLEMENT_URL,
    _candidate_links,
    audit,
    inspect_payload,
)


class OhtsukiOfficialSupplementAuditTests(unittest.TestCase):
    def test_html_manifest_discovers_only_approved_supplement_links(self):
        html="""
        <html><body>
        <a href="/content/suppl/2011/10/12/dmd.111.042259.DC1/supplement.pdf">supp PDF</a>
        <a href="https://ars.els-cdn.com/content/image/1-s2.0-S0090955624225861-mmc1.xlsx">data</a>
        <a href="https://example.com/private.xlsx">bad</a>
        <a href="/ordinary/article">ordinary</a>
        </body></html>
        """
        out=_candidate_links(LEGACY_SUPPLEMENT_URL,html)
        self.assertEqual(len(out),2)
        self.assertEqual({x["extension"] for x in out},{".pdf",".xlsx"})
        self.assertNotIn("example.com",repr(out))

    def test_payload_inspection_never_parses_donor_rows(self):
        html=(
            "<html><body>Supplemental material for 10.1124/dmd.111.042259 "
            '<a href="/lookup/suppl/doi:10.1124/dmd.111.042259/-/DC1/file.xlsx">'
            "download</a></body></html>"
        ).encode()
        out=inspect_payload(LEGACY_SUPPLEMENT_URL,"text/html; charset=utf-8",html)
        self.assertEqual(out["payload_kind"],"html")
        self.assertTrue(out["mentions_supplement"])
        self.assertTrue(out["mentions_doi"])
        self.assertFalse(out["donor_level_content_parsed"])
        self.assertEqual(out["candidate_supplement_links"][0]["extension"],".xlsx")

    def test_audit_abstains_when_official_routes_are_blocked(self):
        def fake(url):
            return {
                "state":"HTTP_ERROR",
                "http_status":403,
                "requested_url":{"host":"dmd.aspetjournals.org","path":"/blocked"},
                "donor_level_content_parsed":False,
            }
        out=audit(fetch=fake)
        self.assertEqual(out["status"],"ABSTAIN_OFFICIAL_SUPPLEMENT_FILE_NOT_RETRIEVED")
        self.assertEqual(out["candidate_link_count"],0)
        self.assertFalse(out["source_rows_acquired"])
        self.assertFalse(out["model_fitted"])

    def test_audit_reports_manifest_links_without_claiming_source_rows(self):
        def fake(url):
            if url==LEGACY_SUPPLEMENT_URL:
                return {
                    "state":"RETRIEVED",
                    "payload_kind":"html",
                    "candidate_supplement_links":[{
                        "scheme":"https","host":"dmd.aspetjournals.org",
                        "path":"/supp/file.xlsx","extension":".xlsx"
                    }],
                    "donor_level_content_parsed":False,
                }
            self.assertEqual(url,CURRENT_ARTICLE_URL)
            return {
                "state":"HTTP_ERROR","http_status":403,
                "donor_level_content_parsed":False,
            }
        out=audit(fetch=fake)
        self.assertEqual(out["status"],"OFFICIAL_SUPPLEMENT_LINKS_DISCOVERED")
        self.assertEqual(out["candidate_link_count"],1)
        self.assertFalse(out["source_rows_acquired"])
        self.assertFalse(out["rna_protein_activity_crosswalk_qualified"])


if __name__=="__main__":
    unittest.main()
