"""Synthetic tests for Zhang 2016 official supplemental-route audit."""
import unittest

from scripts.audit_zhang2016_official_supplement import (
    CURRENT_ARTICLE_URL,
    LEGACY_SUPPLEMENT_URL,
    _candidate_links,
    audit,
    inspect_payload,
)


class ZhangOfficialSupplementAuditTests(unittest.TestCase):
    def test_discovers_only_approved_supplement_links(self):
        html="""
        <html><body>
        <a href="/content/suppl/2016/05/09/jpet.116.233635.DC1/table.xlsx">table</a>
        <a href="https://ars.els-cdn.com/content/image/1-s2.0-S0022356524193494-mmc1.zip">archive</a>
        <a href="https://example.com/untrusted.xlsx">bad</a>
        </body></html>
        """
        out=_candidate_links(LEGACY_SUPPLEMENT_URL,html)
        self.assertEqual(len(out),2)
        self.assertEqual({x["extension"] for x in out},{".xlsx",".zip"})
        self.assertNotIn("example.com",repr(out))

    def test_html_payload_is_manifest_only_not_source_rows(self):
        html=(
            "<html><body>Supplemental material 10.1124/jpet.116.233635 "
            '<a href="/lookup/suppl/doi:10.1124/jpet.116.233635/-/DC1/data.xlsx">'
            "data</a></body></html>"
        ).encode()
        out=inspect_payload(LEGACY_SUPPLEMENT_URL,"text/html",html)
        self.assertEqual(out["payload_kind"],"html")
        self.assertTrue(out["mentions_supplement"])
        self.assertTrue(out["mentions_doi"])
        self.assertFalse(out["donor_level_content_parsed"])

    def test_blocked_routes_abstain(self):
        def fake(url):
            return {
                "state":"HTTP_ERROR","http_status":403,
                "requested_url":{"host":"jpet.aspetjournals.org","path":"/blocked"},
                "donor_level_content_parsed":False,
            }
        out=audit(fetch=fake)
        self.assertEqual(out["status"],"ABSTAIN_OFFICIAL_SUPPLEMENT_FILE_NOT_RETRIEVED")
        self.assertEqual(out["published_n_human_livers"],100)
        self.assertTrue(out["article_reports_same_cohort_cyp_mrna_protein_activity"])
        self.assertFalse(out["source_rows_acquired"])
        self.assertFalse(out["model_fitted"])

    def test_manifest_link_is_not_crosswalk_qualification(self):
        def fake(url):
            if url==LEGACY_SUPPLEMENT_URL:
                return {
                    "state":"RETRIEVED","payload_kind":"html",
                    "candidate_supplement_links":[{
                        "scheme":"https","host":"jpet.aspetjournals.org",
                        "path":"/supp/data.xlsx","extension":".xlsx"
                    }],
                    "donor_level_content_parsed":False,
                }
            self.assertEqual(url,CURRENT_ARTICLE_URL)
            return {"state":"HTTP_ERROR","http_status":403,
                    "donor_level_content_parsed":False}
        out=audit(fetch=fake)
        self.assertEqual(out["status"],"OFFICIAL_SUPPLEMENT_LINKS_DISCOVERED")
        self.assertEqual(out["candidate_link_count"],1)
        self.assertFalse(out["rna_protein_activity_crosswalk_qualified"])


if __name__=="__main__":
    unittest.main()
