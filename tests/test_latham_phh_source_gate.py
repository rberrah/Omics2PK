"""Synthetic-only source-gate tests; no real person/PHH rows checked in."""
import unittest
from xml.etree import ElementTree as ET

from omics2pk.latham_phh_source_gate import (
    SOURCE_DOI, LathamPHHError, benchmark_latham_phh_table4,
)
from scripts.benchmark_latham_phh_original import (
    _get_table4_xml, PublishedTableUnavailable,
)


def fixture():
    out = []
    cyp2c8 = ["*1/*1"] * 8 + ["*1/*2"] * 2 + ["*1/*3"] * 3 + ["*1/*4"]
    for i in range(14):
        geno = "*1/*1" if i == 0 else "*1/*3" if i in (1, 2, 3) else (
            "*3/*6" if i == 4 else "*3/*7" if i == 5 else "*3/*3")
        category = "NM" if i == 0 else "IM" if i < 4 else "PM"
        protein = None if i == 6 else 2.0 + .2 * i
        prot3a5 = None if i == 6 else 1.0 + .1 * i
        out.append({
            "source_doi": SOURCE_DOI,
            "donor_id": "X" + chr(65 + i // 26) + chr(65 + i % 26),
            "cyp2c8_diplotype": cyp2c8[i],
            "cyp3a5_diplotype": geno,
            "published_cyp3a5_category": category,
            "cyp2c8_activity": 2.0 + i,
            "cyp3a_probe_activity": 4.0 + i * .4,
            "cyp2c8_protein": 1.0 + i * .2,
            "cyp3a4_protein": protein,
            "cyp3a5_protein": prot3a5,
            "cyp3a_total_protein": None if protein is None else protein + prot3a5,
            "imatinib_apparent_cl": 1.0 + i * .34,
            "n_desmethyl_imatinib_formation": .001 * (i + 1),
        })
    return out


def jats(rows, label="Table 4.", doi=SOURCE_DOI):
    root = ET.Element("article")
    front = ET.SubElement(root, "article-meta")
    ET.SubElement(front, "article-id", {"pub-id-type": "doi"}).text = doi
    wrap = ET.SubElement(root, "table-wrap")
    ET.SubElement(wrap, "label").text = label
    table = ET.SubElement(wrap, "table")
    body = ET.SubElement(table, "tbody")
    def c(value):
        return "−" if value is None else str(value)
    for row in rows:
        values = [
            row["donor_id"], "42", "F", "W", "27.5", c(row["cyp2c8_activity"]),
            c(row["cyp3a_probe_activity"]), "4.0", row["cyp2c8_diplotype"],
            "*1/*1", row["cyp3a5_diplotype"], row["published_cyp3a5_category"],
            c(row["cyp2c8_protein"]), c(row["cyp3a4_protein"]),
            c(row["cyp3a5_protein"]), c(row["cyp3a_total_protein"]),
            c(row["imatinib_apparent_cl"]), c(row["n_desmethyl_imatinib_formation"]),
        ]
        tr = ET.SubElement(body, "tr")
        for value in values:
            ET.SubElement(tr, "td").text = value
    return ET.tostring(root, encoding="utf-8", xml_declaration=True)


class LathamSourceGateTests(unittest.TestCase):
    def test_four_arms_same_thirteen_donors_without_source_rows(self):
        r = benchmark_latham_phh_table4(fixture())
        self.assertEqual(r["source_rows"], 14)
        self.assertEqual(r["n_same_protein_complete_donors_all_arms"], 13)
        self.assertEqual(r["cyp3a5_PM_phenotype_protein_complete"], 9)
        self.assertEqual(r["cyp3a5_PM_with_positive_targeted_peptide_quantitation"], 9)
        self.assertEqual(r["cyp3a5_PM_with_positive_shared_CYP3A_probe_activity"], 10)
        self.assertTrue(r["positive_peptide_or_shared_probe_is_not_CYP3A5_specific_function"])
        self.assertEqual(r["cyp3a5_phenotype_group_counts_original_14"],
                         {"NM": 1, "IM": 3, "PM": 10})
        self.assertEqual(len(r["heldout_MAE_by_predeclared_arm_source_table_scale"]), 4)
        self.assertEqual(set(r["genotype_stratified_MAE_exploratory"]),
                         {"CYP3A5_nonexpresser_PM","CYP3A5_IM_or_NM",
                          "CYP2C8_reference_star1_star1","CYP2C8_variant_diplotype"})
        self.assertEqual(r["genotype_stratified_MAE_exploratory"]
                         ["CYP3A5_IM_or_NM"]["status"],
                         "ABSTAIN_SMALL_ORIGINAL_GENOTYPE_STRATUM")
        self.assertEqual(r["genotype_stratified_MAE_exploratory"]
                         ["CYP3A5_nonexpresser_PM"]["n_original_donors"],9)
        self.assertTrue(r["no_RNA_or_endogenous_metabolomics_in_original_PHH_table4"])
        self.assertTrue(r["n_desmethyl_formation_is_outcome_not_predose_predictor"])
        self.assertFalse(r["in_vivo_PopPK_or_tacrolimus_prediction_demonstrated"])
        self.assertNotIn("XAA", str(r))

    def test_metabolite_outcome_changes_never_enter_predictor_vectors(self):
        a = fixture()
        original = benchmark_latham_phh_table4(a)
        for row in a:
            row["n_desmethyl_imatinib_formation"] *= 1000
        again = benchmark_latham_phh_table4(a)
        self.assertEqual(original["heldout_MAE_by_predeclared_arm_source_table_scale"],
                         again["heldout_MAE_by_predeclared_arm_source_table_scale"])

    def test_duplicate_or_missing_donor_fails_closed(self):
        a = fixture()
        a[2]["donor_id"] = a[1]["donor_id"]
        with self.assertRaisesRegex(LathamPHHError, "duplicate"):
            benchmark_latham_phh_table4(a)
        with self.assertRaisesRegex(LathamPHHError, "exactly 14"):
            benchmark_latham_phh_table4(fixture()[:-1])

    def test_wrong_phenotype_cannot_override_genotype(self):
        a = fixture()
        a[4]["published_cyp3a5_category"] = "NM"
        with self.assertRaisesRegex(LathamPHHError, "classification"):
            benchmark_latham_phh_table4(a)

    def test_missing_protein_not_imputed(self):
        a = fixture()
        a[6]["cyp3a4_protein"] = 2.0
        with self.assertRaisesRegex(LathamPHHError, "Incomplete correlated"):
            benchmark_latham_phh_table4(a)

    def test_cyp3a_total_consistency(self):
        a = fixture()
        a[2]["cyp3a_total_protein"] *= 5
        with self.assertRaisesRegex(LathamPHHError, "does not reproduce total"):
            benchmark_latham_phh_table4(a)

    def test_xml_reader_round_trip_and_no_row_order_dependence(self):
        a = fixture()
        parsed = _get_table4_xml(jats(a))
        out = benchmark_latham_phh_table4(parsed)
        reversed_out = benchmark_latham_phh_table4(list(reversed(parsed)))
        self.assertEqual(out["heldout_MAE_by_predeclared_arm_source_table_scale"],
                         reversed_out["heldout_MAE_by_predeclared_arm_source_table_scale"])
        self.assertEqual(out["n_same_protein_complete_donors_all_arms"], 13)

    def test_xml_wrong_doi_and_wrong_label_fail_closed(self):
        with self.assertRaisesRegex(PublishedTableUnavailable, "DOI"):
            _get_table4_xml(jats(fixture(), doi="10.0000/not-the-source"))
        with self.assertRaisesRegex(PublishedTableUnavailable, "Table 4"):
            _get_table4_xml(jats(fixture(), label="Table 1."))

    def test_xml_row_layout_change_fails_closed(self):
        doc = ET.fromstring(jats(fixture()))
        table = doc.find(".//tbody")
        table[0].remove(table[0][-1])
        with self.assertRaisesRegex(PublishedTableUnavailable, "layout changed"):
            _get_table4_xml(ET.tostring(doc))

    def test_html_parser_if_beautifulsoup_is_available(self):
        try:
            from bs4 import BeautifulSoup  # noqa: F401
        except ImportError:
            self.skipTest("HTML fallback dependency only in external-source CI")
        from scripts.benchmark_latham_phh_original import _get_table4_html
        from html import escape
        rows = fixture()
        from xml.etree import ElementTree as ET
        root = ET.fromstring(jats(rows))
        headers = "".join(f"<th>{s}</th>" for s in
                          ("Donor ID", "CYP2C8 Genotype", "CYP3A5 Genotype", "Imatinib Clearance"))
        html_rows = []
        for tr in root.findall(".//tbody/tr"):
            html_rows.append("<tr>" + "".join(
                "<td>" + escape(td.text or "") + "</td>" for td in tr.findall("td")) + "</tr>")
        html = ("<html><body>" + SOURCE_DOI + "<table><thead><tr>" + headers +
                "</tr></thead><tbody>" + "".join(html_rows) + "</tbody></table></body></html>")
        output = _get_table4_html(html.encode())
        self.assertEqual(len(output), 14)


if __name__ == "__main__":
    unittest.main()
