"""Synthetic JATS source identity, donor-row and publication checks; no originals."""
from __future__ import annotations
import unittest
import xml.etree.ElementTree as ET

from omics2pk.smx_nat_original_jats import (
    rows_from_original_jats, published_source_validation, SOURCE_DOI,
)
from omics2pk.smx_nat_local_validation import NATSourceError, evaluate


HEADERS=[
    "Donor","Gender","Race","Age (y)","BMI (kg/m2)",
    "NAT1 Genotype","NAT2 Genotype","NAT2 Genotype-Inferred Phenotype",
    "NA-SMX (pmol/min per mill cells)",
    "NA-SMZ (pmol/min per mill cells)",
    "NA-PABA (pmol/min per mill cells)",
    "NAT2 Protein Concentration (pmol/mg) (NAT2.EFLNSHLLPK)",
    "NAT2 Protein Concentration (pmol/mg) (NAT2.DNTDLVEFK)",
]


def synthetic_source(n=26):
    article=ET.Element("article")
    front=ET.SubElement(article,"front")
    metadata=ET.SubElement(front,"article-meta")
    ET.SubElement(metadata,"article-id",{"pub-id-type":"doi"}).text=SOURCE_DOI
    ET.SubElement(metadata,"article-id",{"pub-id-type":"pmc"}).text="12861138"
    body=ET.SubElement(article,"body")
    wrap=ET.SubElement(body,"table-wrap")
    table=ET.SubElement(wrap,"table")
    header=ET.SubElement(ET.SubElement(table,"thead"),"tr")
    for h in HEADERS:
        ET.SubElement(header,"th").text=h
    phenotypes=["Slow"]*9 + ["Intermediate"]*6 + ["Rapid"]*5 + ["-"]*6
    for i in range(n):
        tr=ET.SubElement(ET.SubElement(table,"tbody"),"tr")
        fields=[
            "X%02d"%i,"Female" if i%2 else "Male","WH",str(19+i),"24.5",
            "*4/*4","*4/*4",phenotypes[i],
            "%.3f"%(0.12+i*.006),
            "%.3f"%(7+i*.5) if i<21 else "-",
            "%.3f"%(15+i*.2) if i<21 else "-",
            "%.3f"%(0.5+i*.03),
            "%.3f"%(0.7+i*.015),
        ]
        for val in fields:
            ET.SubElement(tr,"td").text=val
    return ET.tostring(article,encoding="utf-8")


class OriginalSourceJATSTests(unittest.TestCase):
    def test_source_identity_schema_and_cardinality(self):
        rows=rows_from_original_jats(synthetic_source())
        self.assertEqual(len(rows),26)
        result=evaluate(rows)
        self.assertEqual(result["n_assayed_with_both_NAT_probes"],21)
        self.assertEqual(result["n_genotyped"],20)
        self.assertTrue(result["original_rows_redistributed"] is False)

    def test_wrong_doi_abstain(self):
        xml=synthetic_source().replace(SOURCE_DOI.encode(),b"10.1016/bogus")
        with self.assertRaisesRegex(NATSourceError,"DOI/PMC"):
            rows_from_original_jats(xml)

    def test_wrong_column_protein_peptide_abstain(self):
        xml=synthetic_source().replace(b"NAT2.EFLNSHLLPK",b"NAT2.UNKNOWN")
        with self.assertRaisesRegex(NATSourceError,"peptide"):
            rows_from_original_jats(xml)

    def test_25_donor_rows_abstain(self):
        with self.assertRaisesRegex(NATSourceError,"26 donor rows"):
            rows_from_original_jats(synthetic_source(25))

    def test_published_correlations_must_not_be_automatic_from_shape(self):
        rows=rows_from_original_jats(synthetic_source())
        with self.assertRaisesRegex(NATSourceError,"Spearman"):
            published_source_validation(rows)

    def test_html_or_doctype_is_not_original_xml(self):
        with self.assertRaises(NATSourceError):
            rows_from_original_jats(b"<!DOCTYPE foo><article></article>")
        with self.assertRaises(NATSourceError):
            rows_from_original_jats(b"<html>RECAPTCHA</html>")


if __name__=="__main__":
    unittest.main()
