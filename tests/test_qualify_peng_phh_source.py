"""Unit tests use only six fake donor metadata and fake DOCX/GEO archives."""
from io import BytesIO
import gzip
import json
import tarfile
import unittest
import zipfile

from scripts.qualify_peng_phh_source import (
    parse_aws_metadata, geo_sample_metadata, geo_tar_inventory, word_table_counts,
    EXPECTED_BATCHES,
)

def soft_fixture():
    sections = ["^SERIES = GSE289636", "!Series_supplementary_file = ftp://GSE289636_RAW.tar"]
    for i,b in enumerate(EXPECTED_BATCHES):
        sections += [f"^SAMPLE = GSM879564{i+2}", f"!Sample_title = primary human hepatocyte, {b}, rep1"]
    return gzip.compress(("\n".join(sections)+"\n").encode())

class SourceTests(unittest.TestCase):
    def test_six_original_geo_donor_groups(self):
        result=geo_sample_metadata(soft_fixture())
        self.assertEqual(result["n_geo_samples"],6)
        self.assertEqual(result["donor_batches_verified"],6)
        self.assertTrue(result["has_raw_supplement"])

    def test_fail_closed_wrong_donor_batch(self):
        altered=gzip.decompress(soft_fixture()).replace(b"416, rep1",b"999, rep1")
        with self.assertRaisesRegex(ValueError,"recover original GEO donor"):
            geo_sample_metadata(gzip.compress(altered))

    def test_approved_pmc_document_only(self):
        source="s3://pmc-oa-opendata/PMC12532551.1/13619_2025_258_MOESM1_ESM.docx?md5=abc"
        meta={"pmcid":"PMC12532551","version":"1","media_urls":[source]}
        url,license_code = parse_aws_metadata(json.dumps(meta).encode())
        self.assertTrue(url.startswith("https://pmc-oa-opendata.s3.amazonaws.com/"))
        self.assertEqual(license_code,"None")
        meta["pmcid"]="WRONG"
        with self.assertRaisesRegex(ValueError,"Unexpected article"):
            parse_aws_metadata(json.dumps(meta).encode())

    def test_docx_table_counts_no_individual_numerical_values(self):
        doc=(
          '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
          '<w:body><w:tbl><w:tr><w:tc><w:p><w:r><w:t>CYP3A4</w:t></w:r></w:p></w:tc>'
          '<w:tc><w:p><w:r><w:t>SCR</w:t></w:r></w:p></w:tc></w:tr>'
          '<w:tr><w:tc><w:p><w:r><w:t>PHH025</w:t></w:r></w:p></w:tc>'
          '<w:tc><w:p><w:r><w:t>123.4</w:t></w:r></w:p></w:tc></w:tr>'
          '</w:tbl></w:body></w:document>')
        z=BytesIO()
        with zipfile.ZipFile(z,"w") as f:
            f.writestr("word/document.xml",doc)
        audit=word_table_counts(z.getvalue())
        self.assertEqual(audit["table_count"],1)
        self.assertTrue(audit["tables"][0]["contains_cyp_marker"])
        self.assertTrue(audit["tables"][0]["contains_metabolic_terms"])
        self.assertNotIn("123.4", repr(audit))

    def test_raw_tar_inventory_and_checksums(self):
        z=BytesIO()
        with tarfile.open(fileobj=z,mode="w") as t:
            for b in EXPECTED_BATCHES:
                info=tarfile.TarInfo("GSM_synthetic_"+b+".csv")
                blob=b"CYP3A4,1\n"
                info.size=len(blob)
                t.addfile(info,BytesIO(blob))
        report=geo_tar_inventory(z.getvalue())
        self.assertEqual(report["n_files"],6)
        self.assertTrue(all(v==1 for v in report["batch_filename_presence"].values()))
        self.assertNotIn("CYP3A4,1",repr(report))

if __name__=="__main__":
    unittest.main()
