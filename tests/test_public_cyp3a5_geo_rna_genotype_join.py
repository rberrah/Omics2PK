"""GEO original paired-source metadata parser; synthetic records only."""
import csv
import gzip
import io
import unittest
from scripts.audit_public_cyp3a5_geo_rna_genotype_join import (
    parse_geo_series_matrix,GEOQualificationError,
)
def fixture(n=60,prefix="MP-"):
    out=io.StringIO()
    w=csv.writer(out,delimiter="\t")
    w.writerow(["!Sample_geo_accession"]+[f"GSM{3500000+i}" for i in range(n)])
    w.writerow(["!Sample_title"]+[f"{prefix}{i}" for i in range(n)])
    return gzip.compress(out.getvalue().encode())
class GEOSourcePairingTest(unittest.TestCase):
    def test_valid_source_metadata_only_is_not_functional_genotype(self):
        a=parse_geo_series_matrix(fixture(),series="RNA_seq_GSE124074")
        b=parse_geo_series_matrix(fixture(n=96),series="genotyping_GSE147628")
        self.assertEqual(a["n_GEO_sample_records"],60)
        self.assertEqual(b["n_GEO_sample_records"],96)
        self.assertEqual(len(set(a["original_titles"])&set(b["original_titles"])),60)

    def test_mismatch_no_false_same_donor_claim(self):
        a=parse_geo_series_matrix(fixture(),series="RNA")
        b=parse_geo_series_matrix(fixture(prefix="other-"),series="DNA")
        self.assertEqual(len(set(a["original_titles"])&set(b["original_titles"])),0)

    def test_failed_source_abstains(self):
        with self.assertRaises(GEOQualificationError):
            parse_geo_series_matrix(b"bad_gzip"*200,series="RNA")
        blob=gzip.compress(b"!Sample_title\t\"MP0\"\n")
        with self.assertRaisesRegex(GEOQualificationError,"pairing"):
            parse_geo_series_matrix(blob,series="RNA")
        repeated=io.StringIO()
        w=csv.writer(repeated,delimiter="\t")
        w.writerow(["!Sample_geo_accession","GSM1","GSM1"])
        w.writerow(["!Sample_title","A","B"])
        with self.assertRaisesRegex(GEOQualificationError,"duplicated"):
            parse_geo_series_matrix(gzip.compress(repeated.getvalue().encode()),series="RNA")
if __name__=="__main__":
    unittest.main()
