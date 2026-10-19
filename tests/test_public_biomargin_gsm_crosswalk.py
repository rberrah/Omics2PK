"""No real GEO human identifiers in CI: fictional sample-join contracts only."""
import gzip
import unittest
from scripts.qualify_public_biomargin_geo_gsm_bios_crosswalk import (
    inspect_series_matrix,peek_counts_columns,PublicBIOMARGINAuditError)


def gz(lines):
    return gzip.compress("\n".join(lines).encode())


class GSESourceJoinFictionalTests(unittest.TestCase):
    def test_exact_GSM_and_BIOS_crosswalk_exists(self):
        authors={"BIOS.03.0001","BIOS.06.0002"}
        matrix=gz([
            '!Sample_geo_accession\t"GSM_FAKE1"\t"GSM_FAKE2"',
            '!Sample_title\t"specimen1"\t"specimen2"',
            '!Sample_characteristics_ch1\t"participant: BIOS.03.0001"\t"participant: BIOS.06.0002"',
            '!series_matrix_table_begin',
        ])
        report=inspect_series_matrix(matrix,authors,{"GSM_FAKE1","GSM_FAKE2"},
                                     production=False)
        self.assertEqual(report["n_exact_three_way_GSM_count_and_author_BIOS_matches"],2)
        self.assertTrue(report["identifiable_external_join_without_fuzzy_or_row_order"])
        self.assertNotIn("BIOS.",str(report))
        self.assertNotIn("GSM_FAKE",str(report))

    def test_GEO_titles_without_author_crosswalk_do_not_count(self):
        matrix=gz([
            '!Sample_geo_accession\t"GSM_FAKE1"\t"GSM_FAKE2"',
            '!Sample_title\t"specimen1"\t"specimen2"',
            '!Sample_characteristics_ch1\t"tissue: blood"\t"tissue: blood"',
            '!series_matrix_table_begin'])
        report=inspect_series_matrix(matrix,{"BIOS.03.0001"},
                                      {"GSM_FAKE1","GSM_FAKE2"},production=False)
        self.assertFalse(report["identifiable_external_join_without_fuzzy_or_row_order"])
        self.assertEqual(report["n_exact_three_way_GSM_count_and_author_BIOS_matches"],0)

    def test_one_GEO_sample_multiple_biological_owners_rejected(self):
        matrix=gz([
            '!Sample_geo_accession\t"GSM_FAKE1"\t"GSM_FAKE2"',
            '!Sample_characteristics_ch1\t"BIOS.03.0001 and BIOS.06.0002"\t"none"',
            '!series_matrix_table_begin'])
        with self.assertRaisesRegex(PublicBIOMARGINAuditError,"multiple original people"):
            inspect_series_matrix(matrix,set(),{"GSM_FAKE1"},production=False)


if __name__=="__main__":
    unittest.main()
