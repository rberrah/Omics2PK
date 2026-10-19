"""Pure fictional BIOMARGIN source contract tests; no real IDs or patient data."""
from gzip import compress
from hashlib import sha1
import json
from unittest.mock import patch
import unittest

from scripts.qualify_public_biomargin_131_original_source import (
    PublicBIOMARGINAuditError,
    audit_author_notebook, audit_public_geo_counts_gzip, ORIGINAL_MOFA_VIEWS,
)


def synthetic_notebook(n=131, duplicate=False):
    centers=("03","06","09","13")
    names=[f"BIOS.{centers[i%4]}.{i+1:04d}" for i in range(n)]
    if duplicate:
        names[-1]=names[0]
    source='samples.order = c('+','.join('"'+s+'"' for s in names)+')\n'
    source+='\n'.join('X%d <- read.csv2(file="%s")' % (i, f)
                      for i,f in enumerate(ORIGINAL_MOFA_VIEWS))
    fake={"cells":[{"cell_type":"code","source":[source]}]}
    payload=json.dumps(fake).encode()
    blob=sha1(b"blob "+str(len(payload)).encode()+b"\0"+payload).hexdigest()
    return payload,blob


class BIOMARGINSourceSyntheticTests(unittest.TestCase):
    def test_131_original_author_identifier_joins_and_six_views(self):
        payload,blob=synthetic_notebook()
        with patch("scripts.qualify_public_biomargin_131_original_source.AUTHOR_NOTEBOOK_BLOB",blob):
            subjects, report=audit_author_notebook(payload)
        self.assertEqual(len(subjects),131)
        self.assertEqual(report["n_original_centers"],4)
        self.assertEqual(report["n_source_molecular_views"],6)
        self.assertEqual(report["n_complete_six_view_subjects_published_not_reverified_in_this_audit"],35)
        self.assertTrue(report["no_drug_PK_linked_or_evaluated"])
        self.assertNotIn("BIOS.",json.dumps(report))

    def test_duplicate_original_person_abstains(self):
        payload,blob=synthetic_notebook(duplicate=True)
        with patch("scripts.qualify_public_biomargin_131_original_source.AUTHOR_NOTEBOOK_BLOB",blob):
            with self.assertRaisesRegex(PublicBIOMARGINAuditError,"131 independent"):
                audit_author_notebook(payload)

    def test_modified_original_source_hash_fails_closed(self):
        payload,_=synthetic_notebook()
        with self.assertRaisesRegex(PublicBIOMARGINAuditError,"blob identity"):
            audit_author_notebook(payload)

    def test_genuine_no_row_order_assumptions(self):
        payload,blob=synthetic_notebook()
        with patch("scripts.qualify_public_biomargin_131_original_source.AUTHOR_NOTEBOOK_BLOB",blob):
            people,_=audit_author_notebook(payload)
        # Deliberately GEO-style IDs not equal to authors' original biological IDs.
        gz=compress(("gene,GSM_fake_1,GSM_fake_2,GSM_fake_3\n"
                     "FAKE_GENE_1,1,2,3\nFAKE_GENE_2,0,0,0\n").encode())
        qc=audit_public_geo_counts_gzip(gz,people,production=False)
        self.assertEqual(qc["n_direct_exact_BIOS_column_matches"],0)
        self.assertTrue(qc["join_status"].startswith("CROSSWALK_REQUIRED"))
        self.assertEqual(qc["number_gene_feature_rows"],2)
        self.assertEqual(qc["number_assay_columns"],3)
        self.assertNotIn("GSM_fake_",json.dumps(qc))

    def test_fake_negative_count_blocked(self):
        _, blob=synthetic_notebook()
        with self.assertRaisesRegex(PublicBIOMARGINAuditError,"negative"):
            audit_public_geo_counts_gzip(
                compress(b"gene,sample1,sample2\nG1,-1,1\n"),set(),production=False)


if __name__=="__main__":
    unittest.main()
