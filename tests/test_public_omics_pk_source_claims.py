"""Synthetic-only 93-source contract; no SRA download during unit tests."""
import io
import unittest
from unittest.mock import patch

from scripts import audit_degraeve_93_kidney_tac_public_sources as source


class MockResponse:
    def __init__(self, text):self.data=text.encode()
    def __enter__(self):return self
    def __exit__(self,*args):return False
    def read(self,n):return self.data[:n]


def fake_index(n):
    rows=["run_accession\tsample_accession\tsample_alias\tstudy_accession"]
    rows.extend(f"SYNTH_RUN{i}\tSYNTH_SAMPLE{i}\tSYNTH_ALIAS{i}\tSYNTH_PROJECT"
                for i in range(n))
    return "\n".join(rows)+"\n"


class PublicOmicsPKSourceTests(unittest.TestCase):
    def test_93_original_biological_samples_not_93_PK_profiles(self):
        with patch.object(source,"urlopen",return_value=MockResponse(fake_index(93))):
            report=source.ena_index()
        self.assertEqual(report["n_original_distinct_sequence_samples"],93)
        self.assertEqual(report["n_original_16s_runs"],93)
        self.assertNotIn("SYNTH_ALIAS",repr(report))

    def test_source_count_is_not_silently_imputed(self):
        with patch.object(source,"urlopen",return_value=MockResponse(fake_index(92))):
            with self.assertRaisesRegex(ValueError,"Unexpected original 93-person"):
                source.ena_index()

    def test_AUC_estimate_cannot_masquerade_as_observed_NCA(self):
        self.assertEqual(source.EXPECTED_HUMAN_SUBJECTS,93)
        self.assertEqual(source.PROJECT,"PRJNA929904")
        self.assertEqual(source.DOI,"10.1002/cpt.3077")


if __name__=="__main__":
    unittest.main()
