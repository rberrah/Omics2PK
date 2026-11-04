"""Synthetic-only tests for public MSP BinDiscover source-ID crosswalk gates."""
import tempfile
from pathlib import Path
import unittest
from scripts.audit_rotroff_vs_bindiscover_public_msp import scan_msp,TOP
class OriginalMSPSourceContractTests(unittest.TestCase):
    def test_source_bin_is_metadata_not_mz(self):
        sample=('NAME: unknown bin 222169\n'
                'DB#: 222169\n'
                'RETENTIONINDEX: 230000\n'
                'NUM PEAKS: 6\n'
                '73 100\n104 90\n147 70\n182 45\n200 15\n218 14\n\n')
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'synthetic.msp'
            path.write_text(sample*1001,encoding='utf8')
            output=scan_msp(path,set(TOP),path.stat().st_size)
        self.assertEqual(output['n_reference_MSP_records'],1001)
        self.assertEqual(output['n_exact_identity_field_bin_matches'],1)
        self.assertEqual(output['matched_source_bins']['222169']['n_candidate_records'],1001)
        self.assertFalse(output['full_spectra_ion_pairs_exported'])
        self.assertNotIn('73 100',str(output))
    def test_number_in_retention_index_cannot_reidentify_bin(self):
        sample=('NAME: biological unknown\nRETENTIONINDEX: 222169\n'
                'NUM PEAKS: 6\n73 100\n85 50\n90 20\n100 10\n110 8\n120 9\n\n')
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'synthetic.msp'
            path.write_text(sample*1001,encoding='utf8')
            out=scan_msp(path,set(TOP),path.stat().st_size)
        self.assertEqual(out['n_exact_identity_field_bin_matches'],0)
    def test_source_size_contract_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)/'nonsense.msp'
            p.write_text('fake')
            with self.assertRaises(ValueError):
                scan_msp(p,set(TOP),40000000)
if __name__=='__main__':unittest.main()
