"""Fictional-only GNPS spectral-match and bile proxy tests."""
import unittest

try:
    import numpy as np
    from scripts.validate_private_jarmusch_bile_ms1_proxy import (
        BileProxyError, candidate_scans_from_spectral_matches,
        compare_proxy)
except ModuleNotFoundError as exc:
    if exc.name not in ("numpy","sklearn"):
        raise
    np = None


@unittest.skipIf(np is None,"Optional NumPy benchmark dependencies required")
class OriginalSpectralProxyTests(unittest.TestCase):
    def test_no_unrelated_metabolites_or_drug_peaks_called_bile_acids(self):
        rows=[
            {"Compound_Name":"Cholic acid","#Scan#":"2"},
            {"Compound_Name":"Spectral Match to Cholic acid from NIST14","#Scan#":"3"},
            {"Compound_Name":"Chenodeoxycholic acid","#Scan#":"4"},
            {"Compound_Name":"Midazolam","#Scan#":"5"},
            {"Compound_Name":"5-hydroxyomeprazole","#Scan#":"6"},
        ]
        self.assertEqual(candidate_scans_from_spectral_matches(rows),[
            "GNPS_MS1:2","GNPS_MS1:3","GNPS_MS1:4"])

    def test_false_scan_identity_abstains(self):
        with self.assertRaisesRegex(BileProxyError,"scan index"):
            candidate_scans_from_spectral_matches([{
                "Compound_Name":"Cholic acid","#Scan#":"private_bad_scan"}])

    def test_no_overinterpretation_synthetic_proxy_and_fold_groups(self):
        rng=np.random.default_rng(20261019)
        x=np.exp(rng.normal(0,1,(13,7)))
        clinical=np.column_stack((rng.uniform(50,100,13),rng.integers(0,2,13)))
        genotype=np.asarray([["X/X","Y/Y"] if i%2 else ["X/Z","Y/Z"]
                             for i in range(13)],dtype=object)
        pk=np.exp(rng.normal(2,0.2,(13,2)))
        pk[-1,1]=float("nan")
        genes=["GNPS_MS1:"+str(i) for i in range(1,8)]
        report=compare_proxy(x,clinical,genotype,pk,genes,genes[:3])
        self.assertEqual(report["day1"]["n_independent_original_people"],13)
        self.assertEqual(report["day9"]["n_independent_original_people"],12)
        self.assertFalse(report["day1"]["independent_external_validation"])
        self.assertTrue(np.isfinite(report["day1"]["delta_MAE_incremental"]))


if __name__=="__main__":
    unittest.main()
