"""Synthetic-only identity/analytical guards for original NCI60 threeomic benchmark."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
try:
    import numpy as np
    from benchmark_original_nci60_threeomics_swath import canon,distinct_map,run,SourceNotQualified
    _BENCHMARK_AVAILABLE=True
except ImportError:
    np=None
    _BENCHMARK_AVAILABLE=False
@unittest.skipUnless(_BENCHMARK_AVAILABLE,"NCI-60 benchmark NumPy/scikit-learn optional extra not installed in offline-only CI")
class TestOriginalThreeOmicsGuards(unittest.TestCase):
    def test_case_punctuation_identity(self):
        self.assertEqual(canon("MDA-MB-231"),"MDAMB231")
        self.assertEqual(canon("CCRF_CEM"),"CCRFCEM")
        self.assertRaises(SourceNotQualified,distinct_map,["A-B","AB"],"different aliases")
    def test_refuse_missing_origin_and_mismatched_lineage(self):
        rng=np.random.default_rng(22)
        ids=tuple("line_"+str(i) for i in range(47))
        cls=tuple("BREAST" for _ in ids)
        y=rng.normal(size=47)
        r=rng.normal(size=(47,310))
        m=rng.normal(size=(47,125))
        p=rng.normal(size=(60,315))
        all_ids=ids+tuple("extra_"+str(i) for i in range(13))
        self.assertRaises(SourceNotQualified,run,ids,cls,y,r,m,all_ids,("NSCLC",)+("BREAST",)*59,p)
        self.assertRaises(SourceNotQualified,run,ids,cls,y,r,m,tuple("ineligible_"+str(i) for i in range(60)),("BREAST",)*60,p)
    def test_heldout_fixed_arms_complete(self):
        rng=np.random.default_rng(63)
        ids=tuple("original_line_"+str(i) for i in range(47))
        cls=tuple("BREAST" for _ in ids)
        y=rng.normal(size=47)
        a=rng.normal(size=(47,310))
        b=rng.normal(size=(47,125))
        p=rng.normal(size=(60,315))
        all_ids=ids+tuple("extra_"+str(i) for i in range(13))
        out=run(ids,cls,y,a,b,all_ids,("BREAST",)*60,p)
        self.assertEqual(out["n_matched_unique_cell_lines"],47)
        self.assertEqual(len(out["mae_same_source_GI50_units"]),11)
        self.assertEqual(out["independent_biological_observation_count"],47)
        self.assertFalse(out["same_contemporaneous_specimen"])
if __name__=="__main__":unittest.main()
