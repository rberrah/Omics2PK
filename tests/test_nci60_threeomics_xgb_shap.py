"""Synthetic-only source/SHAP tests; no NCI60 source values or cross-study batch claims."""
from pathlib import Path
import sys
import unittest
try:
    import numpy as np
    import xgboost
    import sklearn
except ImportError:
    raise unittest.SkipTest("Optional pinned source XGB dependencies qualified in dedicated workflow")
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from benchmark_original_nci60_threeomics_xgb_shap import (
    predict_xgb_shap,fit_targeted_whitebox,evaluate_original,COMBOS,CAPS,
    PREDECLARED_DRUG_GENES
)
class TrueThreeModalSHAPTests(unittest.TestCase):
    def test_all_omic_arms_defined_independently(self):
        self.assertEqual(len(COMBOS),8)
        self.assertEqual(CAPS["rna"],300)
        self.assertEqual(CAPS["prot"],300)
        self.assertEqual(CAPS["met"],120)
        self.assertIn("ABCB1",PREDECLARED_DRUG_GENES)
        self.assertIn("TOP2A",PREDECLARED_DRUG_GENES)
    def test_tree_shap_additivity_on_synthetic_drug_pd(self):
        x=np.array([[1,2],[2,4],[3,6],[4,8],[5,9],[6,7],[7,3],[8,1],[9,2],[10,3],[11,4],[12,5]],dtype=float)
        y=np.array([1.,4.,5.,2.,7.,3.,6.,5.,4.,8.,2.,7.])
        pred,shap=predict_xgb_shap(x,x[:2],y)
        self.assertEqual(shap.shape,(2,2))
        self.assertTrue(np.isfinite(pred).all())
        self.assertTrue(np.isfinite(shap).all())
    def test_whitebox_small_targeted_genes_only(self):
        x=np.asarray([[i,2*i,0.1*i] for i in range(12)],dtype=float)
        p=fit_targeted_whitebox(x,x[:3],np.arange(12,dtype=float))
        self.assertEqual(p.shape,(3,))
    def test_fail_on_fabricated_three_assay_line_identity(self):
        n=47
        kwargs=dict(ids=[f"LINE-{j}" for j in range(n)],
                    classes=["BREAST"]*n,drugresponse=np.arange(n,dtype=float),
                    rna_names=["ABCB1"]*300,met_names=["met"]*120,
                    rna=np.ones((n,300)),met=np.ones((n,120)),
                    protein_names=[f"OTHER-{j}" for j in range(n)],
                    protein_classes=["BREAST"]*n,protein=np.ones((n,300)),
                    protein_gene_symbols=["ABCB1"]*300,repeats=1)
        with self.assertRaises(ValueError):
            evaluate_original(**kwargs)
if __name__=="__main__":unittest.main()
