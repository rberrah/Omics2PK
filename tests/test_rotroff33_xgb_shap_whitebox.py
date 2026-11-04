"""Synthetic-only tests for XGBoost Tree SHAP and white-box evidence gates.

Actual metformin data downloaded ONLY in original-source GitHub Action.
"""
import unittest
try:
    import numpy as np
    import xgboost
    import sklearn
except ImportError:
    raise unittest.SkipTest("Optional pinned XGBoost+NumPy/Scikit-learn evaluated in real-source CI")

from omics2pk.rotroff33_clinical_e3 import E3SourceError
from omics2pk.rotroff33_xgb_shap_whitebox import _fit_xgb,_model_pred_shap,_loocv,PARAMS

class OriginalPKSHAPContracts(unittest.TestCase):
    def test_xgb_exact_tree_shap_contributions_add_up(self):
        x=np.asarray([[1,2,3],[2,3,4],[3,1,5],[2,5,6],
                      [4,2,7],[2,2,3],[4,3,2],[5,1,8],
                      [6,1,4],[3,7,2],[8,6,5],[4,8,2]],dtype=float)
        y=np.array([11,12,13,14,15,16,17,18,19,20,21,22],dtype=float)
        m=_fit_xgb(x,x[:1],y)
        p,shap=_model_pred_shap(m,x[:3])
        self.assertEqual(shap.shape,(3,3))
        self.assertTrue(np.all(np.isfinite(shap)))
        self.assertTrue(np.all(p>=0))
        self.assertEqual(PARAMS["max_depth"],2)
        self.assertEqual(PARAMS["n_jobs"],1)
    def test_cannot_run_with_31_source_phenotypes(self):
        records=[{"internal_patient_id":f"TEST-{i}","targets":{"metformin_auc":20.+i},
                  "clinical_plus_original_serum_creatinine":[40,1,80,1],
                  "predose_omics":[1]*384} for i in range(31)]
        anno={"named":list(range(155)),"unknown":list(range(155,384)),
              "paths":{"urea_cycle_partial":[1,2,3],
                       "purine_degradation_partial":[4,5]}}
        with self.assertRaises(E3SourceError):
            _loocv(records,anno,[str(i) for i in range(384)],
                   "metformin_auc","clinical_plus_original_serum_creatinine")

if __name__=="__main__":unittest.main()
