"""Artificial-only tests for nested original Jarmusch PGx vs MS1 workflow."""
import unittest
try:
    import numpy as np
    from scripts.benchmark_original_jarmusch_nested_pk_germline import (
        _nested_cv,NestedPKValidationError,
        _nested_predict_once,ALPHAS,OMIC_K)
except ModuleNotFoundError as err:
    if err.name not in ("numpy","sklearn","scipy","h5py"):
        raise
    np=None
    _nested_cv=None


@unittest.skipUnless(_nested_cv is not None,"Scientific optional extras not installed")
class JarmuschNestedResearchContract(unittest.TestCase):
    def setup_source(self):
        rng=np.random.default_rng(8)
        n=13
        x=np.abs(rng.normal(size=(n,30)))
        # Every synthetic individual is unique and x0 encodes phenotype.
        c=np.column_stack((rng.uniform(50,95,n),rng.integers(0,2,n)))
        g=np.array([["*1/*1" if i%2 else "*1/*22",
                     "*1/*3" if i%3 else "*3/*3"] for i in range(n)],object)
        logy=2+.5*np.log1p(x[:,0])+rng.normal(0,.03,n)
        y=np.exp(logy)
        candidate={
           "putative_steroid_bile_feature_indices":np.array([0,1],int),
           "spectrally_annotated_qc_feature_indices":np.array([0,1,2,3],int),
        }
        return x,c,g,y,candidate

    def test_baseline_and_incremental_original_design_no_people_exposed(self):
        x,c,g,y,candidates=self.setup_source()
        r=_nested_cv(x,c,g,y,candidates=candidates,source_label="artificial measured CLF")
        self.assertEqual(r["n_source_people"],13)
        for name,v in r["models"].items():
            self.assertEqual(v["n_original_independent_people"],13)
            self.assertTrue(v["no_source_individual_outcomes_or_predictions_printed"])
            self.assertEqual(sum(v["k_selected_counts"].values()),13)
            self.assertEqual(sum(v["alpha_selected_counts"].values()),13)
            self.assertNotIn("_losses_private",v)
        self.assertEqual(r,_nested_cv(x,c,g,y,candidates=candidates,source_label="artificial measured CLF"))
        self.assertIn("paired_delta_MAE_vs_clinical_PGx_negative_better",
                      r["models"]["PGx_plus_putative_CYP3A_steroid_bile"])

    def test_abstain_when_no_pathway_source(self):
        x,c,g,y,a=self.setup_source()
        a["putative_steroid_bile_feature_indices"]=np.array([],int)
        r=_nested_cv(x,c,g,y,candidates=a,source_label="artificial")
        self.assertEqual(r["models"]["PGx_plus_putative_CYP3A_steroid_bile"]["status"],
                         "ABSTAIN_NO_SOURCE_GROUNDED_FEATURES")
        self.assertEqual(r["models"]["age_sex_weight_PGx_only"]["status"],
                         "MEASURED_ORIGINAL_PATIENT_HELDOUT_NESTED")

    def test_folds_and_invalid_values_rejected(self):
        x,c,g,y,a=self.setup_source()
        with self.assertRaises(NestedPKValidationError):
            _nested_cv(x,c,g,np.ones(9),candidates=a,source_label="artificial")
        with self.assertRaises(NestedPKValidationError):
            _nested_cv(-x,c,g,y,candidates=a,source_label="artificial")
        with self.assertRaises(NestedPKValidationError):
            _nested_predict_once(x,c,g,np.log(y),np.arange(10),np.array([3]),
                                 alpha=1.,k=1,allowed=np.arange(3))
        with self.assertRaises(NestedPKValidationError):
            _nested_predict_once(x,c,g,np.log(y),np.arange(10),np.array([12]),
                                 alpha=.0001,k=1,allowed=np.arange(3))


if __name__=="__main__":
    unittest.main()
