"""Synthetic fold-safe mechanistic source-panel / random RNA contracts."""
import unittest
try:
    import numpy as np
    from omics2pk.mechanism_matched_control import (
        DRUG_GENES, FunctionalControlError,
        benchmark_prior_against_controls, matched_source_control_panels,
    )
except ImportError:  # lightweight mandatory offline suite: numerical dependencies optional
    np = None



@unittest.skipUnless(np is not None, "NumPy/scikit-learn optional scientific extra not installed")
class MatchedSourceGeneNull(unittest.TestCase):
    def setUp(self):
        rng=np.random.default_rng(146)
        self.x=rng.normal(size=(35,90))
        self.genes=DRUG_GENES+tuple("G"+str(i) for i in range(83))
        self.ids=tuple("cell_"+str(i) for i in range(35))
        self.lineage=tuple("A" if i%2 else "B" for i in range(35))
        self.y=0.1*self.x[:,0]+rng.normal(0,.1,size=35)

    def test_train_only_matched_controls_never_reuse_prior_gene(self):
        a=matched_source_control_panels(self.x[:27],tuple(range(7)),
                                        n_controls=5,pool_size=6,seed=88)
        b=matched_source_control_panels(self.x[:27],tuple(range(7)),
                                        n_controls=5,pool_size=6,seed=88)
        self.assertEqual(a,b)
        self.assertEqual(len(a),5)
        self.assertTrue(all(len(set(p))==7 and not(set(p)&set(range(7))) for p in a))
        # Test rows cannot influence the controls: they were never passed.
        self.x[27:]=1e7
        c=matched_source_control_panels(self.x[:27],tuple(range(7)),
                                        n_controls=5,pool_size=6,seed=88)
        self.assertEqual(a,c)

    def test_fixed_predeclared_prior_and_matched_random_loss_only(self):
        out=benchmark_prior_against_controls(self.ids,self.lineage,self.y,
                self.genes,self.x,repeats=1,n_controls=3)
        self.assertEqual(out["n_source_cell_lines"],35)
        self.assertEqual(out["n_predeclared_genes_measured"],7)
        self.assertEqual(out["n_control_pipelines"],3)
        self.assertTrue(out["same_original_cell_line_heldout_test_for_each_arm"])
        self.assertTrue(out["original_cell_line_GI50_is_not_human_PK"])
        self.assertFalse(out["mechanism_function_confirmed"])
        self.assertEqual(len(out["random_gene_control_MAE_5_50_95_percentiles"]),3)
        self.assertNotIn("cell_1",str(out))

    def test_incomplete_prior_cannot_claim_whitebox_function(self):
        with self.assertRaises(FunctionalControlError):
            benchmark_prior_against_controls(
                self.ids,self.lineage,self.y,self.genes[1:],
                self.x[:,1:],repeats=1,n_controls=2)
        with self.assertRaises(FunctionalControlError):
            matched_source_control_panels(self.x[:27],(0,0),n_controls=2)
        with self.assertRaises(FunctionalControlError):
            matched_source_control_panels(self.x[:27],(0,200),n_controls=2)


if __name__=="__main__":
    unittest.main()
