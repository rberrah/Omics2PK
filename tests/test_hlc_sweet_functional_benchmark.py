"""Synthetic qualification for the site-disjoint HLC functional benchmark."""
import unittest

try:
    import numpy as np
    from omics2pk.hlc_sweet_functional_benchmark import (
        HLCFunctionalBenchmarkError,
        evaluate_site_disjoint_hlc,
    )
except ModuleNotFoundError as exc:
    if exc.name!="numpy":
        raise
    np=None


@unittest.skipIf(np is None,"NumPy optional extra not installed")
class HLCSiteDisjointBenchmarkTests(unittest.TestCase):
    def setUp(self):
        rng=np.random.default_rng(20261022)
        self.genes=tuple(["CYP3A4"]+[f"G{i:02}" for i in range(11)])
        self.ids=tuple(f"D{i:02}" for i in range(18))
        self.sites=tuple(["A"]*6+["B"]*6+["C"]*6)
        self.x=rng.normal(size=(len(self.genes),len(self.ids)))
        # Make target transcript genuinely predictive while preserving site shifts.
        site_shift=np.array([0.0]*6+[0.8]*6+[-0.6]*6)
        self.x[0]+=site_shift
        latent=1.2*self.x[0]+0.35*self.x[1]-0.2*self.x[2]
        self.y=np.exp(2.0+0.28*latent)+rng.normal(0,0.3,size=18)
        self.y=np.maximum(self.y,0.01)
        self.edges=self.genes[:5]

    def run_benchmark(self,x=None,y=None):
        return evaluate_site_disjoint_hlc(
            self.x if x is None else x,
            genes=self.genes,subject_ids=self.ids,sites=self.sites,
            activity=self.y if y is None else y,
            target_gene="CYP3A4",edge_genes=self.edges,
            pca_components=3,
        )

    def test_site_disjoint_protocol_and_privacy_flags(self):
        out=self.run_benchmark()
        self.assertEqual(out["n_subjects"],18)
        self.assertEqual(out["n_sites"],3)
        self.assertEqual(set(out["site_results"]),{"A","B","C"})
        self.assertTrue(out["preprocessing_and_network_fit_on_training_sites_only"])
        self.assertTrue(out["training_SWEET_features_cross_fitted"])
        self.assertFalse(out["pk_pd_target_used_for_representation_learning"])
        self.assertTrue(out["not_clinical_PK_or_PopPK"])
        self.assertFalse(out["individual_donor_predictions_exported"])
        self.assertNotIn("D00",str(out))

    def test_each_outer_fold_holds_complete_site(self):
        out=self.run_benchmark()
        for site,row in out["site_results"].items():
            self.assertEqual(row["n_heldout_donors"],6)
            self.assertEqual(row["n_training_donors"],12)
            self.assertEqual(len(row["sweet_pipeline_fingerprint"]),64)

    def test_held_site_activity_cannot_change_other_site_pipeline_fingerprint(self):
        a=self.run_benchmark()
        y=self.y.copy()
        y[:6]*=1000
        b=self.run_benchmark(y=y)
        # Fold holding site A has identical molecular training data and therefore
        # identical SWEET representation despite arbitrary held-out outcomes.
        self.assertEqual(
            a["site_results"]["A"]["sweet_pipeline_fingerprint"],
            b["site_results"]["A"]["sweet_pipeline_fingerprint"])
        # Representation learning never uses outcomes anywhere.
        for site in ("A","B","C"):
            self.assertEqual(
                a["site_results"][site]["sweet_pipeline_fingerprint"],
                b["site_results"][site]["sweet_pipeline_fingerprint"])

    def test_held_site_expression_cannot_change_its_training_pipeline(self):
        a=self.run_benchmark()
        changed=self.x.copy()
        changed[:,0:6]=changed[:,0:6]*100+50
        b=self.run_benchmark(x=changed)
        self.assertEqual(
            a["site_results"]["A"]["sweet_pipeline_fingerprint"],
            b["site_results"]["A"]["sweet_pipeline_fingerprint"])

    def test_invalid_site_design_fails_closed(self):
        with self.assertRaisesRegex(HLCFunctionalBenchmarkError,"three sites"):
            evaluate_site_disjoint_hlc(
                self.x,genes=self.genes,subject_ids=self.ids,
                sites=("A",)*9+("B",)*9,activity=self.y,
                target_gene="CYP3A4",edge_genes=self.edges)
        bad=list(self.sites); bad[-1]="D"
        with self.assertRaisesRegex(HLCFunctionalBenchmarkError,"three donors"):
            evaluate_site_disjoint_hlc(
                self.x,genes=self.genes,subject_ids=self.ids,
                sites=tuple(bad),activity=self.y,
                target_gene="CYP3A4",edge_genes=self.edges)

    def test_missing_target_or_bad_edge_panel_fails_closed(self):
        with self.assertRaisesRegex(HLCFunctionalBenchmarkError,"Target CYP"):
            evaluate_site_disjoint_hlc(
                self.x,genes=self.genes,subject_ids=self.ids,sites=self.sites,
                activity=self.y,target_gene="CYP2D6",edge_genes=self.edges)
        with self.assertRaisesRegex(HLCFunctionalBenchmarkError,"Edge genes"):
            evaluate_site_disjoint_hlc(
                self.x,genes=self.genes,subject_ids=self.ids,sites=self.sites,
                activity=self.y,target_gene="CYP3A4",
                edge_genes=("CYP3A4","CYP3A4","G01"))


if __name__=="__main__":
    unittest.main()
