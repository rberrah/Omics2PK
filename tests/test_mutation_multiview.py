"""Synthetic-only mutation-aware per-omic selection and abstention guards."""
import unittest
try:
    import numpy as np
    from omics2pk.mutation_multiview import (
        MODES,MutationMultiViewError,evaluate_mutation_modality_budget,_view,
    )
except ImportError:
    np=None


@unittest.skipUnless(np is not None,"optional numpy/scikit-learn research extra")
class MutationViews(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng=np.random.default_rng(200)
        cls.n=35
        tp53=np.array([0,1]*17+[1],dtype=float)
        mutation_burden=rng.integers(0,15,35)
        cls.mutation=np.column_stack((tp53,rng.integers(0,2,35),mutation_burden))
        cls.labels=("mut_TP53","mut_TOP2A","mutation_burden")
        cls.ids=tuple("S"+str(i) for i in range(35))
        cls.lineage=tuple(("A","B","C")[i%3] for i in range(35))
        cls.assays={
            "rna":rng.normal(size=(35,40)),
            "prot":rng.normal(size=(35,35)),
            "met":rng.normal(size=(35,30)),
        }
        cls.y=2+cls.assays["rna"][:,0]*(1-tp53)+cls.assays["met"][:,0]*tp53+rng.normal(0,.1,35)

    def test_full_mutation_and_threeomics_nested_fold_disjoint(self):
        a=evaluate_mutation_modality_budget(
            self.ids,self.lineage,self.y,self.labels,self.mutation,self.assays,
            repeats=1)
        self.assertEqual(a["n_unique_biological_cell_lines"],35)
        self.assertEqual(set(a["original_in_vitro_GI50_observed_holdout"]["fixed"]),set(MODES))
        self.assertEqual(set(a["original_in_vitro_GI50_observed_holdout"]["fixed_without_original_somatic_mutation_features"]),set(MODES))
        self.assertEqual(set(a["original_in_vitro_GI50_observed_holdout"]["observed_mutation_incremental_delta_MAE_negative_is_better"]),set(MODES))
        self.assertEqual(set(a["original_in_vitro_GI50_observed_holdout"]["nested_TP53_guided_80_90"]),{"80","90"})
        objective=a["original_in_vitro_GI50_observed_holdout"]["nested_TP53_guided_80_90_using_RMSE_inner"]
        self.assertEqual(set(objective),{"80","90"})
        self.assertTrue(all(v["criterion_chosen_using_only_inner_train_validation"]=="RMSE"
                            and v["never_selected_from_outer_MAE_or_RMSE"]
                            for v in objective.values()))
        self.assertEqual(set(a["multi_metric_exploratory_heldout"]
                             ["nested_TP53_RMSE_inner_modality_selection"]),{"80","90"})
        self.assertIn("RMSE_selected_nested80_vs_MAE_selected_nested80",
                      a["prespecified_paired_heldout_MAE_contrasts"])
        self.assertTrue(a["genotype_subgroup_model_chosen_using_outer_response"] is False)
        self.assertTrue(a["no_externally_replicated_modality_by_genotype_interaction"])
        self.assertTrue(a["no_patient_or_original_cell_line_feature_rows_exported"])
        metrics=a["multi_metric_exploratory_heldout"]["mutation_adjusted_fixed"]
        self.assertEqual(set(metrics),set(MODES))
        self.assertAlmostEqual(
            metrics["mutation_rna"]["MAE_original_source_scale"],
            a["original_in_vitro_GI50_observed_holdout"]["fixed"]["mutation_rna"]["mae"])
        self.assertEqual(a["prespecified_paired_heldout_MAE_contrasts"]
                         ["mutations_added_to_RNA"]["n_independent_line_units_resampled"],35)
        self.assertTrue(a["metric_criterion_sensitivity_descriptive_not_model_selection"]
                        ["no_outer_test_metric_used_for_training_or_selection"])
        self.assertEqual(set(a["prespecified_paired_heldout_rank_Spearman_contrasts"]),
                         {"mutations_added_to_RNA","all_three_vs_mutation_plus_RNA",
                          "mutation_plus_RNA_vs_mutation_only"})
        for p in a["original_in_vitro_GI50_observed_holdout"]["fixed"].values():
            self.assertGreaterEqual(p["mae"],0)

    def test_fold_transform_does_not_use_outer_test_values(self):
        train=self.assays["rna"][:25]
        test=self.assays["rna"][25:].copy()
        first,_=_view(train,test,20)
        altered=test.copy();altered[:]=1e10
        second,_=_view(train,altered,20)
        np.testing.assert_array_equal(first,second)

    def test_no_pgx_fabrication_and_stratum_sizes(self):
        with self.assertRaises(MutationMultiViewError):
            evaluate_mutation_modality_budget(self.ids,self.lineage,self.y,
                ("mutation_burden",),self.mutation[:,-1:],self.assays,repeats=1)
        m=self.mutation.copy();m[:,0]=0
        with self.assertRaises(MutationMultiViewError):
            evaluate_mutation_modality_budget(
                self.ids,self.lineage,self.y,self.labels,m,self.assays,repeats=1)
        m=self.mutation.copy();m[0,0]=np.nan
        with self.assertRaises(MutationMultiViewError):
            evaluate_mutation_modality_budget(
                self.ids,self.lineage,self.y,self.labels,m,self.assays,repeats=1)

if __name__=="__main__":unittest.main()
