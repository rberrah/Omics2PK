"""Synthetic controls on RNA, protein and chemical-name hypotheses by TP53."""
import unittest
try:
    import numpy as np
    from omics2pk.mutation_mechanism_controls import (
        PRIOR,MutationMechanismError,compare_source_priors_by_mutation
    )
except ImportError:
    np=None


@unittest.skipUnless(np is not None,"NumPy/scikit-learn optional scientific extra")
class MutationConditionalMechanism(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng=np.random.default_rng(900)
        cls.ids=tuple("C"+str(i) for i in range(35))
        cls.lineage=tuple(["L" if i%2 else "B" for i in range(35)])
        cls.mut=np.column_stack((np.array([0,1]*17+[1],dtype=float),
                                 rng.integers(0,2,35),rng.integers(0,40,35)))
        cls.labels=("mut_TP53","mut_TOP2B","mutation_burden")
        cls.names={k:tuple(v)+tuple(k+"_CTRL_"+str(i) for i in range(50))
                   for k,v in PRIOR.items()}
        cls.assays={k:rng.normal(size=(35,len(v))) for k,v in cls.names.items()}
        cls.y=rng.normal(size=35)

    def test_all_original_modalities_prior_and_mutation_aware_controls(self):
        r=compare_source_priors_by_mutation(self.ids,self.lineage,self.y,self.mut,
                                           self.labels,self.names,self.assays,
                                           repeats=1,controls=3)
        self.assertEqual(r["n_original_biologically_distinct_cell_lines"],35)
        self.assertEqual(set(r["matched_random_control_comparison_by_layer"]),{"rna","prot","met"})
        self.assertTrue(r["mutations_and_all_assay_preprocessing_train_only"])
        self.assertFalse(r["original_human_drug_pharmacokinetics_assayed"])
        self.assertFalse(r["confirmatory_mechanistic_or_genotype_effect_claim"])
        self.assertEqual(set(r["exploratory_TP53_stratum_MAE"]),{"0","1"})

    def test_missing_real_source_protein_is_not_filled(self):
        p={**self.names,"prot":self.names["prot"][1:]}
        a={**self.assays,"prot":self.assays["prot"][:,1:]}
        with self.assertRaises(MutationMechanismError):
            compare_source_priors_by_mutation(self.ids,self.lineage,self.y,
                self.mut,self.labels,p,a,repeats=1,controls=2)

    def test_duplicated_original_feature_identity_abstains(self):
        bad={**self.names,"rna":tuple(list(self.names["rna"][:1])*2+
            list(self.names["rna"][2:]))}
        with self.assertRaises(MutationMechanismError):
            compare_source_priors_by_mutation(self.ids,self.lineage,self.y,
                self.mut,self.labels,bad,self.assays,repeats=1,controls=2)

if __name__=="__main__":unittest.main()
