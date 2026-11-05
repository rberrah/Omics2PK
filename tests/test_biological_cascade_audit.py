"""Original-source labels, no inferred protein/activity/chemical identity."""
import unittest
try:
    import numpy as np
    from omics2pk.biological_cascade_audit import (
        DIRECT,METS,CascadeAuditError,audit_original_measured_cascade)
except ImportError:
    np=None


@unittest.skipUnless(np is not None, "NumPy/scikit-learn optional")
class CascadeEvidenceTests(unittest.TestCase):
    def setUp(self):
        r=np.random.default_rng(459)
        self.ids=tuple("S"+str(i) for i in range(35))
        self.y=r.normal(size=35)
        self.mut=np.column_stack((
            np.array([0,1]*17+[1],dtype=float),
            np.array([0]*31+[1]*4,dtype=float),
            r.integers(0,80,size=35)))
        self.labels=("mut_TP53","mut_ABCB1","mutation_burden")
        self.source_names={
            "rna":tuple(DIRECT),
            "prot":("TOP2B","ABCB1","CBR1"),
            "met":tuple(METS),
        }
        self.assays={
            "rna":r.normal(size=(35,7)),
            "prot":r.normal(size=(35,3)),
            "met":r.normal(size=(35,len(METS))),
        }
        self.assays["prot"][0,0]=float("nan")

    def test_explicit_available_and_missing_direct_measured_paths(self):
        r=audit_original_measured_cascade(self.ids,self.y,self.mut,
                                          self.labels,self.source_names,self.assays)
        self.assertEqual(r["n_biologically_distinct_cell_lines"],35)
        self.assertTrue(r["biological_interpretation_nodes"]["ABCB1"]
                        ["RNA_measured_in_original_source"])
        self.assertTrue(r["biological_interpretation_nodes"]["ABCB1"]
                        ["SWATH_protein_measured_in_original_source"])
        self.assertEqual(r["biological_interpretation_nodes"]["ABCB1"]
                        ["mutation_event_subgroup_comparison_status"],
                         "ABSTAIN_INSUFFICIENT_MUTATION_CARRIERS_FOR_SUBGROUP")
        self.assertFalse(r["biological_interpretation_nodes"]["TOP2A"]
                         ["SWATH_protein_measured_in_original_source"])
        self.assertEqual(r["biological_interpretation_nodes"]["TOP2B"]
                         ["n_original_RNA_and_SWATH_both_observed"],34)
        self.assertFalse(r["direct_transcript_to_active_protein_or_metabolite_mechanism_proven"])
        self.assertFalse(r["source_metabolites_MSI_level1_structure_verified"])
        self.assertNotIn("S0",str(r))

    def test_absent_source_redox_label_is_not_filled(self):
        names={**self.source_names,"met":self.source_names["met"][:-1]}
        assay={**self.assays,"met":self.assays["met"][:,:-1]}
        with self.assertRaisesRegex(CascadeAuditError,"source label missing"):
            audit_original_measured_cascade(self.ids,self.y,self.mut,
                                            self.labels,names,assay)

    def test_gene_protein_misidentification_does_not_pass(self):
        names={**self.source_names,"prot":("TOP2B","ABCB1","ABCB1")}
        with self.assertRaisesRegex(CascadeAuditError,"Ambiguous original"):
            audit_original_measured_cascade(self.ids,self.y,self.mut,self.labels,
                                            names,self.assays)

    def test_original_rare_subgroup_cannot_supply_functional_mutation_effect(self):
        m=self.mut.copy()
        m[:,0]=0
        with self.assertRaisesRegex(CascadeAuditError,"Insufficient TP53"):
            audit_original_measured_cascade(self.ids,self.y,m,self.labels,
                                            self.source_names,self.assays)


if __name__=="__main__":unittest.main()
