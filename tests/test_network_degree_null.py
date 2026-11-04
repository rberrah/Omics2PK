"""Explicitly limited degree-conditioned selected-network null, synthetic."""
import unittest
from omics2pk.network_degree_null import (
    NetworkNullError, degree_conditioned_bridge_null,
)

GENES = ("A","B","C","D","E","F")
SHAP = ("A","B","C")
PRIOR = ("D","E","F")
EDGES = [
    {"gene_A":"A","gene_B":"D","confidence_functional_association":0.9},
    {"gene_A":"A","gene_B":"B","confidence_functional_association":0.8},
    {"gene_A":"D","gene_B":"E","confidence_functional_association":0.7},
    {"gene_A":"C","gene_B":"F","confidence_functional_association":0.6},
]


class DegreeMatchedNetworkControls(unittest.TestCase):
    def test_condition_exact_degree_and_fixed_selected_nodes(self):
        result = degree_conditioned_bridge_null(GENES, SHAP, PRIOR, EDGES)
        self.assertEqual(result["n_conditioned_role_assignments"], 12)
        self.assertEqual(result["observed_shap_to_prior_edges"], 2)
        self.assertAlmostEqual(result["observed_shap_to_prior_edge_weight_sum"],1.5)
        self.assertEqual(sum(result["conditional_cross_role_edge_count_histogram"].values()),12)
        self.assertTrue(result["null_preserves_exact_node_degree_histogram_by_role"])
        self.assertTrue(result["null_uses_only_original_selected_genes_not_all_measured_omics"])
        self.assertFalse(result["direct_transport_or_pk_causality_demonstrated"])
    def test_empty_network_does_not_become_transport(self):
        r=degree_conditioned_bridge_null(GENES,SHAP,PRIOR,[])
        self.assertEqual(r["observed_shap_to_prior_edges"],0)
        self.assertEqual(r["expected_number_of_cross_role_edges_under_conditional_relabeling"],0)
        self.assertEqual(r["fraction_of_conditional_assignments_with_at_least_observed_bridges"],1)
    def test_source_and_score_guards(self):
        with self.assertRaises(NetworkNullError):
            degree_conditioned_bridge_null(GENES,SHAP,PRIOR,EDGES+EDGES[:1])
        with self.assertRaises(NetworkNullError):
            degree_conditioned_bridge_null(GENES,SHAP,PRIOR,[{**EDGES[0],"confidence_functional_association":1.5}])
        with self.assertRaises(NetworkNullError):
            degree_conditioned_bridge_null(GENES,("A","D"),PRIOR,EDGES)
        with self.assertRaises(NetworkNullError):
            degree_conditioned_bridge_null(GENES,SHAP,PRIOR,EDGES,max_exact_assignments=2)


if __name__=="__main__":
    unittest.main()
