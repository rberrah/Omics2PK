"""Synthetic-only omics→drug-PK typed network evidence gates."""
import unittest
from omics2pk.physiopharmaco_cascade import (
    PKCascadeError,build_layered_pk_network,
)

class NetPKStyleMolecularCascade(unittest.TestCase):
    def test_unknown_metabolite_bin_has_no_false_transport_substrate_edge(self):
        graph=build_layered_pk_network("metformin",[
            {"source_id":"222169","omic":"metabolomics",
             "is_original_measured":True,"has_model_attribution":True,
             "model_run_url":"https://github.com/rberrah/Omics2PK/actions/runs/41912428806"}
        ],molecular_pk_endpoint="AUClast_0_24")
        self.assertEqual(graph["n_automatically_identified_unknown_metabolites"],0)
        self.assertTrue(any(n["id"]=="measured_metabolomics:222169" for n in graph["nodes"]))
        edges=[e for e in graph["edges"] if e["from"]=="measured_metabolomics:222169"]
        self.assertEqual(len(edges),1)
        self.assertEqual(edges[0]["evidence_level"],"model_association_not_chemistry")
        self.assertNotIn("SLC22A2",edges[0]["relation"])
    def test_source_measured_protein_can_join_preexisting_gene_node_not_function(self):
        graph=build_layered_pk_network("tacrolimus",[
          {"source_id":"UNIPROT:ABCB1_TEST","omic":"proteomics",
           "is_original_measured":True,"gene_symbol":"ABCB1"}])
        self.assertTrue(any(e["relation"]=="original_assay_mapping_to_gene"
                           and e["to"]=="gene_or_protein:ABCB1" for e in graph["edges"]))
        self.assertTrue(graph["no_clinical_function_inferred_from_STRING_SHAP_or_feature_correlations"])
    def test_string_edge_only_measured_original_gene_and_correct_type(self):
        inputs=[
          {"source_id":"RNA_IL18","omic":"transcriptomics","is_original_measured":True,
           "gene_symbol":"IL18"},
          {"source_id":"PROT_LGALS3","omic":"proteomics","is_original_measured":True,
           "gene_symbol":"LGALS3"}]
        graph=build_layered_pk_network("tacrolimus",inputs,string_interactions=[
          {"source":"IL18","target":"LGALS3","score":0.74,"network_type":"functional"}])
        matches=[e for e in graph["edges"] if e["relation"].startswith("STRING_")]
        self.assertEqual(len(matches),1)
        self.assertEqual(matches[0]["evidence_level"],"STRING_database_prior_not_PK_causality")
        with self.assertRaises(PKCascadeError):
            build_layered_pk_network("tacrolimus",inputs,string_interactions=[
              {"source":"IL18","target":"CYP3A4","score":0.999,"network_type":"functional"}])
        with self.assertRaises(PKCascadeError):
            build_layered_pk_network("tacrolimus",inputs,string_interactions=[
              {"source":"IL18","target":"LGALS3","score":1.5,"network_type":"physical"}])
    def test_unknown_mass_cannot_be_assigned_from_numeric_bin_label(self):
        with self.assertRaises(PKCascadeError):
            build_layered_pk_network("metformin",[
              {"source_id":"222169","omic":"metabolomics","is_original_measured":True,
               "chemical_id":"CHEBI:16797","chemical_grade":"source_name_only"}])
    def test_gene_rna_is_not_chemical_and_no_pk_when_endpoint_not_assayed(self):
        with self.assertRaises(PKCascadeError):
            build_layered_pk_network("metformin",[
              {"source_id":"RNA_SLC22A2","omic":"transcriptomics",
               "chemical_id":"CHEBI:123","is_original_measured":True}])
        with self.assertRaises(PKCascadeError):
            build_layered_pk_network("metformin",[
              {"source_id":"S1","omic":"metabolomics","is_original_measured":True,
               "has_model_attribution":True}],molecular_pk_endpoint="not_measured")
if __name__=="__main__":unittest.main()
