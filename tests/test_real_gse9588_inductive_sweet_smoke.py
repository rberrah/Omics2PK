"""All data fabricated, no original GEO samples or clinical donor labels."""
import importlib.util
import math
import unittest
from io import BytesIO

from scripts.smoke_real_gse9588_inductive_sweet import (
    TRAIN,PANEL,WEIGHT_PANEL,extract_train_variable_array_probes,
    extract_train_variable_array_probe_sets,real_gse_frozen_network_smoke,
    reference_indices
)

N=427

def toy(nrows=1100):
    payload=BytesIO()
    payload.write(b"!series_matrix_table_begin\n")
    payload.write(b"ID_REF\t"+b"\t".join(
        f"GSM{700000+i}".encode() for i in range(N))+b"\n")
    for j in range(nrows):
        nums=[str(round(3.0+(j%31)*0.1+0.9*math.sin((j+1)*(i+2)/47),6)).encode()
              for i in range(N)]
        payload.write(b"FAKE_PROBE\t"+b"\t".join(nums)+b"\n")
    payload.write(b"!series_matrix_table_end\n")
    payload.seek(0)
    return payload

class SourceFrozenReferenceTests(unittest.TestCase):
    def test_toy_427_sample_extracted_from_source_all_available(self):
        names,probes,n=extract_train_variable_array_probes(toy())
        self.assertEqual(len(names),N)
        self.assertEqual(len(probes),PANEL)
        self.assertEqual(n,1100)
        self.assertNotIn("GSM7000",repr([PANEL,n]))


    def test_toy_extracts_nested_edge_and_weight_probe_universes(self):
        names,edges,weights,n=extract_train_variable_array_probe_sets(toy())
        self.assertEqual(len(names),N)
        self.assertEqual(len(edges),PANEL)
        self.assertEqual(len(weights),WEIGHT_PANEL)
        self.assertEqual(edges,weights[:PANEL])
        self.assertEqual(n,1100)


    def test_hash_reference_selection_is_order_independent_and_outcome_free(self):
        gsm=[f"GSM{700000+i}" for i in range(N)]
        selected=reference_indices(gsm,"sha256_gsm75")
        self.assertEqual(len(selected),TRAIN)
        selected_ids={gsm[i] for i in selected}
        reversed_gsm=list(reversed(gsm))
        reversed_idx=reference_indices(reversed_gsm,"sha256_gsm75")
        self.assertEqual(selected_ids,{reversed_gsm[i] for i in reversed_idx})
        self.assertNotEqual(selected,list(range(TRAIN)))


    def test_hash_reference_sizes_are_nested_and_prespecified(self):
        gsm=[f"GSM{700000+i}" for i in range(N)]
        ids_by_n={}
        for n in (50,75,100,150):
            idx=reference_indices(gsm,f"sha256_gsm{n}")
            self.assertEqual(len(idx),n)
            ids_by_n[n]={gsm[i] for i in idx}
        self.assertTrue(ids_by_n[50] < ids_by_n[75])
        self.assertTrue(ids_by_n[75] < ids_by_n[100])
        self.assertTrue(ids_by_n[100] < ids_by_n[150])
        with self.assertRaises(ValueError):
            reference_indices(gsm,"sha256_gsm10")

    def test_hash_reference_extracts_probe_sets_using_its_training_only_variance(self):
        names,edges,weights,n=extract_train_variable_array_probe_sets(
            toy(),reference_strategy="sha256_gsm75")
        self.assertEqual(len(names),N)
        self.assertEqual(len(edges),PANEL)
        self.assertEqual(len(weights),WEIGHT_PANEL)
        self.assertEqual(edges,weights[:PANEL])
        self.assertEqual(n,1100)

    @unittest.skipUnless(importlib.util.find_spec("numpy"),"Optional numpy not installed")
    def test_toy_427_original_style_frozen_reference_does_not_claim_PK(self):
        ids,probes,n=extract_train_variable_array_probes(toy())
        result=real_gse_frozen_network_smoke(ids,probes,n)
        self.assertEqual(result["fixed_reference_training_GSM_samples"],TRAIN)
        self.assertEqual(result["unseen_GSM_samples_projected"],352)
        self.assertFalse(result["predictive_CYP_activity_test_executed"])
        self.assertFalse(result["equivalent_to_official_SWEET"])
        self.assertFalse(result["gene_symbol_annotation_provenance"])
        self.assertFalse(result["same_original_human_liver_donors_proven_independent"])
        self.assertFalse(result["clinical_or_invitro_PK_performance_claim"])
        self.assertEqual(result["unseen_GSM_samples_projected"] - result["projected_networks_distinct_hash_count"], result["n_projected_individual_graph_hash_collisions"])
        self.assertEqual(result["legacy_v0"]["unseen_GSM_samples_projected"],352)
        self.assertEqual(result["query_anchored_v1"]["unseen_GSM_samples_projected"],352)
        self.assertLessEqual(
            result["query_anchored_v1"]["n_projected_individual_graph_hash_collisions"],
            result["legacy_v0"]["n_projected_individual_graph_hash_collisions"])
        for key in ("legacy_v0","query_anchored_v1"):
            q=result[key]
            self.assertEqual(
                q["n_original_queries_low_clamped"]+
                q["n_original_queries_high_clamped"],
                q["n_original_queries_clamped_reference_weight"])
            self.assertLessEqual(
                q["min_weight_position_unclipped"],
                q["median_weight_position_unclipped"])
            self.assertLessEqual(
                q["median_weight_position_unclipped"],
                q["max_weight_position_unclipped"])
            self.assertLessEqual(q["min_similarity_z_vs_reference"],
                                 q["max_similarity_z_vs_reference"])
        self.assertEqual(
            result["v1_collision_reduction"],
            result["legacy_v0"]["n_projected_individual_graph_hash_collisions"]
            - result["query_anchored_v1"]["n_projected_individual_graph_hash_collisions"])
        self.assertTrue(result["v1_preserves_or_improves_graph_uniqueness"])
        self.assertFalse(result["network_features_sufficiently_person_specific_for_individual_PK"])
        self.assertIn("network_features_sufficiently_person_specific_for_individual_PK",result)
        self.assertNotIn("GSM700",str(result))


    @unittest.skipUnless(importlib.util.find_spec("numpy"),"Optional numpy not installed")
    def test_broader_weight_universe_preserves_no_pk_claim(self):
        ids,edges,weights,n=extract_train_variable_array_probe_sets(toy())
        result=real_gse_frozen_network_smoke(
            ids,edges,n,weight_probe_values=weights)
        self.assertEqual(result["sample_weight_expression_probe_count"],WEIGHT_PANEL)
        self.assertEqual(
            result["sample_weight_universe"],
            "broader_train_variable_probe_universe")
        self.assertFalse(result["predictive_CYP_activity_test_executed"])
        self.assertFalse(result["network_features_sufficiently_person_specific_for_individual_PK"])
        self.assertEqual(result["query_anchored_v1"]["unseen_GSM_samples_projected"],352)


    @unittest.skipUnless(importlib.util.find_spec("numpy"),"Optional numpy not installed")
    def test_hash_reference_smoke_preserves_no_pk_claim(self):
        ids,edges,weights,n=extract_train_variable_array_probe_sets(
            toy(),reference_strategy="sha256_gsm75")
        idx=reference_indices(ids,"sha256_gsm75")
        result=real_gse_frozen_network_smoke(
            ids,edges,n,weight_probe_values=weights,
            train_indices=idx,reference_strategy="sha256_gsm75")
        self.assertEqual(result["reference_selection_strategy"],"sha256_gsm75")
        self.assertFalse(result["reference_selection_uses_expression_values"])
        self.assertFalse(result["reference_selection_uses_PK_or_outcomes"])
        self.assertEqual(result["fixed_reference_training_GSM_samples"],TRAIN)
        self.assertEqual(result["query_anchored_v1"]["unseen_GSM_samples_projected"],352)
        self.assertFalse(result["predictive_CYP_activity_test_executed"])


    @unittest.skipUnless(importlib.util.find_spec("numpy"),"Optional numpy not installed")
    def test_hash_reference_50_has_dynamic_query_count(self):
        ids,edges,weights,n=extract_train_variable_array_probe_sets(
            toy(),reference_strategy="sha256_gsm50")
        idx=reference_indices(ids,"sha256_gsm50")
        result=real_gse_frozen_network_smoke(
            ids,edges,n,weight_probe_values=weights,
            train_indices=idx,reference_strategy="sha256_gsm50")
        self.assertEqual(result["fixed_reference_training_GSM_samples"],50)
        self.assertEqual(result["query_anchored_v1"]["unseen_GSM_samples_projected"],377)
        self.assertFalse(result["predictive_CYP_activity_test_executed"])

if __name__=="__main__":unittest.main()
