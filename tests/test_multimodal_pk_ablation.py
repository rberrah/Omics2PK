"""All benchmarks use fictional measurements. Never commit patient values."""
import copy
import json
import unittest

try:
    import numpy as np
    from omics2pk.multimodal_pk_ablation import (
        compare_modalities_pk, MultimodalPKBenchmarkError)
except ModuleNotFoundError as exc:
    if exc.name not in ("numpy","sklearn"):
        raise
    np=None


@unittest.skipIf(np is None,"Optional NumPy/scikit-learn benchmark not installed")
class MultimodalPatientPKTests(unittest.TestCase):
    def setUp(self):
        rng=np.random.default_rng(20261019)
        n=18
        self.d=dict(
            subject_ids=[f"FAKE_PERSON_{i:03d}" for i in range(n)],
            clinical=rng.normal(size=(n,3)),
            genotype_dosages=rng.integers(0,3,(n,2)),
            measured_pk=np.exp(rng.normal(2.8,0.2,n)),
            index_dose_at=["2025-01-03T09:00:00+00:00"]*n,
            pk_target_name="FICTIONAL_AUC", pk_target_unit="fake_area_units",
            source_kind="synthetic", source_sha256="a"*64,
            fixed_scores={"transcriptomics":{"FAKE_RNA0":0.8,"FAKE_RNA1":-0.2},
                          "proteomics":{"FAKE_PRO0":1}},
        )
        self.d["modalities"]={
            "transcriptomics":{
                "feature_ids":[f"FAKE_RNA{i}" for i in range(9)],
                "values":rng.normal(4,1,(n,9)),
                "collected_at":["2025-01-02T09:00:00Z"]*n,
                "tissue":"liver","value_type":"normalized_log2_expression"},
            "proteomics":{
                "feature_ids":[f"FAKE_PRO{i}" for i in range(7)],
                "values":rng.normal(6,1,(n,7)),
                "collected_at":["2025-01-02T08:00:00Z"]*n,
                "tissue":"plasma","value_type":"normalized_protein_abundance"},
            "metabolomics":{
                "feature_ids":[f"FAKE_MS1:{i}" for i in range(5)],
                "values":rng.uniform(0,1e5,(n,5)),
                "collected_at":["2025-01-02T09:00:00Z"]*n,
                "tissue":"blood_plasma","value_type":"ms1_peak_area"},
        }

    def test_same_patient_same_fold_all_methods_and_privacy(self):
        r=compare_modalities_pk(**self.d)
        self.assertEqual(r["primary_reference"],"clinical_pgx")
        self.assertEqual(r["n_independent_people"],18)
        self.assertEqual(r["status"],"SYNTHETIC_METHOD_CONTRACT_NO_HUMAN_PK_EVIDENCE")
        self.assertIn("early_fusion_pca2",r["models"])
        self.assertIn("late_fusion_fixed_average",r["models"])
        self.assertIn("transcriptomics_fixed_score",r["models"])
        self.assertIn("proteomics_fixed_score",r["models"])
        self.assertNotIn("metabolomics_fixed_score",r["models"])
        self.assertTrue(all(x["n_heldout_biological_people"]==18
                            for x in r["models"].values()))
        self.assertTrue(all(np.isfinite(x["mae"]) for x in r["models"].values()))
        self.assertNotIn("FAKE_PERSON",json.dumps(r))
        self.assertNotIn("FAKE_RNA0",json.dumps(r))
        self.assertNotIn(str(self.d["measured_pk"][0]),json.dumps(r))
        self.assertFalse(r["external_validation"])

    def test_exact_order_independence_after_patient_matching(self):
        before=compare_modalities_pk(**self.d)
        rng=np.random.default_rng(16)
        perm=rng.permutation(len(self.d["subject_ids"]))
        d=copy.deepcopy(self.d)
        for name in ("subject_ids","clinical","genotype_dosages",
                     "measured_pk","index_dose_at"):
            d[name]=np.asarray(d[name])[perm].tolist()
        for spec in d["modalities"].values():
            spec["values"]=np.asarray(spec["values"])[perm].tolist()
            spec["collected_at"]=np.asarray(spec["collected_at"])[perm].tolist()
        after=compare_modalities_pk(**d)
        self.assertEqual(before,after)

    def test_one_modality_no_misleading_fusion(self):
        d=copy.deepcopy(self.d)
        d["modalities"]={"metabolomics":d["modalities"]["metabolomics"]}
        d["fixed_scores"]={}
        r=compare_modalities_pk(**d)
        self.assertNotIn("early_fusion_pca2",r["models"])
        self.assertNotIn("late_fusion_fixed_average",r["models"])
        self.assertIn("metabolomics_pca2",r["models"])

    def test_missing_or_fabricated_genotypes_rejected(self):
        d=copy.deepcopy(self.d)
        d["genotype_dosages"][0,1]=4
        with self.assertRaisesRegex(MultimodalPKBenchmarkError,"0/1/2"):
            compare_modalities_pk(**d)

    def test_raw_rna_counts_not_relabelled_normalized(self):
        d=copy.deepcopy(self.d)
        d["modalities"]["transcriptomics"]["value_type"]="raw_counts"
        with self.assertRaisesRegex(MultimodalPKBenchmarkError,"Unnormalized RNA"):
            compare_modalities_pk(**d)

    def test_postdose_molecular_data_rejected_before_predict(self):
        d=copy.deepcopy(self.d)
        d["modalities"]["proteomics"]["collected_at"][1]="2025-01-04T00:00:00+00:00"
        with self.assertRaisesRegex(MultimodalPKBenchmarkError,"Post-landmark"):
            compare_modalities_pk(**d)

    def test_time_zone_cannot_be_inferred_from_T0(self):
        d=copy.deepcopy(self.d)
        d["modalities"]["metabolomics"]["collected_at"][0]="T0"
        with self.assertRaisesRegex(MultimodalPKBenchmarkError,"ISO"):
            compare_modalities_pk(**d)

    def test_unknown_score_and_mixed_names_rejected(self):
        d=copy.deepcopy(self.d)
        d["fixed_scores"]={"transcriptomics":{"FAKE_NOT_MEASURED":1}}
        with self.assertRaisesRegex(MultimodalPKBenchmarkError,"measured feature"):
            compare_modalities_pk(**d)

    def test_xgboost_sample_size_gate(self):
        d=copy.deepcopy(self.d)
        d["include_xgboost"]=True
        with self.assertRaisesRegex(MultimodalPKBenchmarkError,">=80"):
            compare_modalities_pk(**d)

    def test_no_genetic_source_status_not_genotype_justification(self):
        d=copy.deepcopy(self.d)
        d["genotype_dosages"]=None
        r=compare_modalities_pk(**d)
        self.assertEqual(r["primary_reference"],"clinical")
        self.assertFalse(r["genotype_available"])


if __name__=="__main__":
    unittest.main()
