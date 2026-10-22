"""Synthetic qualification of the prospective frozen SWEET-inspired artifact."""
import json
import tempfile
import unittest
from pathlib import Path

try:
    import numpy as np
    from omics2pk.frozen_sweet_pipeline import (
        FrozenSweetPipelineError,
        _fingerprint_pipeline,
        fit_frozen_inductive_sweet_pipeline,
        frozen_pipeline_manifest,
        load_frozen_inductive_sweet_pipeline,
        predict_frozen_inductive_sweet_pipeline,
        save_frozen_inductive_sweet_pipeline,
    )
    from omics2pk.inductive_sweet import CONSTRUCTION_V1
    from omics2pk.transcriptome_preprocessing import (
        EDGE_SELECTION_PRESPECIFIED,
        EDGE_SELECTION_TRAIN_VARIANCE,
        MODE_PASSTHROUGH,
        MODE_TRAIN_ZSCORE,
        build_preprocessed_inductive_outer_split,
    )
except ModuleNotFoundError as exc:
    if exc.name != "numpy":
        raise
    np = None
    MODE_PASSTHROUGH = None
    MODE_TRAIN_ZSCORE = None


@unittest.skipIf(np is None, "NumPy optional extra not installed")
class FrozenSweetPipelineTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(20261022)
        self.genes = tuple(f"G{i:02}" for i in range(24))
        self.edge_genes = self.genes[:8]
        self.train_ids = tuple(f"TRAIN{i:02}" for i in range(12))
        self.x = rng.normal(size=(len(self.genes), len(self.train_ids)))
        self.query_a = rng.normal(size=len(self.genes))
        self.query_b = rng.normal(size=len(self.genes))

    def fit(self, *, mode=MODE_PASSTHROUGH):
        return fit_frozen_inductive_sweet_pipeline(
            self.x,
            genes=self.genes,
            subject_ids=self.train_ids,
            edge_genes=self.edge_genes,
            tissue="liver",
            mode=mode,
        )

    def test_fit_freezes_complete_training_artifact_and_defaults_to_v2(self):
        pipe = self.fit()
        self.assertEqual(pipe.reference.subject_ids, self.train_ids)
        self.assertEqual(pipe.genes, self.genes)
        self.assertEqual(pipe.edge_genes, self.edge_genes)
        self.assertEqual(pipe.construction, CONSTRUCTION_V1)
        self.assertEqual(pipe.reference.weight_genes, self.genes)
        self.assertEqual(pipe.preprocessor.train_subject_ids, self.train_ids)
        self.assertEqual(pipe.fingerprint_version, 2)
        self.assertEqual(pipe.edge_selection_mode, EDGE_SELECTION_PRESPECIFIED)
        self.assertEqual(pipe.edge_selection_top_k, len(self.edge_genes))
        manifest = frozen_pipeline_manifest(pipe)
        self.assertEqual(manifest["schema"], "omics2pk.frozen-inductive-sweet.v2")
        self.assertTrue(manifest["training_only_fit"])
        self.assertTrue(manifest["edge_selection_fitted_on_training_only"])
        self.assertFalse(manifest["pk_pd_target_used"])
        self.assertEqual(manifest["pipeline_fingerprint"], pipe.fingerprint)
        self.assertEqual(
            manifest["preprocessing_fingerprint"], pipe.preprocessor.fingerprint
        )
        self.assertEqual(
            manifest["reference_fingerprint"], pipe.reference.fingerprint
        )

    def test_prediction_does_not_refit_or_depend_on_another_query(self):
        pipe = self.fit()
        first = predict_frozen_inductive_sweet_pipeline(
            pipe, self.query_a, subject_id="NEW_A", genes=self.genes
        )
        _ = predict_frozen_inductive_sweet_pipeline(
            pipe, self.query_b, subject_id="NEW_B", genes=self.genes
        )
        second = predict_frozen_inductive_sweet_pipeline(
            pipe, self.query_a, subject_id="NEW_A", genes=self.genes
        )
        self.assertEqual(first, second)
        self.assertEqual(first["pipeline_fingerprint"], pipe.fingerprint)
        self.assertEqual(
            first["reference_fingerprint"], pipe.reference.fingerprint
        )
        self.assertFalse(first["preprocessor_refit_for_query"])
        self.assertFalse(first["reference_refit_for_query"])
        self.assertFalse(first["other_query_subjects_used"])
        self.assertFalse(first["edge_selector_refit_for_query"])
        self.assertEqual(first["edge_selection_fingerprint"],
                         pipe.edge_selection_fingerprint)

    def test_query_b_values_cannot_change_query_a_prediction(self):
        pipe = self.fit()
        a = predict_frozen_inductive_sweet_pipeline(
            pipe, self.query_a, subject_id="NEW_A", genes=self.genes
        )
        extreme_b = self.query_b * 1000 + np.linspace(
            -500, 500, len(self.genes)
        )
        _ = predict_frozen_inductive_sweet_pipeline(
            pipe, extreme_b, subject_id="NEW_B", genes=self.genes
        )
        a_again = predict_frozen_inductive_sweet_pipeline(
            pipe, self.query_a, subject_id="NEW_A", genes=self.genes
        )
        self.assertEqual(a, a_again)

    def test_training_subject_cannot_be_reused_as_new_query(self):
        pipe = self.fit()
        with self.assertRaisesRegex(FrozenSweetPipelineError, "new"):
            predict_frozen_inductive_sweet_pipeline(
                pipe,
                self.query_a,
                subject_id=self.train_ids[0],
                genes=self.genes,
            )

    def test_gene_order_and_single_subject_shape_are_strict(self):
        pipe = self.fit()
        with self.assertRaisesRegex(FrozenSweetPipelineError, "Gene identity/order"):
            predict_frozen_inductive_sweet_pipeline(
                pipe,
                self.query_a,
                subject_id="NEW",
                genes=tuple(reversed(self.genes)),
            )
        with self.assertRaisesRegex(FrozenSweetPipelineError, "exactly one"):
            predict_frozen_inductive_sweet_pipeline(
                pipe,
                np.column_stack((self.query_a, self.query_b)),
                subject_id="NEW",
                genes=self.genes,
            )

    def test_edge_panel_must_be_prespecified_unique_subset(self):
        with self.assertRaisesRegex(FrozenSweetPipelineError, "unique"):
            fit_frozen_inductive_sweet_pipeline(
                self.x,
                genes=self.genes,
                subject_ids=self.train_ids,
                edge_genes=(self.genes[0], self.genes[0], self.genes[1]),
                tissue="liver",
                mode=MODE_PASSTHROUGH,
            )
        with self.assertRaisesRegex(FrozenSweetPipelineError, "subset"):
            fit_frozen_inductive_sweet_pipeline(
                self.x,
                genes=self.genes,
                subject_ids=self.train_ids,
                edge_genes=(self.genes[0], self.genes[1], "NOT_IN_TRAIN"),
                tissue="liver",
                mode=MODE_PASSTHROUGH,
            )

    def test_train_zscore_query_imputation_is_frozen(self):
        x = self.x.copy()
        x[3, 2] = np.nan
        pipe = fit_frozen_inductive_sweet_pipeline(
            x,
            genes=self.genes,
            subject_ids=self.train_ids,
            edge_genes=self.edge_genes,
            tissue="liver",
            mode=MODE_TRAIN_ZSCORE,
        )
        query = self.query_a.copy()
        query[3] = np.nan
        prediction = predict_frozen_inductive_sweet_pipeline(
            pipe, query, subject_id="NEW", genes=self.genes
        )
        self.assertEqual(
            prediction["preprocessing_fingerprint"],
            pipe.preprocessor.fingerprint,
        )
        self.assertFalse(prediction["preprocessor_refit_for_query"])
        self.assertTrue(
            all(np.isfinite(v) for v in prediction["features"].values())
        )

    def test_prospective_prediction_is_identical_to_outer_split_heldout_path(self):
        ids = self.train_ids + ("TEST_A", "TEST_B")
        values = np.column_stack((self.x, self.query_a, self.query_b))
        outer = build_preprocessed_inductive_outer_split(
            values,
            genes=self.genes,
            edge_genes=self.edge_genes,
            subject_ids=ids,
            train_ids=self.train_ids,
            test_ids=("TEST_A", "TEST_B"),
            tissue="liver",
            mode=MODE_PASSTHROUGH,
            construction=CONSTRUCTION_V1,
        )
        pipe = self.fit()
        prospective = predict_frozen_inductive_sweet_pipeline(
            pipe, self.query_a, subject_id="TEST_A", genes=self.genes
        )
        heldout = outer["test"][0]
        self.assertEqual(
            pipe.preprocessor.fingerprint, outer["preprocessing_fingerprint"]
        )
        self.assertEqual(
            pipe.reference.fingerprint, outer["outer_reference_fingerprint"]
        )
        for gene in self.edge_genes:
            self.assertAlmostEqual(
                prospective["features"][gene], heldout["features"][gene], places=12
            )
        self.assertAlmostEqual(
            prospective["network_qc"]["weight"], heldout["weight"], places=12
        )
        self.assertAlmostEqual(
            prospective["network_qc"]["sample_similarity_z_vs_reference"],
            heldout["sample_similarity_z_vs_reference"],
            places=12,
        )
        self.assertAlmostEqual(
            prospective["network_qc"]["weight_position_unclipped"],
            heldout["weight_position_unclipped"],
            places=12,
        )
        self.assertEqual(
            prospective["network_qc"]["weight_clamped_outside_reference"],
            heldout["weight_clamped_outside_reference"],
        )
        self.assertEqual(
            prospective["network_qc"]["weight_clamp_direction"],
            heldout["weight_clamp_direction"],
        )

    def test_train_variance_selector_is_frozen_in_prospective_artifact(self):
        x = self.x.copy()
        x[2, :] *= 25
        x[9, :] *= 18
        x[15, :] *= 12
        pipe = fit_frozen_inductive_sweet_pipeline(
            x,
            genes=self.genes,
            subject_ids=self.train_ids,
            edge_genes=None,
            tissue="liver",
            mode=MODE_PASSTHROUGH,
            edge_selection_mode=EDGE_SELECTION_TRAIN_VARIANCE,
            edge_top_k=8,
        )
        self.assertEqual(pipe.edge_selection_mode, EDGE_SELECTION_TRAIN_VARIANCE)
        self.assertEqual(pipe.edge_selection_top_k, 8)
        self.assertEqual(len(pipe.edge_genes), 8)
        self.assertTrue({"G02", "G09", "G15"}.issubset(set(pipe.edge_genes)))
        before = predict_frozen_inductive_sweet_pipeline(
            pipe, self.query_a, subject_id="NEW_A", genes=self.genes
        )
        _ = predict_frozen_inductive_sweet_pipeline(
            pipe, self.query_b * 1000, subject_id="NEW_B", genes=self.genes
        )
        after = predict_frozen_inductive_sweet_pipeline(
            pipe, self.query_a, subject_id="NEW_A", genes=self.genes
        )
        self.assertEqual(before, after)
        self.assertFalse(before["edge_selector_refit_for_query"])

    def test_train_variance_prospective_artifact_matches_outer_split_path(self):
        x = self.x.copy()
        x[3, :] *= 20
        x[11, :] *= 15
        ids = self.train_ids + ("TEST_A", "TEST_B")
        values = np.column_stack((x, self.query_a, self.query_b))
        outer = build_preprocessed_inductive_outer_split(
            values,
            genes=self.genes,
            edge_genes=None,
            subject_ids=ids,
            train_ids=self.train_ids,
            test_ids=("TEST_A", "TEST_B"),
            tissue="liver",
            mode=MODE_PASSTHROUGH,
            construction=CONSTRUCTION_V1,
            edge_selection_mode=EDGE_SELECTION_TRAIN_VARIANCE,
            edge_top_k=8,
        )
        pipe = fit_frozen_inductive_sweet_pipeline(
            x,
            genes=self.genes,
            subject_ids=self.train_ids,
            edge_genes=None,
            tissue="liver",
            mode=MODE_PASSTHROUGH,
            construction=CONSTRUCTION_V1,
            edge_selection_mode=EDGE_SELECTION_TRAIN_VARIANCE,
            edge_top_k=8,
        )
        self.assertEqual(pipe.edge_genes, outer["edge_genes"])
        self.assertEqual(
            pipe.edge_selection_fingerprint, outer["edge_selection_fingerprint"]
        )
        self.assertEqual(
            pipe.reference.fingerprint, outer["outer_reference_fingerprint"]
        )
        prospective = predict_frozen_inductive_sweet_pipeline(
            pipe, self.query_a, subject_id="TEST_A", genes=self.genes
        )
        heldout = outer["test"][0]
        for gene in pipe.edge_genes:
            self.assertAlmostEqual(
                prospective["features"][gene], heldout["features"][gene], places=12
            )
        self.assertAlmostEqual(
            prospective["network_qc"]["weight"], heldout["weight"], places=12
        )

    def test_train_variance_artifact_roundtrip_preserves_selector_and_prediction(self):
        x = self.x.copy()
        x[4, :] *= 30
        pipe = fit_frozen_inductive_sweet_pipeline(
            x,
            genes=self.genes,
            subject_ids=self.train_ids,
            tissue="liver",
            mode=MODE_PASSTHROUGH,
            edge_selection_mode=EDGE_SELECTION_TRAIN_VARIANCE,
            edge_top_k=7,
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "frozen"
            save_frozen_inductive_sweet_pipeline(pipe, root)
            manifest = json.loads((root / "manifest.json").read_text())
            self.assertEqual(
                manifest["schema"], "omics2pk.frozen-inductive-sweet.artifact.v2"
            )
            self.assertEqual(
                manifest["pipeline"]["edge_selection_mode"],
                EDGE_SELECTION_TRAIN_VARIANCE,
            )
            self.assertEqual(
                manifest["pipeline"]["edge_selection_fingerprint"],
                pipe.edge_selection_fingerprint,
            )
            loaded = load_frozen_inductive_sweet_pipeline(root)
            self.assertEqual(loaded.edge_genes, pipe.edge_genes)
            self.assertEqual(
                loaded.edge_selection_fingerprint, pipe.edge_selection_fingerprint
            )
            self.assertEqual(loaded.fingerprint, pipe.fingerprint)
            self.assertEqual(loaded.fingerprint_version, 2)
            before = predict_frozen_inductive_sweet_pipeline(
                pipe, self.query_a, subject_id="NEW", genes=self.genes
            )
            after = predict_frozen_inductive_sweet_pipeline(
                loaded, self.query_a, subject_id="NEW", genes=self.genes
            )
            self.assertEqual(before, after)

    def test_v2_selector_fingerprint_tampering_is_detected(self):
        pipe = self.fit()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "frozen"
            save_frozen_inductive_sweet_pipeline(pipe, root)
            path = root / "manifest.json"
            manifest = json.loads(path.read_text())
            manifest["pipeline"]["edge_selection_fingerprint"] = "0" * 64
            path.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(
                FrozenSweetPipelineError, "selector fingerprint|pipeline fingerprint"
            ):
                load_frozen_inductive_sweet_pipeline(root)

    def test_legacy_v1_prespecified_artifact_remains_loadable(self):
        pipe = self.fit()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "frozen"
            save_frozen_inductive_sweet_pipeline(pipe, root)
            path = root / "manifest.json"
            manifest = json.loads(path.read_text())
            old_fingerprint = _fingerprint_pipeline(
                genes=pipe.genes,
                edge_genes=pipe.edge_genes,
                tissue=pipe.tissue,
                preprocessing_mode=pipe.preprocessing_mode,
                construction=pipe.construction,
                k=pipe.k,
                preprocessor_fingerprint=pipe.preprocessor.fingerprint,
                reference_fingerprint=pipe.reference.fingerprint,
            )
            manifest["schema"] = "omics2pk.frozen-inductive-sweet.artifact.v1"
            manifest["pipeline"]["schema"] = "omics2pk.frozen-inductive-sweet.v1"
            manifest["pipeline"]["pipeline_fingerprint"] = old_fingerprint
            for key in (
                "edge_selection_mode",
                "edge_selection_top_k",
                "edge_selection_fingerprint",
                "edge_selection_fitted_on_training_only",
            ):
                manifest["pipeline"].pop(key, None)
            path.write_text(json.dumps(manifest))
            loaded = load_frozen_inductive_sweet_pipeline(root)
            self.assertEqual(loaded.fingerprint_version, 1)
            self.assertEqual(
                loaded.edge_selection_mode, EDGE_SELECTION_PRESPECIFIED
            )
            self.assertEqual(loaded.edge_genes, pipe.edge_genes)
            before = predict_frozen_inductive_sweet_pipeline(
                pipe, self.query_a, subject_id="NEW", genes=self.genes
            )
            after = predict_frozen_inductive_sweet_pipeline(
                loaded, self.query_a, subject_id="NEW", genes=self.genes
            )
            self.assertNotEqual(before["pipeline_fingerprint"],
                                after["pipeline_fingerprint"])
            for gene in pipe.edge_genes:
                self.assertAlmostEqual(
                    before["features"][gene], after["features"][gene], places=12
                )

    def test_persisted_artifact_roundtrip_preserves_prediction_and_fingerprints(self):
        pipe = self.fit(mode=MODE_TRAIN_ZSCORE)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "frozen"
            saved = save_frozen_inductive_sweet_pipeline(pipe, root)
            manifest = json.loads((root / "manifest.json").read_text())
            self.assertTrue(manifest["contains_training_derived_molecular_data"])
            self.assertTrue(manifest["contains_training_identifiers"])
            self.assertTrue(manifest["redistribution_requires_source_permission"])
            self.assertFalse(manifest["arrays"]["allow_pickle"])
            self.assertEqual(saved["pipeline_fingerprint"], pipe.fingerprint)
            with np.load(root / "arrays.npz", allow_pickle=False) as arrays:
                self.assertEqual(arrays["reference_weight_expression"].dtype.kind, "f")

            loaded = load_frozen_inductive_sweet_pipeline(root)
            self.assertEqual(loaded.fingerprint, pipe.fingerprint)
            self.assertEqual(
                loaded.preprocessor.fingerprint, pipe.preprocessor.fingerprint
            )
            self.assertEqual(
                loaded.reference.fingerprint, pipe.reference.fingerprint
            )
            before = predict_frozen_inductive_sweet_pipeline(
                pipe, self.query_a, subject_id="NEW", genes=self.genes
            )
            after = predict_frozen_inductive_sweet_pipeline(
                loaded, self.query_a, subject_id="NEW", genes=self.genes
            )
            self.assertEqual(before["pipeline_fingerprint"], after["pipeline_fingerprint"])
            for gene in self.edge_genes:
                self.assertAlmostEqual(
                    before["features"][gene], after["features"][gene], places=12
                )
            self.assertAlmostEqual(
                before["network_qc"]["weight"],
                after["network_qc"]["weight"],
                places=12,
            )

    def test_persisted_artifact_rejects_array_tampering(self):
        pipe = self.fit()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "frozen"
            save_frozen_inductive_sweet_pipeline(pipe, root)
            arrays = root / "arrays.npz"
            arrays.write_bytes(arrays.read_bytes() + b"tamper")
            with self.assertRaisesRegex(FrozenSweetPipelineError, "checksum"):
                load_frozen_inductive_sweet_pipeline(root)

    def test_persisted_artifact_rejects_manifest_fingerprint_tampering(self):
        pipe = self.fit()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "frozen"
            save_frozen_inductive_sweet_pipeline(pipe, root)
            path = root / "manifest.json"
            manifest = json.loads(path.read_text())
            manifest["pipeline"]["pipeline_fingerprint"] = "0" * 64
            path.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(FrozenSweetPipelineError, "pipeline fingerprint"):
                load_frozen_inductive_sweet_pipeline(root)

    def test_persist_refuses_nonempty_destination(self):
        pipe = self.fit()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "frozen"
            root.mkdir()
            (root / "unrelated.txt").write_text("do not overwrite")
            with self.assertRaisesRegex(FrozenSweetPipelineError, "must be empty"):
                save_frozen_inductive_sweet_pipeline(pipe, root)

    def test_pipeline_fingerprint_changes_when_frozen_protocol_changes(self):
        a = self.fit()
        b = fit_frozen_inductive_sweet_pipeline(
            self.x,
            genes=self.genes,
            subject_ids=self.train_ids,
            edge_genes=self.genes[1:9],
            tissue="liver",
            mode=MODE_PASSTHROUGH,
        )
        c = fit_frozen_inductive_sweet_pipeline(
            self.x,
            genes=self.genes,
            subject_ids=self.train_ids,
            edge_genes=self.edge_genes,
            tissue="liver",
            mode=MODE_PASSTHROUGH,
            k=0.2,
        )
        self.assertNotEqual(a.fingerprint, b.fingerprint)
        self.assertNotEqual(a.fingerprint, c.fingerprint)


if __name__ == "__main__":
    unittest.main()
