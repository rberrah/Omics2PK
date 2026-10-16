"""Synthetic integrity tests for 48-line cross-study, real-endpoint benchmark."""
import csv
import tempfile
import unittest
from pathlib import Path

try:
    import numpy as np
    from scripts.benchmark_public_nci60_twoomics import (
        RNA_CAP, METAB_CAP, fit_view, load_dataset, evaluate
    )
except ModuleNotFoundError as exc:
    if exc.name not in {"numpy", "sklearn", "xgboost"}:
        raise
    np = None


@unittest.skipIf(np is None, "Optional benchmark dependencies are not installed")
class TwoOmicsTests(unittest.TestCase):
    def test_holdout_never_affects_training_normalization_selection_or_pca(self):
        rng = np.random.default_rng(17)
        data = rng.normal(size=(48, 320))
        data[:3, 0] = np.nan
        initial = fit_view(data[:36], data[36:], RNA_CAP)
        changed = data[36:].copy()
        changed[:, :25] += 100000
        after = fit_view(data[:36], changed, RNA_CAP)
        np.testing.assert_allclose(initial[0], after[0], atol=1e-12)
        np.testing.assert_allclose(initial[2], after[2], atol=1e-12)
        np.testing.assert_array_equal(initial[4], after[4])

    def test_missingness_and_train_fold_only_median(self):
        rng = np.random.default_rng(14)
        data = rng.normal(size=(48, 150))
        data[0, 0] = np.nan
        data[20:, 0] = 300
        data[:, 1] = np.nan
        data[:, 2] = 5
        parts = fit_view(data[:20], data[20:], METAB_CAP)
        self.assertTrue(np.isfinite(parts[0]).all())
        self.assertTrue(np.isfinite(parts[1]).all())
        self.assertNotIn(1, parts[4])
        self.assertNotIn(2, parts[4])

    def test_rejects_misaligned_modality_subjects(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            ids = [f"LINE{i:02d}" for i in range(48)]
            with (folder / "response.csv").open("w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["subject_id", "cancer_class", "y"])
                for i, subject in enumerate(ids):
                    w.writerow([subject, "CNS" if i % 2 else "RENAL", float(i)/20])
            for name, width in (("rna.csv", RNA_CAP + 1), ("metabolite.csv", METAB_CAP + 1)):
                with (folder / name).open("w", newline="") as f:
                    w = csv.writer(f)
                    w.writerow(["subject_id"]+[f"feature_{k}" for k in range(width)])
                    for i, subject in enumerate(ids):
                        w.writerow([subject]+[float(i+k) for k in range(width)])
            load_dataset(folder)
            # A modality row permutation must fail, not implicitly re-pair.
            path = folder / "metabolite.csv"
            with path.open(newline="") as f:
                rows = list(csv.reader(f))
            rows[1], rows[2] = rows[2], rows[1]
            with path.open("w", newline="") as f:
                csv.writer(f).writerows(rows)
            with self.assertRaisesRegex(ValueError, "misaligned"):
                load_dataset(folder)

    def test_all_comparators_evaluate_same_48_and_report_only_aggregates(self):
        rng = np.random.default_rng(108)
        ids = tuple(f"L{i:02d}" for i in range(48))
        classes = tuple(["CNS"]*24 + ["RENAL"]*24)
        rna = rng.normal(size=(48, 320))
        metab = rng.normal(size=(48, 130))
        y = rng.normal(size=48)
        report = evaluate(ids, classes, y, rna, metab, repeats=1)
        self.assertEqual(report["n_biological_cell_lines"], 48)
        self.assertEqual(len(report["metrics"]), 16)
        self.assertEqual(report["selected_feature_count_per_outer_fold"]["rna"], [300]*5)
        self.assertNotIn("L00", str(report))
        for metrics in report["metrics"].values():
            self.assertTrue(np.isfinite(metrics["mae"]))
            self.assertTrue(np.isfinite(metrics["rmse"]))


if __name__ == "__main__":
    unittest.main()
