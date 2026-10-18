"""Synthetic-only source/PK benchmark tests. No public original patient data."""
import json
import tempfile
import unittest
from pathlib import Path
from zipfile import ZipFile

try:
    import numpy as np
    import sklearn
    from scripts.validate_local_jarmusch14_gnps_biom_pk import (
        GNPSOriginalSourceError, leave_one_person_out, load_original_source,
        summarize_models)
except ModuleNotFoundError as exc:
    if exc.name not in ("numpy", "sklearn"):
        raise
    np = None


@unittest.skipIf(np is None, "Requires optional GNPS benchmark numpy and sklearn")
class JarmuschGNPSOriginalTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(20261019)
        self.x = np.exp(rng.normal(4, 1.2, size=(13, 21)))
        self.x[rng.random(self.x.shape) < 0.15] = 0
        self.c = np.column_stack((rng.uniform(50, 100, size=13),
                                  rng.integers(0, 2, size=13))).astype(float)
        self.g = np.array([["A/A" if i % 2 else "A/G",
                            "G/G" if i % 3 else "G/A"] for i in range(13)],
                          dtype=object)
        self.y = np.exp(rng.normal(2.8, 0.3, size=13))

    def test_out_of_fold_models_share_exact_individuals(self):
        result = summarize_models(self.x, self.c, self.g, self.y, shuffles=2)
        self.assertEqual(result["n_unique_original_volunteers"], 13)
        self.assertEqual(len(result["models"]), 3)
        self.assertEqual(result["negative_control_n_shuffled_subject_mappings"], 2)
        self.assertTrue(all(np.isfinite(x["mae"]) for x in result["models"].values()))
        self.assertEqual(result, summarize_models(self.x, self.c, self.g, self.y, shuffles=2))
        serialized = json.dumps(result)
        self.assertNotIn("A/A", serialized)
        self.assertNotIn(str(self.y[0]), serialized)
        self.assertTrue(result["independent_external_validation"] is False)

    def test_followup_12_persons_does_not_claim_13(self):
        result = summarize_models(self.x[:12], self.c[:12], self.g[:12],
                                  self.y[:12], shuffles=0)
        self.assertEqual(result["n_unique_original_volunteers"], 12)

    def test_nonfinite_or_nonpositive_observations_rejected(self):
        y = self.y.copy()
        y[0] = 0
        with self.assertRaisesRegex(GNPSOriginalSourceError, "Matched"):
            leave_one_person_out(self.x, self.c, self.g, y)
        x = self.x.copy()
        x[1, 4] = np.nan
        with self.assertRaisesRegex(GNPSOriginalSourceError, "Matched"):
            leave_one_person_out(x, self.c, self.g, self.y)

    def test_fail_closed_on_non_original_archive(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "not_original_source.zip"
            with ZipFile(p, "w") as z:
                z.writestr("private_fake_original.txt", "invented")
            with self.assertRaisesRegex(GNPSOriginalSourceError, "intact"):
                load_original_source(p)


if __name__ == "__main__":
    unittest.main()
