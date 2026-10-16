"""Synthetic SWEET output parsing; does not replicate upstream SWEET inference."""
import csv
import tempfile
import unittest
from pathlib import Path

from omics2pk.sweet_features import SWEETOutputError, summarize_networks


def write_edges(root, name, rows):
    root.mkdir(parents=True, exist_ok=True)
    with (root / f"{name}.txt").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["gene1", "gene2", "raw_edge_score"])
        w.writerows(rows)


class SWEETFeatureTests(unittest.TestCase):
    def _case(self, p):
        edges = p / "edges"
        write_edges(edges, "P1", [("A", "B", "0.2"), ("A", "C", "-0.4"), ("B", "C", "0.1")])
        write_edges(edges, "P2", [("A", "B", "-0.1"), ("A", "C", "0.5"), ("B", "C", "0.2")])
        (p / "patients.txt").write_text("P1\nP2\n")
        return edges, p / "patients.txt"

    def test_correct_unbiased_feature_calculation(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)
            edges, patients = self._case(p)
            result = summarize_networks(edges, patients, p / "out")
            self.assertEqual((result["n_patients"], result["n_genes"]), (2, 3))
            with (p / "out" / "gene_absolute_strength.csv").open() as f:
                records = list(csv.DictReader(f))
            self.assertAlmostEqual(float(records[0]["abs_strength::A"]), 0.3)
            self.assertAlmostEqual(float(records[0]["abs_strength::B"]), 0.15)
            self.assertAlmostEqual(float(records[0]["abs_strength::C"]), 0.25)
            with (p / "out" / "network_summary.csv").open() as f:
                table = list(csv.DictReader(f))
            self.assertAlmostEqual(float(table[0]["edge_mean"]), -0.1 / 3)
            self.assertAlmostEqual(float(table[0]["negative_edge_fraction"]), 1 / 3)
            self.assertFalse(result["out_of_sample_sweet_validated"])

    def test_missing_edges_blocked(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)
            edges, patients = self._case(p)
            write_edges(edges, "P2", [("A", "B", ".2"), ("B", "C", ".3")])
            with self.assertRaisesRegex(SWEETOutputError, "Incomplete"):
                summarize_networks(edges, patients, p / "out")
            self.assertFalse((p / "out").exists())

    def test_duplicate_edges_blocked(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)
            edges, patients = self._case(p)
            write_edges(edges, "P2", [("A", "B", ".2"), ("B", "A", ".3"), ("B", "C", ".3")])
            with self.assertRaisesRegex(SWEETOutputError, "Repeated edge"):
                summarize_networks(edges, patients, p / "out")

    def test_nan_blocked(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)
            edges, patients = self._case(p)
            write_edges(edges, "P2", [("A", "B", "nan"), ("A", "C", ".3"), ("B", "C", ".3")])
            with self.assertRaisesRegex(SWEETOutputError, "Non-finite"):
                summarize_networks(edges, patients, p / "out")

    def test_filename_traversal_blocked(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)
            edges, patients = self._case(p)
            patients.write_text("../P1\nP2\n")
            with self.assertRaisesRegex(SWEETOutputError, "unsafe"):
                summarize_networks(edges, patients, p / "out")

    def test_unmatched_gene_universe_blocked(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)
            edges, patients = self._case(p)
            write_edges(edges, "P2", [("A", "B", ".2"), ("A", "D", ".3"), ("B", "D", ".3")])
            with self.assertRaisesRegex(SWEETOutputError, "universes"):
                summarize_networks(edges, patients, p / "out")

    def test_output_can_not_destroy_edge_inputs(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)
            edges, patients = self._case(p)
            with self.assertRaisesRegex(SWEETOutputError, "overwrite"):
                summarize_networks(edges, patients, edges / "results")


if __name__ == "__main__":
    unittest.main()
