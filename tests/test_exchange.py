"""Synthetic-only contract tests; no participants or clinical claims."""
import csv
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from omics2pk.exchange import ExchangeError, inspect_exchange, prepare_sweet


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def fixture(root, value_type="normalized_log2_expression", repeated=False, private=False):
    root.mkdir(parents=True, exist_ok=True)
    write(root / "samples.csv", ["assay_id","sample_id","subject_id","omic","tissue","collection_time"], [
        ["A1","S1","P1","transcriptomics","liver","baseline"],
        ["A2","S2","P1" if repeated else "P2","transcriptomics","liver","baseline"]])
    write(root / "features.csv", ["omic","feature_id","id_namespace"], [
        ["transcriptomics","ENSG000001","Ensembl"], ["transcriptomics","ENSG000002","Ensembl"]])
    write(root / "matrices/transcriptomics.csv", ["feature_id","A1","A2"], [
        ["ENSG000001","1.1","2.0"], ["ENSG000002","3.4","4.0"]])
    ref = lambda p: {"path":p, "sha256":sha(root / p)}
    m = {"schema_version":"1.0.0", "analysis_id":"synthetic-fixture",
         "study":{"organism":"human","biological_unit":"participant","data_access":"private" if private else "public"},
         "files":{"samples":ref("samples.csv"),"features":ref("features.csv")},
         "matrices":[{"omic":"transcriptomics","file":ref("matrices/transcriptomics.csv"),
                      "value_type":value_type,"orientation":"feature_rows_sample_columns","missing_value":"NA"}],
         "software":{"producer":"synthetic fixture","commit":"abcdef1234567",
                     "r_version":None,"python_version":None,"dependency_lockfiles":[]},
         "validation":{"status":"exploratory","methods":[]}}
    (root / "manifest.json").write_text(json.dumps(m))
    return m


class ContractTests(unittest.TestCase):
    def test_valid_exchange_and_sweet_input(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)
            fixture(p/"exchange")
            report = inspect_exchange(p/"exchange")
            self.assertEqual(report["subjects"],2)
            self.assertFalse(report["network_inferred"])
            out = prepare_sweet(p/"exchange",p/"out")
            self.assertEqual(out["gene_count"],2)
            self.assertIn("gene\tA1\tA2", (p/"out/expression.txt").read_text())
            self.assertIn("ENSG000001\t1.1\t2.0",(p/"out/expression.txt").read_text())

    def test_private_local_not_uploaded(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)/"exchange"
            fixture(p, private=True)
            self.assertEqual(inspect_exchange(p)["data_access"],"private")

    def test_tamper_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)/"exchange"
            fixture(p)
            (p/"samples.csv").write_text("corrupt")
            with self.assertRaisesRegex(ExchangeError,"Checksum mismatch"):
                inspect_exchange(p)

    def test_traversal_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)/"exchange"
            m = fixture(p)
            m["files"]["samples"]["path"] = "../../etc/passwd"
            (p/"manifest.json").write_text(json.dumps(m))
            with self.assertRaisesRegex(ExchangeError,"Unsafe"):
                inspect_exchange(p)

    def test_raw_counts_rejected_for_sweet(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)
            fixture(p/"exchange",value_type="raw_counts")
            self.assertEqual(inspect_exchange(p/"exchange")["subjects"],2)
            with self.assertRaisesRegex(ExchangeError,"normalized/transformed"):
                prepare_sweet(p/"exchange",p/"out")

    def test_repeated_subject_rejected_for_sweet(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)
            fixture(p/"exchange",repeated=True)
            with self.assertRaisesRegex(ExchangeError,"Repeated subject"):
                prepare_sweet(p/"exchange",p/"out")

    def test_missing_value_rejected_before_writing(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)
            root=p/"exchange";m=fixture(root)
            write(root/"matrices/transcriptomics.csv",["feature_id","A1","A2"],[
                ["ENSG000001","NA","2"], ["ENSG000002","3","4"]])
            m["matrices"][0]["file"]["sha256"]=sha(root/"matrices/transcriptomics.csv")
            (root/"manifest.json").write_text(json.dumps(m))
            with self.assertRaisesRegex(ExchangeError,"cannot contain NA"):
                prepare_sweet(root,p/"out")
            self.assertFalse((p/"out").exists())


    def _rewrite_expression(self, root, rows):
        m=json.loads((root/"manifest.json").read_text())
        write(root/"matrices/transcriptomics.csv",["feature_id","A1","A2"],rows)
        m["matrices"][0]["file"]["sha256"]=sha(root/"matrices/transcriptomics.csv")
        (root/"manifest.json").write_text(json.dumps(m))

    def test_constant_gene_rejected_before_sweet_export(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td); root=p/"exchange"; fixture(root)
            self._rewrite_expression(root,[
                ["ENSG000001","1","1"],
                ["ENSG000002","3","4"],
            ])
            with self.assertRaisesRegex(ExchangeError,"constant gene"):
                prepare_sweet(root,p/"out")
            self.assertFalse((p/"out").exists())

    def test_constant_sample_rejected_before_sweet_export(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td); root=p/"exchange"; fixture(root)
            self._rewrite_expression(root,[
                ["ENSG000001","1","2"],
                ["ENSG000002","1","3"],
            ])
            with self.assertRaisesRegex(ExchangeError,"constant sample"):
                prepare_sweet(root,p/"out")

    def test_nonconstant_zero_sum_gene_rejected_for_upstream_step2_semantics(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td); root=p/"exchange"; fixture(root)
            self._rewrite_expression(root,[
                ["ENSG000001","-1","1"],
                ["ENSG000002","2","3"],
            ])
            with self.assertRaisesRegex(ExchangeError,"deletes zero-sum gene"):
                prepare_sweet(root,p/"out")

    def test_negative_values_are_not_blanket_rejected_when_upstream_safe(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td); root=p/"exchange"; fixture(root)
            self._rewrite_expression(root,[
                ["ENSG000001","-1","2"],
                ["ENSG000002","3","5"],
            ])
            out=prepare_sweet(root,p/"out")
            self.assertTrue(out["sweet_input_preflight"]["finite_expression"])
            self.assertEqual(out["sweet_input_preflight"]["constant_genes"],0)
            self.assertEqual(out["sweet_input_preflight"]["zero_sum_nonconstant_genes"],0)

    def test_assay_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)/"exchange"
            m=fixture(p)
            write(p/"matrices/transcriptomics.csv",["feature_id","A1","BAD"],[
                ["ENSG000001","1","2"],["ENSG000002","3","4"]])
            m["matrices"][0]["file"]["sha256"]=sha(p/"matrices/transcriptomics.csv")
            (p/"manifest.json").write_text(json.dumps(m))
            with self.assertRaisesRegex(ExchangeError,"Sample/assay mismatch"):
                inspect_exchange(p)


if __name__ == "__main__":
    unittest.main()
