"""Synthetic native MultiOmics ZIP/file contract -> private PK matched bridge."""
import csv
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

try:
    import numpy as np
    import sklearn
    from omics2pk.matched_multiomics_pk import (
        run_matched_rna_exchange_pk, MatchedExchangePKError)
    from omics2pk.predose_representation_benchmark import RepresentationBenchmarkError
except ModuleNotFoundError as exc:
    if exc.name not in ("numpy", "sklearn"):
        raise
    np = None


def write_csv(path, header, records):
    path.parent.mkdir(exist_ok=True, parents=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(records)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@unittest.skipIf(np is None, "Requires optional benchmark dependencies")
class NativeClinicalBridgeTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        rng = np.random.default_rng(20261019)
        n = 12
        subjects = [f"FAKE_PERSON_{j}" for j in range(n)]
        assays = [f"FAKE_ASSAY_{j}" for j in range(n)]
        write_csv(self.root/"samples.csv",
                  ["assay_id", "sample_id", "subject_id", "omic", "tissue"],
                  [[assays[j], f"FAKE_SPECIMEN_{j}", subjects[j], "transcriptomics", "liver"]
                   for j in range(n)])
        genes = [f"ENSG_FICT_{i}" for i in range(7)]
        write_csv(self.root/"features.csv", ["omic","feature_id","id_namespace"],
                  [["transcriptomics",g,"Ensembl"] for g in genes])
        self.expression = rng.normal(5, 1, (7, n))
        write_csv(self.root/"matrices/transcriptomics.csv", ["feature_id"] + assays,
                  [[g]+list(self.expression[i]) for i,g in enumerate(genes)])
        ref = lambda s: {"path":s,"sha256":digest(self.root/s)}
        manifest = {
            "schema_version":"1.0.0","analysis_id":"fictional_training_fixture",
            "study":{"biological_unit":"participant","data_access":"private"},
            "files":{"samples":ref("samples.csv"),"features":ref("features.csv")},
            "matrices":[{"omic":"transcriptomics",
                         "file":ref("matrices/transcriptomics.csv"),
                         "value_type":"normalized_log2_expression",
                         "orientation":"feature_rows_sample_columns"}],
            "software":{"producer":"fictional fixture","commit":"fake_testing_code",
                        "dependency_lockfiles":[]},
            "validation":{"status":"exploratory"},
        }
        (self.root/"manifest.json").write_text(json.dumps(manifest))
        self.sidecar = [
            dict(subject_id=subjects[i], assay_id=assays[i],
                 collected_at="2025-01-02T09:00:00Z",
                 first_dose_at="2025-01-03T09:00:00Z",
                 measured_pk=float(np.exp(rng.normal(2, 0.2))),
                 clinical_values=[float(rng.normal()),float(rng.normal())],
                 genotype_dosages=[int(rng.integers(0,3))])
            for i in range(n)
        ]
        self.kwargs=dict(target_name="synthetic_cl",target_unit="fabricated_L_per_h",
                         source_kind="synthetic",source_sha256="a"*64)

    def call(self):
        return run_matched_rna_exchange_pk(self.root,self.sidecar,**self.kwargs)

    def test_original_ids_joined_not_by_row_order(self):
        first=self.call()
        self.sidecar.reverse()
        second=self.call()
        self.assertEqual(first,second)
        self.assertEqual(first["n_independent_subjects"],12)
        self.assertEqual(first["primary_reference"],"clinical_pgx")
        self.assertNotIn("FAKE_PERSON",json.dumps(first))
        self.assertNotIn("FAKE_ASSAY",json.dumps(first))

    def test_subject_crosswalk_swap_fails_closed(self):
        self.sidecar[0]["assay_id"],self.sidecar[1]["assay_id"] = (
            self.sidecar[1]["assay_id"],self.sidecar[0]["assay_id"])
        with self.assertRaisesRegex(MatchedExchangePKError,"Crossed"):
            self.call()

    def test_repeated_biological_subject_abstains(self):
        self.sidecar[1]["subject_id"]=self.sidecar[0]["subject_id"]
        with self.assertRaisesRegex(MatchedExchangePKError,"Duplicate"):
            self.call()

    def test_bad_source_sha_abstains(self):
        with (self.root/"samples.csv").open("a") as f:f.write("corrupt,record\n")
        from omics2pk.exchange import ExchangeError
        with self.assertRaisesRegex(ExchangeError,"Checksum"):
            self.call()

    def test_unqualified_original_source_not_claimed(self):
        self.kwargs["source_kind"]="original_source_declared"
        with self.assertRaisesRegex(MatchedExchangePKError,"not qualified"):
            self.call()

    def test_postdose_sidecar_rejected(self):
        self.sidecar[0]["collected_at"]="2025-01-04T09:00:00Z"
        with self.assertRaisesRegex(RepresentationBenchmarkError,"Postdose"):
            self.call()


if __name__ == "__main__":
    unittest.main()
