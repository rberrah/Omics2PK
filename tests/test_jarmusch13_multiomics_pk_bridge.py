"""Fictitious subjects and ions only. No real GNPS/MultiOmics data in CI."""
import csv
from hashlib import sha256
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from omics2pk.exchange import ExchangeError
from scripts.validate_private_jarmusch13_multiomics_pipeline import (
    JarmuschMS1PipelineError, check_real_native_multiomics_export,
    read_prepared_private_zip
)


def write(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


class SyntheticNativeMS1Tests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name)
        self.samples = [f"FAKE_PERSON_{i}" for i in range(13)]
        self.assays = [f"FAKE_MS_{i}" for i in range(13)]
        self.features = ["GNPS_MS1:FAKE_1", "GNPS_MS1:FAKE_2", "GNPS_MS1:FAKE_3"]
        self.values = np.arange(39,dtype=float).reshape(13,3)
        self.source = {
            "subject_ids":tuple(self.samples),
            "assay_ids":tuple(self.assays),
            "feature_ids":tuple(self.features),
            "matrix":self.values,
        }
        write(self.path/"samples.csv",
              ["assay_id","sample_id","subject_id","omic","tissue","timepoint"],
              [[self.assays[i],f"FAKE_SPECIMEN_{i}",self.samples[i],
                "metabolomics","blood_plasma","T0"] for i in range(13)])
        write(self.path/"features.csv",
              ["omic","feature_id","id_namespace"],
              [["metabolomics",feature,"GNPS_MS1_unknown_chemical"] for feature in self.features])
        self.write_matrix()
        def record(name):
            return {"path":name,"sha256":digest(self.path/name)}
        self.manifest = {
            "schema_version":"1.0.0",
            "analysis_id":"fake-out-of-sample-test",
            "study":{"biological_unit":"participant","data_access":"private"},
            "software":{"commit":"fake_test_commit","producer":"fictional_fixture",
                        "dependency_lockfiles":[]},
            "validation":{"status":"exploratory"},
            "files":{"samples":record("samples.csv"),
                     "features":record("features.csv")},
            "matrices":[{"omic":"metabolomics",
                         "file":record("matrices/metabolomics.csv"),
                         "value_type":"ms1_peak_area_original",
                         "orientation":"feature_rows_sample_columns"}],
        }
        self.save_manifest()

    def write_matrix(self):
        write(self.path/"matrices/metabolomics.csv", ["feature_id"] + self.assays,
              [[f]+list(self.values[:,i]) for i,f in enumerate(self.features)])

    def save_manifest(self):
        (self.path/"manifest.json").write_text(json.dumps(self.manifest))

    def run_gate(self):
        with patch("scripts.validate_private_jarmusch13_multiomics_pipeline.N_FEATURES",3):
            return check_real_native_multiomics_export(self.path,self.source)

    def test_real_native_schema_and_numeric_alignment(self):
        r=self.run_gate()
        self.assertEqual(r["n_feature_assay_numeric_matches"],39)
        self.assertEqual(r["status"],
                         "ACTUAL_NATIVE_MULTIOMICS_EXPORT_SOURCE_ID_AND_NUMERIC_PARITY_PASS")
        serialized=json.dumps(r)
        self.assertNotIn("FAKE_PERSON",serialized)
        self.assertNotIn("GNPS_MS1:FAKE",serialized)

    def test_crossed_subject_assay_not_accepted(self):
        self.source["assay_ids"]=tuple(reversed(self.assays))
        with self.assertRaisesRegex(JarmuschMS1PipelineError,"join"):
            self.run_gate()

    def test_modified_source_hash_fails_first(self):
        with (self.path/"matrices/metabolomics.csv").open("a") as f:
            f.write("CORRUPT,1,2,3\\n")
        with self.assertRaisesRegex(ExchangeError,"Checksum"):
            self.run_gate()

    def test_recomputed_sha_cannot_hide_modified_numeric_value(self):
        self.values[0,0] += 10
        self.write_matrix()
        self.manifest["matrices"][0]["file"]["sha256"]=digest(
            self.path/"matrices/metabolomics.csv")
        self.save_manifest()
        with self.assertRaisesRegex(JarmuschMS1PipelineError,"changed MS1 values"):
            self.run_gate()

    def test_original_source_must_be_present_and_sha_pinned(self):
        with self.assertRaisesRegex(JarmuschMS1PipelineError,"absent"):
            read_prepared_private_zip(self.path/"NOT_ATTACHED_PRIVATE.zip")


if __name__=="__main__":
    unittest.main()
