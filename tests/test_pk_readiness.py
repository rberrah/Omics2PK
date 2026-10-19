"""Entirely synthetic multimodal PK readiness test cohort. No true people."""
import csv
from contextlib import redirect_stdout
from hashlib import sha256
from io import StringIO
import json
from pathlib import Path
import tempfile
import unittest

from omics2pk.pk_readiness import assess_pk_readiness, PKReadinessError
from omics2pk.cli import main


def _save(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        out=csv.writer(f)
        out.writerow(header)
        out.writerows(rows)


def _digest(path):
    return sha256(path.read_bytes()).hexdigest()


class SyntheticOriginalPKSourceReadinessTests(unittest.TestCase):
    def setUp(self):
        self.temporary=tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root=Path(self.temporary.name)
        self.ids=[f"FICTIONAL_PERSON_{i:02d}" for i in range(12)]
        self.omics=("transcriptomics","proteomics")
        self.assay={}
        self.specimen={}
        samples=[]
        features=[]
        matrices=[]
        for omic in self.omics:
            assays=[]
            for i, pid in enumerate(self.ids):
                assay=f"FICTIONAL_{omic}_ASSAY_{i:02d}"
                sample=f"FICTIONAL_{omic}_SPECIMEN_{i:02d}"
                self.assay[(omic,pid)]=assay
                self.specimen[(omic,pid)]=sample
                assays.append(assay)
                samples.append([assay,sample,pid,omic,"blood_plasma",
                                "2025-01-01T09:00:00+00:00"])
            features.extend([[omic,f"FAKE_{omic}_{k}","fictional"]
                             for k in range(3)])
            name=f"matrices/{omic}.csv"
            _save(self.root/name,["feature_id"]+assays,
                  [[f"FAKE_{omic}_{k}"] + [i+k for i in range(12)]
                   for k in range(3)])
            matrices.append({
                "omic":omic,
                "file":{"path":name,"sha256":_digest(self.root/name)},
                "value_type":"normalized_expression",
                "orientation":"feature_rows_sample_columns"})
        _save(self.root/"samples.csv",
              ["assay_id","sample_id","subject_id","omic","tissue",
               "collection_time"],samples)
        _save(self.root/"features.csv",
              ["omic","feature_id","id_namespace"],features)
        self.manifest={
            "schema_version":"1.0.0","analysis_id":"FICTIONAL_ONLY",
            "study":{"biological_unit":"participant","data_access":"private"},
            "software":{"producer":"FAKE","commit":"FAKE_SHA","dependency_lockfiles":[]},
            "validation":{"status":"exploratory"},
            "files":{
                "samples":{"path":"samples.csv","sha256":_digest(self.root/"samples.csv")},
                "features":{"path":"features.csv","sha256":_digest(self.root/"features.csv")},
            },
            "matrices":matrices,
        }
        (self.root/"manifest.json").write_text(json.dumps(self.manifest))
        self.availability=[]
        for pid in self.ids:
            self.availability.append({
                "subject_id":pid,"has_clinical":True,"has_measured_pk":True,
                "has_pgx":True,"has_dose_timing":True,
                "assay_ids_by_omic":{o:self.assay[(o,pid)] for o in self.omics},
                "sample_ids_by_omic":{o:self.specimen[(o,pid)] for o in self.omics},
                "collected_at_by_omic":{o:"2025-01-01T09:00:00Z"
                                        for o in self.omics},
                "index_dose_at":"2025-01-02T09:00:00Z",
            })

    def test_source_only_cli_no_access_to_patient_pk_or_pgx(self):
        result=assess_pk_readiness(self.root)
        self.assertEqual(result["n_original_biological_people"],12)
        self.assertEqual(result["n_omic_assay_rows"],24)
        self.assertEqual(result["n_intersecting_biological_people_by_modalities"]
                         ["proteomics+transcriptomics"],12)
        self.assertFalse(result["availability_sidecar_supplied"])
        self.assertIsNone(result["n_predose_qualified_by_modalities"])
        out=StringIO()
        with redirect_stdout(out):
            main(["audit-pk-readiness",str(self.root)])
        reported=json.loads(out.getvalue())
        self.assertEqual(reported,result)
        self.assertNotIn("FICTIONAL_PERSON",out.getvalue())

    def test_full_pk_pgx_source_and_predose_exact_join(self):
        report=assess_pk_readiness(self.root,self.availability)
        self.assertEqual(report["n_clinical_pk_pgx_reference_declared"],12)
        self.assertEqual(report["n_reference_and_verified_original_specimen_join_by_modality"]
                         ["transcriptomics"],12)
        self.assertEqual(report["n_predose_qualified_by_modalities"]
                         ["proteomics+transcriptomics"],12)
        self.assertEqual(report["blocking_requirements"],[])
        txt=json.dumps(report)
        self.assertNotIn("FICTIONAL_PERSON",txt)
        self.assertNotIn("FICTIONAL_",txt)

    def test_molecular_specimen_crosswalk_mismatch_abstains_for_one(self):
        r=dict(self.availability[0])
        r["sample_ids_by_omic"]=dict(r["sample_ids_by_omic"])
        r["sample_ids_by_omic"]["proteomics"]=self.specimen[(
            "proteomics",self.ids[1])]
        rows=[r]+self.availability[1:]
        result=assess_pk_readiness(self.root,rows)
        self.assertEqual(result["n_clinical_pk_pgx_reference_declared"],12)
        self.assertEqual(result["n_predose_qualified_by_modalities"]
                         ["proteomics+transcriptomics"],11)

    def test_missing_pgx_flag_excludes_model_reference_not_people(self):
        self.availability[0]["has_pgx"]=False
        result=assess_pk_readiness(self.root,self.availability)
        self.assertEqual(result["n_clinical_pk_without_full_pgx_declared"],1)
        self.assertEqual(result["n_clinical_pk_pgx_reference_declared"],11)
        self.assertEqual(result["n_predose_qualified_by_modalities"]
                         ["proteomics+transcriptomics"],11)

    def test_fake_t0_not_accepted_as_source_timestamp(self):
        self.availability[0]["collected_at_by_omic"]["transcriptomics"]="T0"
        result=assess_pk_readiness(self.root,self.availability)
        self.assertEqual(result["n_predose_qualified_by_modalities"]
                         ["proteomics+transcriptomics"],11)

    def test_future_collected_omics_does_not_qualify(self):
        self.availability[0]["index_dose_at"]="2024-12-31T09:00:00Z"
        result=assess_pk_readiness(self.root,self.availability)
        self.assertEqual(result["n_predose_qualified_by_modalities"]
                         ["proteomics+transcriptomics"],11)

    def test_source_without_original_absolute_collection_clock_never_qualifies(self):
        rows=list(csv.DictReader((self.root/"samples.csv").open(
            encoding="utf-8",newline="")))
        for r in rows:
            r["collection_time"]=""
        _save(self.root/"samples.csv",list(rows[0]),[[r[k] for k in rows[0]]
                                                     for r in rows])
        self.manifest["files"]["samples"]["sha256"]=_digest(
            self.root/"samples.csv")
        (self.root/"manifest.json").write_text(json.dumps(self.manifest))
        report=assess_pk_readiness(self.root,self.availability)
        self.assertEqual(report["n_clinical_pk_pgx_reference_declared"],12)
        self.assertEqual(report["n_predose_qualified_by_modalities"]
                         ["proteomics+transcriptomics"],0)
        self.assertTrue(report["blocking_requirements"])

    def test_missing_one_proteome_specimen_decreases_true_complete_case_n(self):
        rows=list(csv.DictReader((self.root/"samples.csv").open(
            encoding="utf-8",newline="")))
        target=self.assay[("proteomics",self.ids[-1])]
        rows=[r for r in rows if r["assay_id"]!=target]
        _save(self.root/"samples.csv",list(rows[0]),[[r[k] for k in rows[0]]
                                                     for r in rows])
        path=self.root/"matrices/proteomics.csv"
        with path.open(encoding="utf-8",newline="") as f:
            records=list(csv.DictReader(f))
        columns=[k for k in records[0] if k!=target]
        _save(path,columns,[[r[k] for k in columns] for r in records])
        self.manifest["files"]["samples"]["sha256"]=_digest(self.root/"samples.csv")
        for item in self.manifest["matrices"]:
            if item["omic"]=="proteomics":
                item["file"]["sha256"]=_digest(path)
        (self.root/"manifest.json").write_text(json.dumps(self.manifest))
        report=assess_pk_readiness(self.root,self.availability)
        self.assertEqual(report["n_original_biological_people"],12)
        self.assertEqual(report["n_omic_assay_rows"],23)
        self.assertEqual(report["n_unique_one_assay_people_by_modality"]
                         ["proteomics"],11)
        self.assertEqual(report["n_intersecting_biological_people_by_modalities"]
                         ["proteomics+transcriptomics"],11)
        self.assertEqual(report["n_predose_qualified_by_modalities"]
                         ["proteomics+transcriptomics"],11)

    def test_unknown_biological_person_is_rejected(self):
        self.availability[0]["subject_id"]="FICTIONAL_UNKNOWN"
        with self.assertRaisesRegex(PKReadinessError,"unknown participants"):
            assess_pk_readiness(self.root,self.availability)

    def test_missing_boolean_cannot_be_inferred(self):
        self.availability[0]["has_pgx"]="yes"
        with self.assertRaisesRegex(PKReadinessError,"explicit booleans"):
            assess_pk_readiness(self.root,self.availability)

    def test_tampered_source_sha_is_rejected(self):
        with (self.root/"samples.csv").open("a") as f:
            f.write("FICTIONAL_TAMPER")
        from omics2pk.exchange import ExchangeError
        with self.assertRaisesRegex(ExchangeError,"Checksum"):
            assess_pk_readiness(self.root,self.availability)


if __name__=="__main__":
    unittest.main()
