"""Synthetic MultiOmics V1 > same-person clinical PK sidecar, no patient data."""
import csv
from hashlib import sha256
import json
from pathlib import Path
import tempfile
import unittest

try:
    import numpy as np
    from omics2pk.multimodal_exchange_pk import (
        run_multimodal_exchange_pk, MultimodalExchangePKError)
    from omics2pk.multimodal_pk_ablation import MultimodalPKBenchmarkError
except ModuleNotFoundError as exc:
    if exc.name not in ("numpy","sklearn"):
        raise
    np=None


def csvfile(path, header, rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f);w.writerow(header);w.writerows(rows)


def checksum(path):
    return sha256(path.read_bytes()).hexdigest()


@unittest.skipIf(np is None,"Optional source PK estimation dependencies not installed")
class MultimodalOriginalExchangeFictionalTests(unittest.TestCase):
    def setUp(self):
        self.td=tempfile.TemporaryDirectory()
        self.addCleanup(self.td.cleanup)
        self.root=Path(self.td.name)
        rng=np.random.default_rng(20261019)
        n=14
        self.omics=("transcriptomics","proteomics")
        subjects=[f"FAKE_{i:03d}" for i in range(n)]
        samples=[]
        features=[]
        matrices=[]
        self.sidecar=[]
        for i,p in enumerate(subjects):
            assays={m:f"FICTIONAL_{m}_ASSAY_{i}" for m in self.omics}
            self.sidecar.append(dict(
                subject_id=p,assay_ids_by_omic=assays,
                sample_ids_by_omic={m:f"FAKE_SPECIMEN_{i}" for m in self.omics},
                collected_at_by_omic={m:"2025-01-02T08:00:00+00:00"
                                      for m in self.omics},
                index_dose_at="2025-01-03T09:00:00Z",
                clinical_values=rng.normal(size=2).tolist(),
                genotype_dosages=rng.integers(0,3,size=2).tolist(),
                measured_pk=float(np.exp(rng.normal(2.1,0.25))),
            ))
            for m in self.omics:
                samples.append([assays[m],f"FAKE_SPECIMEN_{i}",p,m,"blood",
                                "2025-01-02T08:00:00+00:00"])
        csvfile(self.root/"samples.csv",
                ["assay_id","sample_id","subject_id","omic","tissue",
                 "collection_time"],samples)
        for m in self.omics:
            genes=[f"FAKE_{m}_GENE_{j}" for j in range(4)]
            features += [[m,g,"fictional_original"] for g in genes]
            assay_list=[d["assay_ids_by_omic"][m] for d in self.sidecar]
            arr=rng.normal(5,1,(4,n))
            path=f"matrices/{m}.csv"
            csvfile(self.root/path,["feature_id"]+assay_list,
                    [[g,*arr[j].tolist()] for j,g in enumerate(genes)])
            matrices.append(dict(omic=m,
                     file=dict(path=path,sha256=checksum(self.root/path)),
                     value_type="normalized_expression" if m=="transcriptomics"
                                else "normalized_protein_abundance",
                     orientation="feature_rows_sample_columns"))
        csvfile(self.root/"features.csv",
                ["omic","feature_id","id_namespace"],features)
        self.manifest=dict(
            schema_version="1.0.0",analysis_id="fake_no_actual_people",
            study={"biological_unit":"participant","data_access":"private"},
            files={
                "samples":{"path":"samples.csv","sha256":checksum(self.root/"samples.csv")},
                "features":{"path":"features.csv","sha256":checksum(self.root/"features.csv")},
            },
            matrices=matrices,
            software={"producer":"fictional","commit":"fictional_sha",
                      "dependency_lockfiles":[]},
            validation={"status":"exploratory"},
        )
        self.save()
        self.kwargs=dict(selected_omics=self.omics,
                         pk_target_name="FICTIONAL_PK",
                         pk_target_unit="FICTIONAL_UNITS",
                         source_kind="synthetic",source_sha256="b"*64)

    def save(self):
        (self.root/"manifest.json").write_text(json.dumps(self.manifest),
                                                encoding="utf-8")

    def call(self):
        return run_multimodal_exchange_pk(self.root,self.sidecar,**self.kwargs)

    def test_real_native_manifest_and_patient_level_link(self):
        r=self.call()
        self.assertEqual(r["n_independent_people"],14)
        self.assertEqual(r["n_complete_case_people_used"],14)
        self.assertIn("early_fusion_pca2",r["models"])
        self.assertIn("late_fusion_fixed_average",r["models"])
        self.assertTrue(r["all_test_arms_share_identical_complete_case_people"])
        self.assertNotIn("FAKE_",json.dumps(r))

    def test_private_original_site_labels_drive_whole_center_holdout(self):
        for i,r in enumerate(self.sidecar):
            r["site_id"]=f"FAKE_CENTER_{i%4}"
        self.kwargs["validation_mode"]="leave_one_site_out"
        out=self.call()
        self.assertEqual(out["n_outer_splits"],4)
        self.assertEqual(out["validation_mode"],"leave_one_site_out")
        self.assertTrue(out["site_disjoint_holdout_is_not_independent_external_cohort"])
        self.assertNotIn("FAKE_CENTER_",json.dumps(out))
        self.sidecar.reverse()
        self.assertEqual(out,self.call())

    def test_site_metadata_cannot_be_silently_mixed(self):
        for i,r in enumerate(self.sidecar):
            r["site_id"]=f"FAKE_CENTER_{i%4}"
        with self.assertRaisesRegex(MultimodalExchangePKError,
                                    "do not silently use mixed-site folds"):
            self.call()

    def test_partial_center_identity_rejected(self):
        self.sidecar[0]["site_id"]="FAKE_CENTER_A"
        self.kwargs["validation_mode"]="leave_one_site_out"
        with self.assertRaisesRegex(MultimodalExchangePKError,
                                    "Every PK person needs"):
            self.call()

    def test_permuted_original_clinical_rows_gives_same_report(self):
        report=self.call()
        self.sidecar.reverse()
        self.assertEqual(report,self.call())

    def test_wrong_assay_id_must_not_join_by_position(self):
        a,b=self.sidecar[0],self.sidecar[1]
        m="transcriptomics"
        a["assay_ids_by_omic"][m],b["assay_ids_by_omic"][m]=(
            b["assay_ids_by_omic"][m],a["assay_ids_by_omic"][m])
        with self.assertRaisesRegex(MultimodalExchangePKError,"identity mismatch"):
            self.call()

    def test_crossed_source_specimen_must_fail_even_if_assay_correct(self):
        m="transcriptomics"
        a,b=self.sidecar[0],self.sidecar[1]
        a["sample_ids_by_omic"][m],b["sample_ids_by_omic"][m]=(
            b["sample_ids_by_omic"][m],a["sample_ids_by_omic"][m])
        with self.assertRaisesRegex(MultimodalExchangePKError,
                                    "specimen-to-patient"):
            self.call()

    def test_missing_source_specimen_mapping_cannot_pass(self):
        self.sidecar[0].pop("sample_ids_by_omic")
        with self.assertRaisesRegex(MultimodalExchangePKError,
                                    "specimen/time maps"):
            self.call()

    def test_qualified_original_requires_absolute_source_clock(self):
        self.kwargs["source_kind"]="original_source_declared"
        self.manifest["validation"]["status"]="qualified_for_declared_scope"
        with (self.root/"samples.csv").open("r",newline="") as f:
            records=list(csv.DictReader(f))
        for rec in records:
            rec["collection_time"]="T0"
        csvfile(self.root/"samples.csv",
                ["assay_id","sample_id","subject_id","omic","tissue",
                 "collection_time"],
                [[r[k] for k in ("assay_id","sample_id","subject_id","omic",
                                 "tissue","collection_time")] for r in records])
        self.manifest["files"]["samples"]["sha256"]=checksum(self.root/"samples.csv")
        self.save()
        with self.assertRaisesRegex(MultimodalExchangePKError,
                                    "absolute specimen collection clock"):
            self.call()

    def test_t0_is_not_an_absolute_timestamp(self):
        self.sidecar[0]["collected_at_by_omic"]["transcriptomics"]="T0"
        with self.assertRaises((MultimodalExchangePKError,MultimodalPKBenchmarkError)):
            self.call()

    def test_missing_modality_requires_explicit_complete_case(self):
        self.manifest["matrices"]=[x for x in self.manifest["matrices"]]
        # Study can omit clinical evidence for one modality in one subject,
        # but should explicitly request complete-case restriction.
        self.sidecar[0]["assay_ids_by_omic"].pop("proteomics")
        with self.assertRaisesRegex(MultimodalExchangePKError,"identity mismatch"):
            self.call()

    def test_original_qualified_status_never_silently_assumed(self):
        self.kwargs["source_kind"]="original_source_declared"
        with self.assertRaisesRegex(MultimodalExchangePKError,"upstream source qualification"):
            self.call()

    def test_modified_omics_matrix_sha_fails(self):
        path=self.root/"matrices/proteomics.csv"
        with path.open("a") as f:
            f.write("FAKE_ADDENDUM")
        from omics2pk.exchange import ExchangeError
        with self.assertRaisesRegex(ExchangeError,"Checksum"):
            self.call()


if __name__=="__main__":
    unittest.main()
