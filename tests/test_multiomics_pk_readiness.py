"""Only explicitly synthetic MultiOmics V1 research samples and PK targets."""
import csv
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from omics2pk.exchange import ExchangeError
from omics2pk.multiomics_pk_readiness import inspect_pk_readiness
from omics2pk.nca_targets import (
    NCATargetError, qualify_nca_targets, qualify_baseline_omic_nca_join
)

SYNTH_SHA = "a" * 64


def write_csv(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w",encoding="utf-8",newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def fixture(root, *, time_columns=True, multiple_modalities=False):
    root.mkdir(parents=True,exist_ok=True)
    h=["assay_id","sample_id","subject_id","omic","tissue",
       "collection_time","collection_time_unit"] if time_columns else (
       ["assay_id","sample_id","subject_id","omic","tissue",
        "collection_time","collection_time_unit"] )
    values=[
        ["R1","S0","PAT-001","transcriptomics","blood","0","h" if time_columns else ""],
        ["R2","S1","PAT-001","transcriptomics","blood","14","h" if time_columns else ""],
        ["R3","S2","PAT-002","transcriptomics","blood","0","h" if time_columns else ""],
    ]
    if not time_columns:
        values=[r[:5]+["",""] for r in values]
    if multiple_modalities:
        values.append(["M1","S0","PAT-001","metabolomics","blood","0","h" if time_columns else ""])
    write_csv(root/"samples.csv",h,values)
    features=[["transcriptomics","ENSG000001","Ensembl"]]
    if multiple_modalities:features += [["metabolomics","CHEBI:12","ChEBI"]]
    write_csv(root/"features.csv",["omic","feature_id","id_namespace"],features)
    write_csv(root/"matrices/transcriptomics.csv",
              ["feature_id","R1","R2","R3"],[["ENSG000001",4,10,6]])
    if multiple_modalities:
        write_csv(root/"matrices/metabolomics.csv",["feature_id","M1"],[["CHEBI:12",1.2]])
    def entry(p):
        return {"path":p,"sha256":hashlib.sha256((root/p).read_bytes()).hexdigest()}
    matrices=[{"omic":t,"file":entry("matrices/"+t+".csv"),"value_type":"normalized_log2_expression" if t=="transcriptomics" else "log_intensity","orientation":"feature_rows_sample_columns","missing_value":"NA"} for t in (["transcriptomics","metabolomics"] if multiple_modalities else ["transcriptomics"])]
    m={"schema_version":"1.0.0","analysis_id":"synthetic_only",
       "study":{"organism":"human","biological_unit":"participant","data_access":"public"},
       "files":{"samples":entry("samples.csv"),"features":entry("features.csv")},
       "matrices":matrices,
       "software":{"producer":"synthetic","commit":"abcdef01234567","r_version":None,"python_version":None,"dependency_lockfiles":[]},
       "validation":{"status":"exploratory","methods":[]}}
    (root/"manifest.json").write_text(json.dumps(m))
    return root


def targets():
    return [
      {"subject_id":"PAT-001","occasion_id":"O1","dose_id":"D1","drug":"tacrolimus",
       "analyte":"tacrolimus","specimen":"whole_blood",
       "dose_time_h":12,"window_end_h":24,"endpoint":"Cmax","value":11.2,
       "unit":"mg/L","n_observed_concentrations":6,"pk_source_sha256":SYNTH_SHA},
      {"subject_id":"PAT-001","occasion_id":"O1","dose_id":"D1","drug":"tacrolimus",
       "analyte":"tacrolimus","specimen":"whole_blood",
       "dose_time_h":12,"window_end_h":24,"endpoint":"Tmax","value":1.2,
       "unit":"h","n_observed_concentrations":6,"pk_source_sha256":SYNTH_SHA},
      {"subject_id":"PAT-001","occasion_id":"O1","dose_id":"D1","drug":"tacrolimus",
       "analyte":"tacrolimus","specimen":"whole_blood",
       "dose_time_h":12,"window_end_h":24,"endpoint":"AUC0last","value":55.3,
       "unit":"mg*h/L","n_observed_concentrations":6,"pk_source_sha256":SYNTH_SHA,
       "integration_method":"linear_trapezoid"},
    ]


class MultiOmicsPKReadinessTests(unittest.TestCase):
    def test_current_producer_without_collection_time_or_unit_cannot_predict_pk(self):
        with tempfile.TemporaryDirectory() as td:
            source=fixture(Path(td)/"exchange",time_columns=False)
            qc=inspect_pk_readiness(source)
            self.assertEqual(qc["n_independent_biological_subjects"],2)
            self.assertFalse(qc["temporal_numeric_metadata_ready"])
            self.assertEqual(qc["n_omic_assays_missing_collection_time"],3)
            self.assertFalse(qc["has_original_PK_dose_concentration_history"])
            self.assertIn("NUMERIC_COLLECTION_TIME_AND_EXPLICIT_UNIT_PER_ASSAY",qc["missing_for_E3_and_E4"])
            self.assertNotIn("PAT-001",repr(qc))
            with self.assertRaisesRegex(NCATargetError,"Unquantified"):
                qualify_baseline_omic_nca_join(source,targets(),cohort_time_axis_confirmed=True)

    def test_true_temporal_assays_do_not_count_technical_repeats_as_subjects(self):
        with tempfile.TemporaryDirectory() as td:
            p=fixture(Path(td)/"exchange",multiple_modalities=True)
            qc=inspect_pk_readiness(p)
            self.assertTrue(qc["temporal_numeric_metadata_ready"])
            self.assertEqual(qc["n_original_unique_specimens"],3)
            self.assertEqual(qc["n_independent_biological_subjects"],2)
            self.assertEqual(qc["n_subjects_with_two_or_more_distinct_numeric_collection_times"],1)
            self.assertEqual(qc["n_extra_technical_repeat_assays_not_independent_biology"],0)
            self.assertFalse(qc["clinical_or_PK_effect_validated"])

    def test_NCA_physical_units_window_and_source_provenance(self):
        qc=qualify_nca_targets(targets())
        self.assertEqual(qc["n_by_parameter"],{"Cmax":1,"Tmax":1,"AUC0last":1})
        self.assertEqual(qc["n_unique_original_subjects"],1)
        self.assertNotIn("PAT-001",repr(qc))
        bad=targets()
        bad[1]["unit"]="mg/L"
        with self.assertRaisesRegex(NCATargetError,"unit"):
            qualify_nca_targets(bad)
        bad=targets()
        bad[2]["window_end_h"]=20
        with self.assertRaisesRegex(NCATargetError,"incompatible PK time windows"):
            qualify_nca_targets(bad)
        bad=targets()
        bad[0]["pk_source_sha256"]="missing"
        with self.assertRaisesRegex(NCATargetError,"SHA256"):
            qualify_nca_targets(bad)

    def test_genuinely_pre_dose_omic_NCA_join_requires_same_study_origin(self):
        with tempfile.TemporaryDirectory() as td:
            p=fixture(Path(td)/"exchange",time_columns=True)
            with self.assertRaisesRegex(NCATargetError,"crosswalk"):
                qualify_baseline_omic_nca_join(p,targets(),cohort_time_axis_confirmed=False)
            qc=qualify_baseline_omic_nca_join(p,targets(),cohort_time_axis_confirmed=True)
            self.assertEqual(qc["n_nca_endpoints_linked_to_pre_dose_omic"],3)
            self.assertEqual(qc["n_original_independent_subjects_with_nca"],1)
            self.assertTrue(qc["group_effect_or_qvalue_is_not_an_individual_PK_covariate"])

    def test_no_using_omics_collected_after_dose_for_early_Cmax(self):
        with tempfile.TemporaryDirectory() as td:
            p=fixture(Path(td)/"exchange",time_columns=True)
            # The only available original state is time 0: change it to after dose.
            with (p/"samples.csv").open(newline="",encoding="utf8") as f:
                rows=list(csv.reader(f))
            rows[1][5]="20"
            rows[2][5]="21"
            write_csv(p/"samples.csv",rows[0],rows[1:])
            mf=json.loads((p/"manifest.json").read_text())
            mf["files"]["samples"]["sha256"]=hashlib.sha256((p/"samples.csv").read_bytes()).hexdigest()
            (p/"manifest.json").write_text(json.dumps(mf))
            qc=qualify_baseline_omic_nca_join(p,targets(),cohort_time_axis_confirmed=True)
            self.assertEqual(qc["n_nca_endpoints_linked_to_pre_dose_omic"],0)
            self.assertEqual(qc["n_subjects_without_pre_dose_omic"],1)

    def test_cross_subject_same_specimen_collision_refused(self):
        with tempfile.TemporaryDirectory() as td:
            p=fixture(Path(td)/"exchange")
            with (p/"samples.csv").open(newline="",encoding="utf8") as f:
                rows=list(csv.reader(f))
            rows[3][1]="S0"
            write_csv(p/"samples.csv",rows[0],rows[1:])
            mf=json.loads((p/"manifest.json").read_text())
            mf["files"]["samples"]["sha256"]=hashlib.sha256((p/"samples.csv").read_bytes()).hexdigest()
            (p/"manifest.json").write_text(json.dumps(mf))
            with self.assertRaisesRegex(ExchangeError,"multiple subjects"):
                inspect_pk_readiness(p)

if __name__=="__main__":
    unittest.main()
