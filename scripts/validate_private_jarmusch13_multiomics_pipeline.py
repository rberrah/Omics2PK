"""Private Jarmusch MS1 → MultiOmics V1 → Omics2PK numerical parity gate.

The source measurements remain in a local, researcher-authorized ZIP. No
participant IDs, patient measurements, molecular values, genotype calls or
individual predictions are committed, logged, returned or uploaded to CI.

The prepared ZIP is a MultiOmics *input*, NOT by itself a real browser-produced
DOWNSTREAM_EXCHANGE_V1 export. Only --native-exchange checks such an export.
Research-only internal LOPO, not independent PK validation.
"""
from __future__ import annotations

import argparse
from hashlib import file_digest
from pathlib import Path
import csv
import json
from io import TextIOWrapper
from zipfile import ZipFile, BadZipFile

import numpy as np

from omics2pk.exchange import _csv, _file, inspect_exchange

ZIP_SHA256 = "95a556d7ca288da9cf41f45d58212e05bd8e4e531b05e830505ede37cef8f294"
ORIGINAL_GNPS_SHA256 = "2047ea2032d73e3b48e2696fb1f652a7dedaa6f8cd5d40d0926d810bdd46c026"
GNPS_TASK = "b3966445abad4658a3cdc63c8198a666"
FILES = frozenset((
    "README.txt", "manifest.json", "multiomics_ms1_wide.csv",
    "multiomics_ms1_peak_area_long.csv", "multiomics_metadata.csv",
    "ms1_features_annotations.csv", "clinical_midazolam_PRIVATE.csv",
))
N_SUBJECTS, N_FEATURES = 13, 6349


class JarmuschMS1PipelineError(ValueError):
    """Abstain on any failed identity, unit, timing or source-join gate."""


def _require(condition: bool, reason: str):
    if not condition:
        raise JarmuschMS1PipelineError(reason)


def _rows(z, name):
    with z.open(name) as fd:
        reader = csv.DictReader(TextIOWrapper(fd, encoding="utf-8-sig"))
        _require(reader.fieldnames is not None and
                 len(reader.fieldnames) == len(set(reader.fieldnames)),
                 "Invalid or duplicate table column names")
        rows = list(reader)
    _require(bool(rows) and all(None not in row for row in rows),
             "Missing/ragged original input table")
    return rows, reader.fieldnames


def read_prepared_private_zip(filename):
    """Return private arrays in memory and a safe aggregate-only QC record."""
    p = Path(filename)
    _require(p.is_file() and 100_000 < p.stat().st_size < 4_000_000,
             "Prepared source is absent or implausible")
    with p.open("rb") as fd:
        _require(file_digest(fd, "sha256").hexdigest() == ZIP_SHA256,
                 "Unrecognized researcher-supplied original-derived input SHA256")
    try:
        with ZipFile(p) as z:
            _require(set(z.namelist()) == FILES and
                     all(0 < item.file_size < 10_000_000 and
                         ".." not in Path(item.filename).parts
                         for item in z.infolist()),
                     "Original input member inventory changed")
            manifest = json.loads(z.read("manifest.json"))
            _require(manifest.get("source_uploaded_archive_sha256") == ORIGINAL_GNPS_SHA256 and
                     manifest.get("original_GNPS_task") == GNPS_TASK and
                     manifest.get("source_biological_people") == N_SUBJECTS and
                     manifest.get("MS1_features") == N_FEATURES,
                     "Original study provenance cannot be confirmed")
            samples, _ = _rows(z, "multiomics_metadata.csv")
            clinical, _ = _rows(z, "clinical_midazolam_PRIVATE.csv")
            wide, columns = _rows(z, "multiomics_ms1_wide.csv")
            long, _ = _rows(z, "multiomics_ms1_peak_area_long.csv")
            annotations, _ = _rows(z, "ms1_features_annotations.csv")
    except (OSError, UnicodeError, BadZipFile, KeyError, json.JSONDecodeError) as err:
        raise JarmuschMS1PipelineError("Unable to parse exact protected original source") from None

    _require(len(samples) == N_SUBJECTS and len(clinical) == N_SUBJECTS and
             len(wide) == N_FEATURES and len(annotations) == N_FEATURES and
             len(long) == N_FEATURES * N_SUBJECTS and
             columns[0] == "feature_id" and len(columns) == N_SUBJECTS+1,
             "Private source dimensions changed")
    subjects = {r["subject_id"]: r for r in samples}
    _require(len(subjects) == N_SUBJECTS and
             len({r["assay_id"] for r in samples}) == N_SUBJECTS and
             len({r["sample_id"] for r in samples}) == N_SUBJECTS and
             set(columns[1:]) == {r["assay_id"] for r in samples},
             "Source assay/biological identity ambiguous")
    _require(all(r["omic"] == "metabolomics" and
                 r["condition"] == "baseline_day1" and
                 r["timepoint"] == "T0" and
                 r["tissue"] == "blood_plasma" and
                 r["sample_type"] == "biological"
                 for r in samples),
             "Predose plasma study metadata changed")
    source_ids = sorted(subjects)  # deterministic only AFTER explicit key join
    assays = [subjects[s]["assay_id"] for s in source_ids]
    features = [r["feature_id"] for r in wide]
    _require(len(set(features)) == N_FEATURES and
             all(f.startswith("GNPS_MS1:") for f in features) and
             {r["feature_id"] for r in annotations} == set(features) and
             all(r["putative_annotation_warning"].startswith("GNPS MS2 tentative")
                 for r in annotations),
             "Unidentified MS1 ions cannot be silently mapped to true metabolites")

    try:
        matrix = np.asarray([[float(r[a]) for a in assays] for r in wide], dtype=float).T
        _require(np.isfinite(matrix).all() and np.all(matrix >= 0),
                 "Nonfinite or negative MS1 areas")
        feature_idx = {f:i for i,f in enumerate(features)}
        assay_idx = {a:i for i,a in enumerate(assays)}
        seen = set()
        for r in long:
            feature, assay = r["feature_id"], r["assay_id"]
            key = (feature, assay)
            _require(key not in seen and feature in feature_idx and assay in assay_idx,
                     "Duplicate or unknown MS1 long observation")
            seen.add(key)
            specimen = subjects[r["subject_id"]]
            _require(specimen["assay_id"] == assay and
                     specimen["sample_id"] == r["sample_id"] and
                     r["timepoint"] == "T0" and r["condition"] == "baseline_day1",
                     "MS1 long-format patient/sample crosswalk failed")
            value = float(r["peak_area"])
            _require(np.isfinite(value) and value >= 0 and
                     np.isclose(value, matrix[assay_idx[assay], feature_idx[feature]],
                                rtol=1e-9, atol=1e-9),
                     "Wide versus long MultiOmics ion-area numeric mismatch")
        _require(len(seen) == N_SUBJECTS*N_FEATURES,
                 "Not all declared original feature-assay pairs are present")

        people = {r["subject_id"]:r for r in clinical}
        _require(len(people) == N_SUBJECTS and set(people) == set(subjects),
                 "Clinical/Pgx cohort differs from molecular cohort")
        covariates, pgx, pk = [], [], []
        for subject in source_ids:
            record = people[subject]
            _require(record["assay_id"] == subjects[subject]["assay_id"] and
                     record["condition"] == "baseline_day1" and
                     record["timepoint"] == "T0",
                     "Original clinical and assay identity disagree")
            sex = record["sex"].strip().lower()
            _require(sex in ("male","female","m","f"),
                     "Original sex metadata unqualified")
            covariates.append([float(record["body_weight_kg"]),
                               float(sex in ("male","m"))])
            pgx.append([record["original_CYP3A4_call"],
                        record["original_CYP3A5_call"]])
            followup = record["day9_measured_midazolam_CLF"].strip()
            pk.append([float(record["day1_measured_midazolam_CLF"]),
                       float(followup) if followup not in ("", "NA", "N/A") else np.nan])
        covariates = np.asarray(covariates, dtype=float)
        pk = np.asarray(pk, dtype=float)
        _require(np.isfinite(covariates).all() and
                 np.all(covariates[:, 0] > 0) and
                 np.isfinite(pk[:, 0]).all() and np.all(pk[:, 0] > 0) and
                 np.isfinite(pk[:, 1]).sum() == 12 and
                 np.all(pk[np.isfinite(pk[:, 1]), 1] > 0),
                 "Original measured midazolam CL/F incomplete or nonpositive")
    except (KeyError, TypeError, ValueError, OverflowError) as err:
        if isinstance(err, JarmuschMS1PipelineError):
            raise
        raise JarmuschMS1PipelineError("Invalid original molecular/PK data") from None

    return {
        "matrix": matrix, "clinical": covariates,
        "genotypes": np.asarray(pgx, dtype=object), "pk": pk,
        "subject_ids": tuple(source_ids), "assay_ids": tuple(assays),
        "feature_ids": tuple(features),
    }, {
        "status": "ORIGINAL_DERIVED_PRIVATE_MULTIOMICS_INPUT_PARITY_CONFIRMED",
        "n_original_independent_volunteers": N_SUBJECTS,
        "n_ms1_features": N_FEATURES,
        "n_total_feature_assay_pairs_numeric_verified": len(seen),
        "n_nonzero_in_some_predose_sample": int((matrix > 0).any(axis=0).sum()),
        "n_day1_midazolam_CLF": int(np.isfinite(pk[:,0]).sum()),
        "n_day9_midazolam_CLF": int(np.isfinite(pk[:,1]).sum()),
        "actual_native_MultiOmics_browser_export_verified": False,
        "relative_day1_T0_not_absolute_date": True,
        "verified_metabolite_compound_identity": False,
        "original_participant_values_exported": False,
    }


def check_real_native_multiomics_export(native_dir, private_source):
    """Require numeric parity with a REAL independently exported V1 directory.

    This function does not fabricate a native exporter result. The caller must
    supply the actual downloaded `exchange/` directory from the MultiOmics UI.
    A valid schema alone does not establish export provenance or permit raw
    vs normalized quantities to be silently interchanged.
    """
    report = inspect_exchange(native_dir)
    _require(report["data_access"] == "private" and
             report["subjects"] == N_SUBJECTS and
             "metabolomics" in report["matrices"] and
             report["matrices"]["metabolomics"]["n_features"] == N_FEATURES,
             "Actual native MultiOmics cohort or privacy declaration inconsistent")
    root = Path(native_dir).resolve()
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    sample_spec = manifest["files"]["samples"]
    _, rows = _csv(_file(root, sample_spec),
                   ["subject_id","sample_id","assay_id","omic"])
    samples = {r["subject_id"]:r for r in rows if r["omic"]=="metabolomics"}
    _require(set(samples)==set(private_source["subject_ids"]) and
             len(samples)==N_SUBJECTS and
             all(samples[s]["assay_id"] == a
                 for s,a in zip(private_source["subject_ids"],
                                private_source["assay_ids"])),
             "Actual native exporter changed the source biological/assay join")
    matrix_records=[r for r in manifest["matrices"] if r["omic"]=="metabolomics"]
    _require(len(matrix_records)==1, "Actual native export has ambiguous metabolomic matrix")
    head, matrixrows = _csv(_file(root, matrix_records[0]["file"]), ["feature_id"])
    _require(set(head[1:])==set(private_source["assay_ids"]) and
             len(matrixrows)==N_FEATURES and
             {r["feature_id"] for r in matrixrows}==set(private_source["feature_ids"]),
             "Actual native metabolite feature/sample universe changed")
    original = {f:i for i,f in enumerate(private_source["feature_ids"])}
    assay = {a:i for i,a in enumerate(private_source["assay_ids"])}
    try:
        for r in matrixrows:
            for a in private_source["assay_ids"]:
                actual=float(r[a]); expect=private_source["matrix"][assay[a],original[r["feature_id"]]]
                _require(np.isfinite(actual) and
                         np.isclose(actual,expect,rtol=1e-9,atol=1e-9),
                         "Actual MultiOmics output changed MS1 values; normalization/source stage needs audit")
    except (ValueError,TypeError,KeyError):
        raise JarmuschMS1PipelineError("Non-numeric or mismatched actual MultiOmics ion areas") from None
    return {"status":"ACTUAL_NATIVE_MULTIOMICS_EXPORT_SOURCE_ID_AND_NUMERIC_PARITY_PASS",
            "n_source_assays_verified":N_SUBJECTS,"n_features_verified":N_FEATURES,
            "n_feature_assay_numeric_matches":N_SUBJECTS*N_FEATURES,
            "original_individual_measurements_exported":False,
            "producer_commit_declared":bool(manifest.get("software",{}).get("commit")),
            "upstream_research_validation_scope":report["validation"]}


def benchmark_real_pk(source):
    """Run the previously frozen genotype vs unsupervised metabolome benchmark."""
    from scripts.validate_local_jarmusch14_gnps_biom_pk import summarize_models
    x,c,g,pk = (source[k] for k in ("matrix","clinical","genotypes","pk"))
    # Real PK outcomes are held out by original PERSON, not ion/sample row.
    day1 = summarize_models(x,c,g,pk[:,0],shuffles=0)
    day9 = summarize_models(x,c,g,pk[:,1],shuffles=0)
    return {"day1":day1,"day9":day9,
            "status":"EXPLORATORY_ORIGINAL_CLINICAL_PK_NO_INCREMENTAL_OMICS_GAIN",
            "independent_external_validation":False}


def main(argv=None):
    cli=argparse.ArgumentParser(description="Private Jarmusch MultiOmics/MS1→PK QC only")
    cli.add_argument("--private-input-zip",required=True)
    cli.add_argument("--native-exchange",default=None,
                     help="Genuine locally exported MultiOmics exchange directory")
    a=cli.parse_args(argv)
    source, qc=read_prepared_private_zip(a.private_input_zip)
    out={"source_input_qc":qc,"pk":benchmark_real_pk(source)}
    if a.native_exchange is not None:
        out["native_export_parity"]=check_real_native_multiomics_export(a.native_exchange,source)
    print(json.dumps(out,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
