"""Actual GNPS archive -> genuine midazolam PK / CYP3A4-5 research benchmark.

The source ZIP MUST be supplied locally by an authorized user. Never commit,
upload, log, or return original individual molecular intensities, names, sample
identifiers, clinical outcomes, genotypes or per-person model predictions.
The pre-dose MS table is a GNPS QIIME2-derived result, NOT verified identical
to the authors' original MZmine GNPS_quant.csv. Research only, not dosing.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
import csv
from hashlib import file_digest
from io import BytesIO, TextIOWrapper
import json
from pathlib import Path
from zipfile import ZipFile, BadZipFile

import numpy as np

GNPS_TASK = "b3966445abad4658a3cdc63c8198a666"
# Exact user-supplied original public output archive; no patient rows are stored.
ARCHIVE_SHA256 = "2047ea2032d73e3b48e2696fb1f652a7dedaa6f8cd5d40d0926d810bdd46c026"
N_FEATURES, N_SAMPLES = 6349, 1198
RANDOM_SEED = 20261019


class GNPSOriginalSourceError(ValueError):
    pass


def _observed(value):
    return str(value or "").strip() not in ("", "NA", "N/A", "nan", "NaN", "None", "not applicable")


def _number(value):
    try:
        result = float(value)
        return result if np.isfinite(result) else None
    except (TypeError, ValueError, OverflowError):
        return None


def load_original_source(archive):
    """Strict exact-source import with aggregate-only QC and private arrays."""
    import h5py
    from scipy.sparse import csr_matrix

    path = Path(archive).resolve()
    if not path.is_file() or not 50_000_000 < path.stat().st_size < 150_000_000:
        raise GNPSOriginalSourceError("Expected intact source archive, not raw clinical records")
    with path.open("rb") as f:
        digest = file_digest(f, "sha256").hexdigest()
    if digest != ARCHIVE_SHA256:
        raise GNPSOriginalSourceError("Unrecognized original GNPS archive SHA-256")
    try:
        with ZipFile(path) as z:
            if (len(z.namelist()) != 12 or any(
                    item.file_size > 100_000_000 or item.is_dir() or ".." in Path(item.filename).parts
                    for item in z.infolist())):
                raise GNPSOriginalSourceError("Original GNPS output member inventory changed")
            params = z.read("params.xml").decode("utf-8")
            if (GNPS_TASK not in params or "GNPS_quant.csv" not in params
                    or "QUANT_FILE_NORM" not in params):
                raise GNPSOriginalSourceError("Original GNPS task/quantification source not verified")
            metadata = list(csv.DictReader(
                TextIOWrapper(z.open("qiime2_output/qiime2_metadata.tsv"),
                              encoding="utf-8-sig"), delimiter="\t"))
            mandatory = (
                "sample_name", "ATTRIBUTE_Subject", "ATTRIBUTE_Study_Day",
                "ATTRIBUTE_Time_Point_Mins", "ATTRIBUTE_Sample_Type", "MS_SampleType",
                "PharmKinetics_CYP3A4", "PharmGeno_CYP3A4", "PharmGeno_CYP3A5",
                "Day0_Vitals_weight_kg", "ATTRIBUTE_sex",
            )
            if len(metadata) != 1440 or not metadata or not set(mandatory) <= set(metadata[0]):
                raise GNPSOriginalSourceError("GNPS original per-assay metadata schema changed")
            metadata_by_sample = defaultdict(list)
            for row in metadata:
                metadata_by_sample[row["sample_name"]].append(row)
            if len(metadata_by_sample) != 1202:
                raise GNPSOriginalSourceError("Original GNPS analytical sample set changed")
            # Metadata repeats several analytical sample names. Never treat
            # copies as biological persons, and require clinical/source agreement.
            for copies in metadata_by_sample.values():
                if any(any(other.get(key) != copies[0].get(key) for key in mandatory)
                       for other in copies[1:]):
                    raise GNPSOriginalSourceError("Conflicting duplicated original sample metadata")
            unique = {key: copies[0] for key, copies in metadata_by_sample.items()}
            original_predose = defaultdict(list)
            for assay, record in unique.items():
                if (record["MS_SampleType"] == "sample"
                        and record["ATTRIBUTE_Sample_Type"] == "blood_plasma"
                        and record["ATTRIBUTE_Study_Day"] == "1"
                        and record["ATTRIBUTE_Time_Point_Mins"] == "0"):
                    original_predose[record["ATTRIBUTE_Subject"]].append((assay, record))
            if len(original_predose) != 13 or any(len(x) != 1 for x in original_predose.values()):
                raise GNPSOriginalSourceError("T0 day1 plasma biological sample crosswalk changed")
            # Subject sorting only establishes stable array order AFTER explicit
            # original specimen-person join. It is never a row-order join.
            people = sorted(original_predose)
            assays = [original_predose[p][0][0] for p in people]
            source_rows = [original_predose[p][0][1] for p in people]
            if any(not all(_observed(row[col]) for col in mandatory[1:]) for row in source_rows):
                raise GNPSOriginalSourceError("Missing genuine baseline CYP/PK/PGx source")
            clinical = np.array([
                [float(r["Day0_Vitals_weight_kg"]),
                 float(r["ATTRIBUTE_sex"].lower() in ("male", "m"))]
                for r in source_rows], dtype=float)
            genotypes = np.array([
                [r["PharmGeno_CYP3A4"], r["PharmGeno_CYP3A5"]]
                for r in source_rows], dtype=object)
            if not np.isfinite(clinical).all() or np.any(clinical[:, 0] <= 0):
                raise GNPSOriginalSourceError("Unqualified actual individual clinical baseline")
            predose_pk = np.array(
                [_number(r["PharmKinetics_CYP3A4"]) for r in source_rows], dtype=float)
            # PK outcome values repeat across original analytical samples.
            # Collapse exactly once per biological individual and day.
            day9 = defaultdict(set)
            for r in metadata:
                if (r.get("MS_SampleType") == "sample"
                        and r.get("ATTRIBUTE_Study_Day") == "9"
                        and _observed(r.get("ATTRIBUTE_Subject"))):
                    y = _number(r.get("PharmKinetics_CYP3A4"))
                    if y is not None and y > 0:
                        day9[r["ATTRIBUTE_Subject"]].add(y)
            if any(len(y) > 1 for y in day9.values()):
                raise GNPSOriginalSourceError("Original same-person followup PK values disagree")
            followup_pk = np.array([
                next(iter(day9[p])) if p in day9 else np.nan for p in people], dtype=float)
            if (np.sum(np.isfinite(followup_pk)) != 12 or
                    not np.isfinite(predose_pk).all() or np.any(predose_pk <= 0)):
                raise GNPSOriginalSourceError("Actual day1/day9 PK endpoints missing/invalid")

            qza = z.read("qiime2_output/qiime2_table.qza")
            if not 50_000_000 < len(qza) < 75_000_000:
                raise GNPSOriginalSourceError("Original QIIME2 table is not intact")
            with ZipFile(BytesIO(qza)) as q:
                members = [name for name in q.namelist() if name.endswith("/data/feature-table.biom")]
                if len(members) != 1:
                    raise GNPSOriginalSourceError("Expected one genuine QIIME2 BIOM feature table")
                biom = q.read(members[0])
                if not 60_000_000 < len(biom) < 85_000_000:
                    raise GNPSOriginalSourceError("Invalid original BIOM uncompressed dimensions")
            with h5py.File(BytesIO(biom), "r") as h:
                original_sample_ids = [s.decode() for s in h["sample/ids"][:]]
                if (len(original_sample_ids) != N_SAMPLES
                        or len(set(original_sample_ids)) != N_SAMPLES
                        or len(h["observation/ids"]) != N_FEATURES):
                    raise GNPSOriginalSourceError("Original GNPS QIIME2 matrix dimensions mismatch")
                sample_to_col = {s: i for i, s in enumerate(original_sample_ids)}
                if not set(assays) <= set(sample_to_col):
                    raise GNPSOriginalSourceError("Original predose assay ID absent from QIIME2 BIOM")
                values = h["sample/matrix/data"][:]
                indices = h["sample/matrix/indices"][:]
                indptr = h["sample/matrix/indptr"][:]
                if (len(values) > 6_000_000 or not np.isfinite(values).all()
                        or np.any(values < 0)):
                    raise GNPSOriginalSourceError("Unqualified peak-area values in BIOM")
                matrix = csr_matrix(
                    (values, indices, indptr), shape=(N_SAMPLES, N_FEATURES))
                measured = matrix[[sample_to_col[s] for s in assays]].toarray()
            if measured.shape != (13, N_FEATURES):
                raise GNPSOriginalSourceError("Original measured molecular matrix mismatch")
            qc = {
                "original_GNPS_task": GNPS_TASK,
                "source_archive_sha256": digest,
                "original_biom_features": N_FEATURES,
                "original_biom_analytical_samples": N_SAMPLES,
                "baseline_plasma_biological_people": len(people),
                "n_available_day1_midazolam_CLF": int(np.isfinite(predose_pk).sum()),
                "n_available_day9_midazolam_CLF": int(np.isfinite(followup_pk).sum()),
                "n_features_detected_in_some_predose_plasma": int((measured > 0).any(axis=0).sum()),
                "original_t0_is_pre_admin_per_publication": True,
                "gnps_normalization": "QUANT_FILE_NORM=RowSum",
                "output_is_original_GNPS_derived_biom_not_authors_exact_GNPS_quant_csv": True,
                "original_participant_records_never_exported": True,
            }
            return measured, clinical, genotypes, {"day1": predose_pk, "day9": followup_pk}, qc
    except (OSError, KeyError, BadZipFile, UnicodeDecodeError) as exc:
        raise GNPSOriginalSourceError("Original GNPS source could not be decoded") from exc


def _one_model(train_matrix, test_matrix, log_target):
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import StandardScaler
    scaler = StandardScaler().fit(train_matrix)
    reg = Ridge(alpha=1.0).fit(scaler.transform(train_matrix), log_target)
    y = float(np.exp(reg.predict(scaler.transform(test_matrix))[0]))
    if not np.isfinite(y) or y <= 0:
        raise GNPSOriginalSourceError("Invalid research PK prediction")
    return y


def leave_one_person_out(original_matrix, clinical, genotypes, pk):
    """Fixed unsupervised PCA(2); outer fold-only fitting, n>=10 true persons.

    No global sample standardization, feature prevalence filter or PCA is
    fitted outside the outer TRAIN. This is a descriptive source benchmark,
    NOT a non-transductive independent clinical validation of MZmine/GNPS.
    """
    from sklearn.preprocessing import OneHotEncoder, StandardScaler
    from sklearn.decomposition import PCA
    x = np.asarray(original_matrix, dtype=float)
    c = np.asarray(clinical, dtype=float)
    g = np.asarray(genotypes, dtype=object)
    y = np.asarray(pk, dtype=float)
    n = len(y)
    if (not 10 <= n <= 30 or x.ndim != 2 or x.shape[0] != n
            or c.shape != (n, 2) or g.shape != (n, 2)
            or not np.isfinite(x).all() or np.any(x < 0)
            or not np.isfinite(c).all() or y.shape != (n,)
            or not np.isfinite(y).all() or np.any(y <= 0)):
        raise GNPSOriginalSourceError("Matched individual original inputs required")
    # Per-individual compositional total-peak normalization; no PK outcome
    # enters it. Ionization effects/batch remain major limitations.
    xx = np.log1p(x / np.maximum(x.sum(axis=1, keepdims=True), 1.0) * 1_000_000)
    result = {name: [] for name in ("clinical", "clinical_pgx", "clinical_pgx_metabolome_pca2")}
    for held in range(n):
        train = np.arange(n) != held
        test = ~train
        encoder = OneHotEncoder(handle_unknown="ignore")
        gt = encoder.fit_transform(g[train]).toarray()
        gh = encoder.transform(g[test]).toarray()
        btrain, btest = np.column_stack((c[train], gt)), np.column_stack((c[test], gh))
        target = np.log(y[train])
        result["clinical"].append(_one_model(c[train], c[test], target))
        result["clinical_pgx"].append(_one_model(btrain, btest, target))
        feature_mask = np.mean(x[train] > 0, axis=0) >= 0.5
        if feature_mask.sum() < 2:
            raise GNPSOriginalSourceError("Insufficient fold-independent predose features")
        norm = StandardScaler().fit(xx[train][:, feature_mask])
        z_train = norm.transform(xx[train][:, feature_mask])
        z_test = norm.transform(xx[test][:, feature_mask])
        pca = PCA(n_components=2, svd_solver="full").fit(z_train)
        result["clinical_pgx_metabolome_pca2"].append(_one_model(
            np.column_stack((btrain, pca.transform(z_train))),
            np.column_stack((btest, pca.transform(z_test))), target))
    return {k: np.asarray(v) for k, v in result.items()}


def summarize_models(measured, clinical, genotypes, target, *, shuffles=25):
    y = np.asarray(target, dtype=float)
    good = np.isfinite(y)
    x = np.asarray(measured)[good]
    c = np.asarray(clinical)[good]
    g = np.asarray(genotypes)[good]
    y = y[good]
    preds = leave_one_person_out(x, c, g, y)
    models = {}
    for model, prediction in preds.items():
        delta = prediction - y
        models[model] = {
            "mae": round(float(np.mean(abs(delta))), 6),
            "rmse": round(float(np.sqrt(np.mean(delta ** 2))), 6),
            "bias": round(float(np.mean(delta)), 6),
            "p20": round(float(np.mean(np.abs(delta/y) <= 0.2)), 6),
        }
    d = (abs(preds["clinical_pgx_metabolome_pca2"]-y)
         - abs(preds["clinical_pgx"]-y))
    rng = np.random.default_rng(RANDOM_SEED)
    bootstrap = np.mean(d[rng.integers(0, len(y), size=(1999, len(y)))], axis=1)
    negative = []
    for _ in range(shuffles):
        permuted = leave_one_person_out(x[rng.permutation(len(y))], c, g, y)
        negative.append(float(np.mean(abs(permuted["clinical_pgx_metabolome_pca2"] - y))))
    actual = models["clinical_pgx_metabolome_pca2"]["mae"]
    return {
        "n_unique_original_volunteers": len(y),
        "models": models,
        "incremental_omic_delta_mae": round(float(d.mean()), 6),
        "incremental_omic_delta_mae_pct": round(
            float(100*d.mean()/np.mean(abs(preds["clinical_pgx"]-y))), 4),
        "paired_bootstrap_95pct_descriptive_no_refit": [
            round(float(z), 6) for z in np.percentile(bootstrap, [2.5, 97.5])],
        "n_original_subjects_with_lower_omic_absolute_error": int((d < 0).sum()),
        "negative_control_n_shuffled_subject_mappings": len(negative),
        "negative_control_median_shuffled_omic_mae": round(float(np.median(negative)), 6),
        "negative_control_fraction_shuffled_better_or_equal": round(
            float(np.mean(np.asarray(negative) <= actual)), 4),
        "source_scale_units_physically_unverified": True,
        "independent_external_validation": False,
        "small_sample_exploratory_only": True,
    }


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", required=True, help="Authorized private original GNPS ZIP")
    args = parser.parse_args(argv)
    matrix, clinical, genotype, pk, qc = load_original_source(args.archive)
    # Only aggregates printed; original records, subject IDs and predictions
    # remain within this local process and are never saved/uploaded.
    out = {
        "status": "ORIGINAL_GNPS_DERIVED_QIIME_BIOM_TO_TRUE_HUMAN_PK_INTERNAL_ONLY",
        "original_source_aggregate_qc": qc,
        "day1_midazolam_oral_CLF": summarize_models(matrix, clinical, genotype, pk["day1"]),
        "day9_midazolam_oral_CLF": summarize_models(matrix, clinical, genotype, pk["day9"]),
        "does_not_establish_improved_PK_from_metabolomics": True,
        "privacy": "aggregate-only, no original subject rows in output",
    }
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
