"""Real-sample GEO GSE49012: miRNA versus measured shared CYP3A activity.

Vuppalanchi et al. 2013, DOI 10.1371/journal.pone.0074471.
Source: family SOFT with each original liver specimen's GEO activity metadata
and source-normalized miRNA feature values. 22 cirrhosis / 12 controls are
miRNA-array samples; the larger paper activity cohort is NOT interchangeable.

No CYP3A5-specific activity, CYP3A5 PGx, full RNA-seq, proteomics or clinical PK
is inferred from this dataset. Disease is a strong confounder, so every miRNA
increment is assessed above the explicit cirrhosis baseline. No raw GEO sample
rows, clinical metadata or individual predictions are ever returned.
"""
from __future__ import annotations

from collections import Counter
from hashlib import sha256
from math import exp, isfinite, log, log1p
from random import Random
from statistics import mean, median
import gzip
import re
from typing import Mapping, Sequence

from omics2pk.functional_activity_benchmark import _ridge_predict


class GEO49012Error(ValueError):
    pass


SERIES = "GSE49012"
DOI = "10.1371/journal.pone.0074471"
SOURCE_PUBMED = "24058572"
EXPECTED_PLATFORM = "GPL17470"
EXPECTED_GROUPS = {"DISEASE": 22, "NORMAL": 12}
FEATURES = {
    "mir27b_206_external_literature": ("hsa-miR-27b", "hsa-miR-206"),
    "mir103_107_external_literature_negative_control": ("hsa-miR-103", "hsa-miR-107"),
    "mir21_130a_142_external_literature": ("hsa-miR-21", "hsa-miR-130a", "hsa-miR-142"),
    "mir27b_only": ("hsa-miR-27b",),
    # GEO 2013 platform makes the correct strand identity critical:
    # these are NOT synonyms for Liu's miR-27b-3p / miR-21-5p.
    "mir142_5p_external_literature": ("hsa-miR-142-5p",),
    "mir142_3p_distinct_strand_control": ("hsa-miR-142-3p",),
    "mir27b_star_distinct_strand_control": ("hsa-miR-27b*",),
    "mir21_star_distinct_strand_control": ("hsa-miR-21*",),
    "mir130a_star_distinct_strand_control": ("hsa-miR-130a*",),
}
FULL_PANEL = (
    "hsa-miR-27b", "hsa-miR-206", "hsa-miR-21", "hsa-miR-130a", "hsa-miR-142",
)
MAX_SOFT_BYTES = 14_000_000


def _require(condition: bool, message: str):
    if not condition:
        raise GEO49012Error(message)


def _decode_source(source: bytes) -> str:
    _require(isinstance(source, bytes) and len(source) <= MAX_SOFT_BYTES,
             "Original GEO SOFT source missing or too large")
    try:
        source = gzip.decompress(source) if source[:2] == b"\x1f\x8b" else source
        return source.decode("utf-8-sig")
    except (OSError, UnicodeDecodeError) as exc:
        raise GEO49012Error("Invalid original GEO SOFT compressed/text bytes") from exc


def parse_original_geo_soft(source: bytes) -> list[dict]:
    """Source-aligned GEO samples; NEVER match donor metadata by row position."""
    text = _decode_source(source)
    _require(len(text) < MAX_SOFT_BYTES * 4, "GEO expansion exceeds size limit")
    _require("^SERIES = GSE49012" in text, "Wrong GEO accession")
    _require("24058572" in text, "Original paper PMID not found in source")
    samples: list[dict] = []
    current: dict | None = None
    table = False
    observed_accessions: set[str] = set()
    for raw in text.splitlines():
        line = raw.rstrip("\r")
        if line.startswith("^SAMPLE = "):
            accession = line.partition("=")[2].strip()
            _require(re.fullmatch(r"GSM\d{6,9}", accession) and
                     accession not in observed_accessions,
                     "Duplicate or invalid original GEO sample accession")
            observed_accessions.add(accession)
            current = {"accession": accession, "metadata": {}, "intensity": {},
                       "title": None, "platform": None, "tissue": None,
                       "organism": None}
            samples.append(current)
            table = False
            continue
        if current is None:
            continue
        if line.startswith("!sample_table_begin"):
            table = True
            continue
        if line.startswith("!sample_table_end"):
            table = False
            continue
        if table:
            if line.startswith("ID_REF\t") or line.startswith("#") or not line.strip():
                continue
            cells=line.split("\t")
            _require(len(cells)>=2, "Invalid original GEO sample miRNA table")
            fid=cells[0].strip().strip('"')
            value=cells[1].strip().strip('"')
            if fid in current["intensity"]:
                raise GEO49012Error("Duplicate original miRNA feature ID")
            if fid:
                current["intensity"][fid] = value
            continue
        for key, field in (
            ("!Sample_title = ", "title"),
            ("!Sample_platform_id = ", "platform"),
            ("!Sample_source_name_ch1 = ", "tissue"),
            ("!Sample_organism_ch1 = ", "organism"),
        ):
            if line.startswith(key):
                current[field] = line[len(key):].strip()
        prefix="!Sample_characteristics_ch1 = "
        if line.startswith(prefix):
            chunk=line[len(prefix):].strip()
            k,sep,v=chunk.partition(":")
            if sep:
                k=k.casefold().strip()
                _require(k not in current["metadata"],
                         "Duplicate original GEO characteristic " + k)
                current["metadata"][k]=v.strip()
    _require(len(samples)==34, "Need exactly 34 original GEO miRNA liver assays")
    titles=[r["title"] for r in samples]
    _require(all(isinstance(v,str) and bool(v.strip()) for v in titles)
             and len(set(titles))==34, "Missing/duplicate biological donor titles")
    _require(all(r["platform"] == EXPECTED_PLATFORM for r in samples),
             "Unexpected GEO miRNA assay platform")
    _require(all(r["organism"] == "Homo sapiens" for r in samples),
             "Non-human or unspecified source sample")
    _require(all("liver" in str(r["tissue"]).lower() for r in samples),
             "GEO tissue origin is not liver")
    status_counts=Counter()
    for r in samples:
        meta=r["metadata"]
        st=meta.get("status","").strip().upper()
        if st in {"CONTROL","HEALTHY"}: st="NORMAL"
        _require(st in EXPECTED_GROUPS,
                 "Source cirrhosis/control status missing or unrecognized")
        r["status"] = st
        status_counts[st]+=1
        outcome = meta.get("hepatic cyp3a activity")
        _require(isinstance(outcome,str) and outcome.strip(),
                 "Original per-donor measured CYP3A activity metadata absent")
        try:
            activity=float(outcome)
        except (ValueError,OverflowError):
            raise GEO49012Error("Invalid original measured CYP3A activity") from None
        _require(isfinite(activity) and activity > 0,
                 "Observed shared CYP3A activity must be strictly positive")
        r["activity"] = activity
        _require(len(r["intensity"]) >= 400, "GEO miRNA assay absent or truncated")
    _require(dict(status_counts) == EXPECTED_GROUPS,
             "Original miRNA cohort is 22 cirrhosis and 12 normal liver donors")
    return samples


def _folds(records: list[dict], n_splits: int = 5) -> list[list[int]]:
    """Predeclared within-disease stratification; no outcome labels used."""
    folds=[[] for _ in range(n_splits)]
    for group in ("DISEASE","NORMAL"):
        idx=sorted((i for i,r in enumerate(records) if r["status"]==group),
                   key=lambda i:records[i]["accession"])
        Random(20261104 if group=="DISEASE" else 20261105).shuffle(idx)
        for j,i in enumerate(idx):
            folds[j % n_splits].append(i)
    return folds


def compare_geo_mirna_vs_cyp3a_activity(
    records: Sequence[Mapping[str, object]],
) -> dict:
    """Outer held-donor MAE, same 34 people, prespecified external miRNA panel.

    The feature provenance is from Liu et al. 2016 (separate liver study).
    Public matrix is source-processed, not prospective freeze-from-raw.
    This is a secondary measured *shared CYP3A* E2 experiment.
    """
    rec=list(records)
    _require(len(rec)==34, "34 original GEO samples required")
    _require(len({r["accession"] for r in rec})==34, "Original GEO sample identities not independent")
    _require(dict(Counter(r["status"] for r in rec))==EXPECTED_GROUPS,
             "Original condition counts mismatch")
    rec.sort(key=lambda r:r["accession"])
    feature_arms={"disease_status_only":(), **FEATURES, "five_literature_miRNAs":FULL_PANEL}
    available={}
    omitted={}
    for arm,features in feature_arms.items():
        if all(f in r["intensity"] for r in rec for f in features):
            try:
                for r in rec:
                    for f in features:
                        x=float(r["intensity"][f])
                        if not isfinite(x) or x < 0:
                            raise ValueError()
            except (TypeError,ValueError,OverflowError):
                omitted[arm]="INVALID_OR_MISSING_ORIGINAL_FEATURE_VALUE"
                continue
            available[arm]=features
        else:
            omitted[arm]="PRESPECIFIED_MIRNA_FEATURE_ABSENT"
    _require("disease_status_only" in available, "Disease status baseline unavailable")
    # A whole biological-miRNA exploratory arm is predeclared independently
    # of activity, requiring exactly source-native mature/probe feature IDs.
    # Exclude S_Control* and any IDs without a human miRNA prefix.
    biological_ids=tuple(sorted(
        f for f in rec[0]["intensity"]
        if re.fullmatch(r"hsa-(?:miR|let)-[A-Za-z0-9*_.-]+", f)))
    X_mirna=None
    optional_pca=None
    pca_arms=("disease_plus_all_miRNA_PCA2","disease_plus_all_miRNA_PCA4")
    if len(biological_ids) < 50:
        for name in pca_arms:
            omitted[name]="ABSTAIN_TOO_FEW_SOURCE_BIOLOGICAL_MIRNA_PROBES"
    else:
        try:
            import numpy as np
            from sklearn.impute import SimpleImputer
            from sklearn.preprocessing import StandardScaler
            from sklearn.decomposition import PCA
        except ImportError:
            for name in pca_arms:
                omitted[name]="OPTIONAL_NUMPY_SKLEARN_DEPENDENCIES_UNAVAILABLE"
        else:
            # No center/outcome-driven feature filtering. Median imputation,
            # standardization and PCA are fit inside each outer TRAIN split.
            def source_value(raw):
                try:
                    val=float(raw)
                    return log1p(val) if isfinite(val) and val >= 0 else float("nan")
                except (ValueError,TypeError,OverflowError):
                    return float("nan")
            X_mirna=np.asarray([
                [source_value(r["intensity"].get(f)) for f in biological_ids]
                for r in rec
            ],dtype=float)
            if np.isfinite(X_mirna).sum() < .6 * X_mirna.size:
                for name in pca_arms:
                    omitted[name]="INSUFFICIENT_SOURCE_MIRNA_VALUE_COVERAGE"
            else:
                for name in pca_arms:
                    available[name]=()
                optional_pca=(np,SimpleImputer,StandardScaler,PCA)
    folds=_folds(rec)
    errors={arm:[] for arm in ("log_geometric_train_mean", *available)}
    by_group={arm:{"NORMAL":[],"DISEASE":[]} for arm in errors}
    n_pred={arm:0 for arm in errors}
    for held in folds:
        held=set(held)
        train=[i for i in range(34) if i not in held]
        test=[i for i in range(34) if i in held]
        ytrain=[log(rec[i]["activity"]) for i in train]
        for arm in errors:
            if arm=="log_geometric_train_mean":
                guess=[exp(mean(ytrain))]*len(test)
            elif arm in pca_arms:
                np,SimpleImputer,StandardScaler,PCA=optional_pca
                rank=2 if arm.endswith("PCA2") else 4
                imputer=SimpleImputer(strategy="median")
                t_train=imputer.fit_transform(X_mirna[train])
                t_test=imputer.transform(X_mirna[test])
                scaler=StandardScaler()
                z_train=scaler.fit_transform(t_train)
                z_test=scaler.transform(t_test)
                pca=PCA(n_components=rank,svd_solver="full")
                p_train=pca.fit_transform(z_train)
                p_test=pca.transform(z_test)
                x_train=[[float(rec[i]["status"]=="DISEASE"), *v]
                         for i,v in zip(train,p_train)]
                x_test=[[float(rec[i]["status"]=="DISEASE"), *v]
                        for i,v in zip(test,p_test)]
                guess=_ridge_predict(x_train,ytrain,x_test)
            else:
                features=available[arm]
                def vector(i):
                    r=rec[i]
                    return [
                        float(r["status"]=="DISEASE"),
                        *(log1p(float(r["intensity"][f])) for f in features),
                    ]
                guess=_ridge_predict([vector(i) for i in train],ytrain,
                                     [vector(i) for i in test])
            for i,pred in zip(test,guess):
                loss=abs(pred-rec[i]["activity"])
                errors[arm].append(loss)
                by_group[arm][rec[i]["status"]].append(loss)
                n_pred[arm]+=1
    _require(all(v==34 for v in n_pred.values()),
             "Fold split lost original biological samples")
    maes={arm:round(mean(loss),8) for arm,loss in errors.items()}
    medae={arm:round(median(loss),8) for arm,loss in errors.items()}
    base=maes["disease_status_only"]
    def quantile(vals, frac):
        vals=sorted(vals)
        pos=(len(vals)-1)*frac
        j=int(pos)
        return vals[j]+(vals[min(j+1,len(vals)-1)]-vals[j])*(pos-j)
    # The heldout loss per disease group is aligned across all arms because
    # the same sample partitions are used. This stratified paired bootstrap
    # resamples precomputed heldout losses only; it is NOT a refitted bootstrap
    # nor a calibrated inference interval for external population performance.
    intervals={}
    for idx,arm in enumerate(errors):
        if arm=="disease_status_only":
            continue
        group_deltas={
            group:[loss-ref for loss,ref in zip(
                by_group[arm][group],by_group["disease_status_only"][group])]
            for group in ("NORMAL","DISEASE")
        }
        rng=Random(20261104+idx)
        boot=[]
        for _ in range(1999):
            sampled=[
                group_deltas[g][rng.randrange(len(group_deltas[g]))]
                for g in ("NORMAL","DISEASE")
                for i in range(len(group_deltas[g]))
            ]
            boot.append(mean(sampled))
        intervals[arm]=[
            round(quantile(boot,.025),8),round(quantile(boot,.975),8)]
    source_distribution={
        group:{
            "n":len([r for r in rec if r["status"]==group]),
            "mean":round(mean(r["activity"] for r in rec
                              if r["status"]==group),6),
            "median":round(median(r["activity"] for r in rec
                                  if r["status"]==group),6),
            "min":round(min(r["activity"] for r in rec
                                if r["status"]==group),6),
            "max":round(max(r["activity"] for r in rec
                                if r["status"]==group),6),
        } for group in ("NORMAL","DISEASE")
    }
    # Mandatory numerical cross-check: the source paper says n=28 cirrhotic
    # livers had mean CYP3A CLint 0.28 (rounded to 2 decimals). GEO provides
    # n=22 cirrhotic miRNA profiles, reported to be a SUBSET of that cohort.
    # Since activity is nonnegative, 22 observed activities cannot sum to
    # more than the upper bound 28*(0.28+0.005)=7.98 if the measure is identical.
    # Failure is a genuine source-provenance discrepancy. NEVER "correct"
    # individual source values or convert units speculatively.
    source_disease_sum=sum(r["activity"] for r in rec if r["status"]=="DISEASE")
    maximum_published_disease_total=28 * (.28 + .005)
    disease_group_source_parity=(source_disease_sum <= maximum_published_disease_total)
    control_source_mean=mean(r["activity"] for r in rec if r["status"]=="NORMAL")
    control_group_source_parity=(abs(control_source_mean-.47) <= .005)
    source_parity_qualified=disease_group_source_parity and control_group_source_parity
    return {
        "status":("RESEARCH_E2_SUPPLIED_GEO_SOURCE_PARITY_POSSIBLE"
                  if source_parity_qualified else
                  "ABSTAIN_GEO_ACTIVITY_PUBLICATION_GROUP_PARITY_FAILURE"),
        "source_geo":SERIES,"source_doi":DOI,"source_pubmed":SOURCE_PUBMED,
        "n_original_unique_liver_biological_samples":34,
        "condition_counts":EXPECTED_GROUPS,
        "source_publication_parity_audit":{
            "paper_cirrhosis_n":28,
            "paper_cirrhosis_mean_rounded":0.28,
            "GEO_cirrhosis_n":22,
            "GEO_cirrhosis_activity_sum":round(source_disease_sum,6),
            "maximum_compatible_source_total_with_mean_rounding":round(maximum_published_disease_total,6),
            "disease_subset_of_paper_nonnegative_activity_possible":disease_group_source_parity,
            "paper_control_n":12,
            "paper_control_mean_rounded":0.47,
            "GEO_control_n":12,
            "GEO_control_mean":round(control_source_mean,6),
            "paper_control_mean_consistent":control_group_source_parity,
            "source_endpoint_paper_parity_qualified":source_parity_qualified,
            "same_endpoint_and_original_donor_subset_assumed_for_parity_check":True,
            "biological_superiority_inference_permitted":False,
            "unverified_outcome_units_or_original_sample_crosswalk_not_silently_repaired":True,
        },
        "n_outer_sample_heldout_folds":len(folds),
        "source_molecular_platform":EXPECTED_PLATFORM,
        "target":"original_sample_metadata_midazolam_shared_CYP3A_in_vitro_activity",
        "target_units":"source_PLoS_reported_mL_min_inverse_mg_protein_inverse",
        "not_CYP3A5_specific":True,"no_CYP3A5_genotype":True,
        "no_proteomics_per_sample":True,"no_human_drug_PK":True,
        "experiment_not_prospective_causal":True,
        "specific_miRNA_features_declared_without_source_outcome_selection":True,
        "mature_and_star_strands_not_interchangeable":True,
        "CYP3A5_isoform_specific_function_not_measured":True,
        "source_n_biological_miRNA_feature_IDs":len(biological_ids),
        "source_control_probes_excluded_from_unsupervised_PCA":True,
        "PCA_imputation_scaling_and_components_fit_within_outer_training_only":True,
        "miRNA_feature_sets_used":{k:list(v) for k,v in available.items()},
        "miRNA_feature_arms_abstained":omitted,
        "heldout_MAE_by_arm_source_activity_units":maes,
        "heldout_median_absolute_error_by_arm_source_activity_units":medae,
        "observed_source_activity_aggregate_by_disease":source_distribution,
        "exploratory_stratified_paired_delta_MAE_bootstrap_95pct_vs_disease":intervals,
        "paired_bootstrap_resamples_heldout_losses_not_refitted_models":1999,
        "external_population_uncertainty_calibrated":False,
        "heldout_MAE_by_arm_and_disease_status":{
            arm:{group:round(mean(loss),8) for group,loss in samples.items()}
            for arm,samples in by_group.items()},
        "delta_MAE_vs_disease_baseline":{
            arm:round(mae-base,8) for arm,mae in maes.items()
            if arm!="disease_status_only"},
        "source_upstream_normalization_not_train_frozen":True,
        "train_only_model_scaling_and_outcome_fits":True,
        "n_original_individual_rows_or_predictions_exported":0,
        "independent_center_or_external_cohort_validation":False,
    }
