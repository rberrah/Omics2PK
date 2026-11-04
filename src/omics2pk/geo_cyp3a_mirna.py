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
from statistics import mean
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
    base=maes["disease_status_only"]
    return {
        "status":"PUBLIC_ORIGINAL_GEO_LIVER_MIRNA_SHARED_CYP3A_ACTIVITY_E2",
        "source_geo":SERIES,"source_doi":DOI,"source_pubmed":SOURCE_PUBMED,
        "n_original_unique_liver_biological_samples":34,
        "condition_counts":EXPECTED_GROUPS,
        "n_outer_sample_heldout_folds":len(folds),
        "source_molecular_platform":EXPECTED_PLATFORM,
        "target":"original_sample_metadata_midazolam_shared_CYP3A_in_vitro_activity",
        "target_units":"source_PLoS_reported_mL_min_inverse_mg_protein_inverse",
        "not_CYP3A5_specific":True,"no_CYP3A5_genotype":True,
        "no_proteomics_per_sample":True,"no_human_drug_PK":True,
        "experiment_not_prospective_causal":True,
        "omic_features_selected_a_priori_from_independent_Liu2016_study":True,
        "miRNA_feature_sets_used":{k:list(v) for k,v in available.items()},
        "miRNA_feature_arms_abstained":omitted,
        "heldout_MAE_by_arm_source_activity_units":maes,
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
