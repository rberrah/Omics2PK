"""Fail-closed multi-CYP same-donor RNA/protein/activity benchmark.

Designed for legitimately obtained ORIGINAL individual liver tables such as
Zhang 2016 (DOI 10.1124/jpet.116.233635; published n=100) or Ohtsuki 2012
(DOI 10.1124/dmd.111.042259; published n=17). No source rows are bundled.

Input is long-form donor × CYP. Every included CYP must contain the complete
published donor count for that source and one authenticated row per donor with
mRNA, protein and measured probe-substrate activity on the same liver.

Research only. Outcome is in-vitro drug metabolism, not clinical plasma PK.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from hashlib import sha256
from math import expm1, isfinite, log1p
from pathlib import Path
from random import Random
from statistics import mean

from omics2pk.epa_liver_cyp_validation import _ridge


SOURCE_PROFILES = {
    "10.1124/jpet.116.233635": {
        "label": "Zhang_2016_10CYP",
        "expected_donors_per_included_cyp": 100,
        "published_scope": "10 CYP mRNA + protein + probe intrinsic activities",
    },
    "10.1124/dmd.111.042259": {
        "label": "Ohtsuki_2012_CYP",
        "expected_donors_per_included_cyp": 17,
        "published_scope": "CYP mRNA + absolute protein + probe activities",
    },
}

LABELS = (
    "donor_id", "source_doi", "row_kind", "cyp_isoform", "probe_substrate",
    "rna_value", "rna_unit", "protein_value", "protein_unit",
    "activity_value", "activity_unit",
)
SEED = 20261022


class MultiCYPSourceNotQualified(ValueError):
    pass


def _positive(row, key):
    try:
        value=float(row[key])
    except (KeyError,TypeError,ValueError):
        raise MultiCYPSourceNotQualified("Missing numeric source field: "+key) from None
    if not isfinite(value) or value<=0:
        raise MultiCYPSourceNotQualified("Source measurement must be strictly positive: "+key)
    return value


def _text(row,key):
    value=str(row.get(key,"")).strip()
    if not value:
        raise MultiCYPSourceNotQualified("Missing source field: "+key)
    return value


def source_contract(rows):
    """Validate an authenticated long-form same-donor source extract.

    A source may contain one or several CYPs. It may be a partial CYP extract
    from the paper, but every included CYP must include the exact published
    donor count so an arbitrary favorable subset cannot be modeled.
    """
    if not rows:
        raise MultiCYPSourceNotQualified("No original source rows")
    if any(set(row)!=set(LABELS) for row in rows):
        raise MultiCYPSourceNotQualified("Require exact donor/CYP assay source schema")
    dois={_text(row,"source_doi") for row in rows}
    if len(dois)!=1:
        raise MultiCYPSourceNotQualified("One original DOI per source file required")
    doi=next(iter(dois))
    if doi not in SOURCE_PROFILES:
        raise MultiCYPSourceNotQualified("Source DOI is not a preregistered same-donor study")
    profile=SOURCE_PROFILES[doi]
    by_cyp=defaultdict(list)
    seen=set()
    units={}
    for row in rows:
        if _text(row,"row_kind").lower()!="individual":
            raise MultiCYPSourceNotQualified("Pooled, summary, simulated or PBPK rows forbidden")
        donor=_text(row,"donor_id")
        cyp=_text(row,"cyp_isoform").upper()
        if not cyp.startswith("CYP") or len(cyp)>12:
            raise MultiCYPSourceNotQualified("Invalid CYP isoform label")
        key=(donor,cyp)
        if key in seen:
            raise MultiCYPSourceNotQualified("Duplicate donor × CYP source row")
        seen.add(key)
        probe=_text(row,"probe_substrate")
        rna=_positive(row,"rna_value")
        protein=_positive(row,"protein_value")
        activity=_positive(row,"activity_value")
        row_units=(
            _text(row,"rna_unit"), _text(row,"protein_unit"), _text(row,"activity_unit")
        )
        if cyp in units and units[cyp]!=row_units:
            raise MultiCYPSourceNotQualified("Units vary within one CYP source assay")
        units[cyp]=row_units
        by_cyp[cyp].append({
            "donor_id":donor, "probe_substrate":probe,
            "rna":rna, "protein":protein, "activity":activity,
        })
    expected=profile["expected_donors_per_included_cyp"]
    if not by_cyp:
        raise MultiCYPSourceNotQualified("No CYP assay rows")
    for cyp,records in by_cyp.items():
        donors={x["donor_id"] for x in records}
        probes={x["probe_substrate"] for x in records}
        if len(records)!=expected or len(donors)!=expected:
            raise MultiCYPSourceNotQualified(
                f"{cyp}: require exactly {expected} original donors; no post-hoc subset")
        if len(probes)!=1:
            raise MultiCYPSourceNotQualified(f"{cyp}: one prespecified probe assay per benchmark")
    return doi,profile,dict(by_cyp),units


def _predict(train, held, arm):
    if arm=="training_mean":
        return mean(x["activity"] for x in train)
    field={"RNA_only":"rna","protein_only":"protein"}.get(arm)
    if field:
        xtrain=[[log1p(x[field])] for x in train]
        xtest=[log1p(held[field])]
    elif arm=="RNA_plus_protein":
        xtrain=[[log1p(x["rna"]),log1p(x["protein"])] for x in train]
        xtest=[log1p(held["rna"]),log1p(held["protein"])]
    else:
        raise MultiCYPSourceNotQualified("Unknown benchmark arm")
    y=[log1p(x["activity"]) for x in train]
    pred=max(0.0,expm1(min(30.0,_ridge(xtrain,y,xtest))))
    if not isfinite(pred):
        raise MultiCYPSourceNotQualified("Nonfinite held-donor prediction")
    return pred


def _paired_ci(values, seed):
    rng=Random(seed)
    n=len(values)
    draws=sorted(mean(values[rng.randrange(n)] for _ in range(n)) for _ in range(1999))
    return [round(draws[49],6),round(draws[1949],6)]


def _evaluate_cyp(records, seed):
    records=sorted(records,key=lambda x:x["donor_id"])
    arms=("training_mean","RNA_only","protein_only","RNA_plus_protein")
    errors={arm:[] for arm in arms}
    for i,held in enumerate(records):
        train=[x for j,x in enumerate(records) if j!=i]
        for arm in arms:
            pred=_predict(train,held,arm)
            errors[arm].append(abs(pred-held["activity"]))
    delta_combo=[a-b for a,b in zip(errors["RNA_plus_protein"],errors["RNA_only"])]
    delta_protein=[a-b for a,b in zip(errors["protein_only"],errors["RNA_only"])]
    return {
        "n_held_donor_folds":len(records),
        "probe_substrate":records[0]["probe_substrate"],
        "mae_by_model":{arm:round(mean(v),6) for arm,v in errors.items()},
        "delta_MAE_RNA_plus_protein_minus_RNA":round(mean(delta_combo),6),
        "delta_MAE_protein_minus_RNA":round(mean(delta_protein),6),
        "paired_loss_bootstrap_CI95_RNA_plus_protein_minus_RNA":_paired_ci(delta_combo,seed),
        "n_donors_combo_lower_absolute_error_than_RNA":sum(x<0 for x in delta_combo),
        "n_donors_combo_higher_absolute_error_than_RNA":sum(x>0 for x in delta_combo),
    }


def evaluate_original_source(rows):
    doi,profile,by_cyp,units=source_contract(rows)
    results={}
    for offset,cyp in enumerate(sorted(by_cyp)):
        result=_evaluate_cyp(by_cyp[cyp],SEED+offset)
        result["rna_unit"]=units[cyp][0]
        result["protein_unit"]=units[cyp][1]
        result["activity_unit"]=units[cyp][2]
        results[cyp]=result
    return {
        "status":"LOCAL_ORIGINAL_MULTI_CYP_SOURCE_GATED_NOT_EXTERNAL_VALIDATION",
        "source_doi":doi,
        "source_profile":profile["label"],
        "published_scope":profile["published_scope"],
        "n_included_cyp_endpoints":len(results),
        "included_cyp_isoforms":sorted(results),
        "expected_complete_donors_per_included_cyp":profile["expected_donors_per_included_cyp"],
        "results_by_cyp":results,
        "primary_question":"does RNA+protein improve held-donor activity prediction over RNA alone",
        "model":"donor-LODO fixed ridge alpha1; log1p measurements; train-only scaling",
        "not_clinical_plasma_PK":True,
        "not_whole_transcriptome_or_whole_proteome":True,
        "not_external_validation":True,
        "no_cross_study_fusion":True,
        "original_individual_rows_or_predictions_exported":False,
    }


def analyze_authorized_source_csv(path):
    blob=Path(path).read_bytes()
    if not 512<=len(blob)<=2_000_000:
        raise MultiCYPSourceNotQualified("Unexpected authorized source file size")
    try:
        text=blob.decode("utf-8-sig")
    except UnicodeDecodeError:
        raise MultiCYPSourceNotQualified("Authorized source extract must be UTF-8 CSV") from None
    rows=list(csv.DictReader(text.splitlines()))
    result=evaluate_original_source(rows)
    result["authorized_local_source_sha256"]=sha256(blob).hexdigest()
    result["source_file_not_redistributed"]=True
    return result
