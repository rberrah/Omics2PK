"""Source-gated same-human-liver CYP3A4 transcriptome vs proteome benchmark.

Use only a legitimately obtained individual donor table from Jamwal 2018,
DOI 10.1021/acs.molpharmaceut.8b00159, n74 human donor livers.

Outcome is ORIGINAL measured MIDAZOLAM 1'-hydroxylation enzyme activity
in HLM (not clinical plasma AUC; not paper's PBPK simulations).
CYP3A4 qPCR normalized to 18S is one targeted transcript measurement,
NOT whole-transcriptome; CYP3A4 protein is one targeted proteomic assay,
NOT whole-proteome. Comparison tests whether the extra *single* protein
improves held-donor-out CYP3A4 activity relative to one RNA abundance.

No original person-level values bundled, predicted, downloaded or printed.
No fabricated "n74 original results" when supplementary does not contain
verified same-donor numerical assay measurements.
"""
from __future__ import annotations
from hashlib import sha256
from math import isfinite,log1p,expm1,sqrt
from pathlib import Path
from random import Random
from statistics import mean
import csv

from omics2pk.epa_liver_cyp_validation import _ridge

SOURCE_DOI="10.1021/acs.molpharmaceut.8b00159"
LABELS=(
    "donor_id","source_doi","row_kind",
    "age_years","sex","histology",
    "CYP3A4_mRNA_relative_18S",
    "CYP3A4_protein_pmol_mg",
    "midazolam_1hydroxylation_pmol_min_mg",
)
CLASS={"normal","nafl","nash"}
ARMS={
    "demographics_age_sex_histology":(0,1,2,3),
    "RNA_only":(4,),
    "protein_only":(5,),
    "RNA_plus_protein":(4,5),
    "clinical_plus_RNA":(0,1,2,3,4),
    "clinical_plus_protein":(0,1,2,3,5),
    "clinical_plus_RNA_and_protein":(0,1,2,3,4,5),
}
EXPECTED_N=74
SEED=20261018

class SourceNotQualified(ValueError):
    pass


def _positive(row,key):
    try: v=float(row[key])
    except (TypeError,ValueError,KeyError):
        raise SourceNotQualified("Missing original numeric "+key) from None
    if not isfinite(v) or v<=0:
        raise SourceNotQualified("Original donor assay not strictly positive: "+key)
    return v


def source_contract(rows):
    if len(rows)!=EXPECTED_N:
        raise SourceNotQualified("Require exactly original 74 independent liver donors")
    keys=set()
    data=[]
    counts={v:0 for v in CLASS}
    for r in rows:
        if set(r)!=set(LABELS):
            raise SourceNotQualified("Require precise original donor/assay/DOI schema")
        if r["source_doi"].strip()!=SOURCE_DOI or r["row_kind"].strip().lower()!="individual":
            raise SourceNotQualified("Source DOI mismatch or pooled/simulated PBPK non-donor row")
        donor=r["donor_id"].strip()
        if not donor or donor in keys:
            raise SourceNotQualified("Missing or duplicated original biological donor")
        keys.add(donor)
        age=_positive(r,"age_years")
        if age<18 or age>110:
            raise SourceNotQualified("Original human donor age outside allowed range")
        sex=r["sex"].strip().lower()
        if sex not in ("male","female","m","f"):
            raise SourceNotQualified("Original donor sex not recognized")
        group=r["histology"].strip().lower()
        if group not in CLASS:
            raise SourceNotQualified("Original NAFLD disease histology not recognized")
        counts[group]+=1
        mrna=_positive(r,"CYP3A4_mRNA_relative_18S")
        protein=_positive(r,"CYP3A4_protein_pmol_mg")
        activity=_positive(r,"midazolam_1hydroxylation_pmol_min_mg")
        # The three histological classes are encoded as 2 dummy variables;
        # never use ordinal class levels as numeric risk score.
        x=(age,float(sex in ("male","m")),float(group=="nafl"),float(group=="nash"),
           log1p(mrna),log1p(protein))
        data.append((donor,x,activity))
    if counts!={"normal":24,"nafl":26,"nash":24}:
        raise SourceNotQualified("Original study groups not 24 normal / 26 NAFL / 24 NASH")
    return sorted(data),counts


def _predict(train,held,name):
    if name=="training_mean":
        return mean(r[2] for r in train)
    fields=ARMS[name]
    xx=[[r[1][k] for k in fields] for r in train]
    test=[held[1][k] for k in fields]
    # Fixed alpha=1 ridge and fold-local scaling; log1p(y) response and
    # reverse conversion to original experimental midazolam assay unit.
    y=[log1p(r[2]) for r in train]
    predicted=_ridge(xx,y,test)
    p=max(0.,expm1(min(30.,predicted)))
    if not isfinite(p):
        raise SourceNotQualified("Nonfinite donor-held-out predictive estimate")
    return p


def _paired_ci(deltas,seed):
    rng=Random(seed)
    n=len(deltas)
    draws=sorted(mean(deltas[rng.randrange(n)] for _ in range(n))
                 for _ in range(1999))
    return [round(draws[49],6),round(draws[1949],6)]


def evaluate_original_source(rows):
    donors,counts=source_contract(rows)
    errors={k:[] for k in ("training_mean",*ARMS)}
    for i,held in enumerate(donors):
        train=[r for j,r in enumerate(donors) if i!=j]
        for name in errors:
            p=_predict(train,held,name)
            errors[name].append(abs(p-held[2]))
    primary=[
        a-b for a,b in zip(errors["clinical_plus_RNA_and_protein"],
                           errors["clinical_plus_RNA"])
    ]
    secondary=[
        a-b for a,b in zip(errors["RNA_plus_protein"],
                           errors["RNA_only"])
    ]
    mae={name:round(mean(v),6) for name,v in errors.items()}
    return {
        "status":"LOCAL_ORIGINAL_DONOR_SOURCE_GATED_NOT_EXTERNAL_VALIDATION",
        "original_publisher_DOI":SOURCE_DOI,
        "n_unique_independent_original_liver_donors":74,
        "n_original_donors_by_histology":counts,
        "measured_drug_substrate":"midazolam CYP3A4 1hydroxylation",
        "outcome_units":"pmol/min/mg microsomal protein",
        "molecular_layers":["qPCR CYP3A4 mRNA/18S rRNA", "CYP3A4 protein targeted quantification"],
        "whole_transcriptome_assayed_in_this_model":False,
        "whole_proteome_assayed_in_this_model":False,
        "model":"74 donor-LODO fixed ridge alpha1, log1p original activity, train-only scaling",
        "mae_by_model_original_activity_units":mae,
        "primary_delta_MAE_clinical_plus_RNA_protein_minus_clinical_plus_RNA":round(mean(primary),6),
        "primary_descriptive_paired_loss_bootstrap_CI_95":_paired_ci(primary,SEED),
        "secondary_delta_MAE_RNA_protein_minus_RNA_only":round(mean(secondary),6),
        "secondary_descriptive_paired_loss_bootstrap_CI_95":_paired_ci(secondary,SEED+1),
        "n_original_people_combination_beats_clinical_plus_RNA":sum(x<0 for x in primary),
        "n_original_people_combination_worse_than_clinical_plus_RNA":sum(x>0 for x in primary),
        "not_measured_clinical_plasma_PK":True,
        "not_original_whole_omics_multiomics":True,
        "not_external_validation":True,
        "original_individual_donor_rows_and_predictions_exported":False,
    }


def analyze_authorized_source_csv(path):
    blob=Path(path).read_bytes()
    if not 1024<=len(blob)<=200000:
        raise SourceNotQualified("Source file inconsistent with 74 original rows")
    try:
        text=blob.decode("utf-8-sig")
    except UnicodeDecodeError:
        raise SourceNotQualified("Original source must be a UTF8 CSV") from None
    source=list(csv.DictReader(text.splitlines()))
    result=evaluate_original_source(source)
    result["original_authorized_local_curation_sha256"]=sha256(blob).hexdigest()
    result["not_asserted_to_be_publicly_available"]=True
    return result
