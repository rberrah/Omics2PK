"""MaveDB Suiter-2020 NUDT15 dual phenotype: abundance versus thiopurine tolerance.

Source original paired score set urn:mavedb:00000055-0-1.
The outcome is independently measured HEK293T thioguanine cell survival
(normalized *tolerance*, low means sensitive), not NUDT15 hydrolase
activity, clinical patient toxicity, drug PK or patient response.
Cell-line constructs are indexed by amino-acid substitution, NOT donors.
Holdout is by original protein position to reduce substitution leakage.

This analysis uses all measured matched variants, then retains confirmed
single amino-acid missense substitutions for AA-level modeling. Other
variants (synonymous, nonsense, stop, complex) remain in source audit
but are not silently mapped to a synthetic single amino acid change.
"""
from __future__ import annotations
import csv
import io
import re
from collections import defaultdict
from math import isfinite, sqrt
from random import Random
from statistics import mean

from .public_cyp2c9_variant_function_benchmark import (
    _features, _train_ridge, _ranks, _correlation, _auroc,
)

URN="urn:mavedb:00000055-0-1"
SOURCE_DOI="10.1073/pnas.1915680117"
EXPECTED_ORIGINAL_PAIRED_ROWS=3015
EXPECTED_PAIRED_MISSENSE=2825
EXPECTED_SOURCE_SHA256="2e248907c0dcf54372d6ebf79c4d8a91b9ec2f2cd50d6c8aa365abec39d85a59"
AMINO_LONG={"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C",
"Gln":"Q","Glu":"E","Gly":"G","His":"H","Ile":"I",
"Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P",
"Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V"}
REQUIRED_HEADERS={"accession","hgvs_nt","hgvs_splice","hgvs_pro",
                  "score","abundance_score","sensitivity_score"}

class NUDT15PublicSourceError(ValueError):
    pass


def _numeric(value):
    if value in (None,"","NA","nan","NaN"):
        return None
    try:
        number=float(value)
    except (ValueError,TypeError,OverflowError) as exc:
        raise NUDT15PublicSourceError("Original MaveDB phenotype has nonnumeric nonmissing value") from exc
    if not isfinite(number):
        raise NUDT15PublicSourceError("Original MaveDB phenotype contains nonfinite observation")
    return number


def _variant(hgvs):
    """Accept only source amino-acid substitutions, never guess splice/LOF."""
    if not isinstance(hgvs,str):
        return None
    text=hgvs.strip()
    text=text.removeprefix("p.")
    if text.startswith("(") and text.endswith(")"):
        text=text[1:-1]
    match=re.fullmatch(r"([A-Z][a-z]{2})([1-9][0-9]*)([A-Z][a-z]{2})",text)
    if match:
        ref,alt=AMINO_LONG.get(match[1]),AMINO_LONG.get(match[3])
        pos=int(match[2])
    else:
        short=re.fullmatch(r"([ACDEFGHIKLMNPQRSTVWY])([1-9][0-9]*)([ACDEFGHIKLMNPQRSTVWY])",text)
        if not short:return None
        ref,pos,alt=short[1],int(short[2]),short[3]
    if not ref or not alt or ref==alt or not 1<=pos<=163:
        return None
    return ref,pos,alt


def parse_original_mavedb(content, *, strict_original=True):
    if not isinstance(content,(bytes,bytearray)) or len(content)>5_000_000:
        raise NUDT15PublicSourceError("Original NUDT15 paired CSV must be bounded original bytes")
    try:
        text=bytes(content).decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise NUDT15PublicSourceError("Original MaveDB NUDT15 CSV encoding mismatch") from exc
    reader=csv.DictReader(io.StringIO(text))
    if not reader.fieldnames or not REQUIRED_HEADERS.issubset(reader.fieldnames):
        raise NUDT15PublicSourceError("Original MaveDB paired score set missing original columns")
    if len(reader.fieldnames)!=len(set(reader.fieldnames)):
        raise NUDT15PublicSourceError("Original MaveDB repeated column name")
    seen=set()
    all_rows=0
    usable=[]
    not_single_aa=0
    missing_scores=0
    for row in reader:
        all_rows+=1
        if None in row:
            raise NUDT15PublicSourceError("MaveDB column count mismatch")
        urn=row["accession"]
        if not isinstance(urn,str) or not urn.startswith(URN+"#") or urn in seen:
            raise NUDT15PublicSourceError("Original MaveDB variant accession missing, duplicated or from other set")
        seen.add(urn)
        ab=_numeric(row["abundance_score"])
        se=_numeric(row["sensitivity_score"])
        if ab is None or se is None:
            missing_scores+=1
            continue
        change=_variant(row["hgvs_pro"])
        if change is None:
            not_single_aa+=1
            continue
        aa,pos,bb=change
        usable.append({
            "variant":aa+str(pos)+bb,
            "start":aa,"end":bb,"position":pos,
            "abundance":ab,
            "thiopurine_tolerance":se,
        })
    if len({r["variant"] for r in usable})!=len(usable):
        raise NUDT15PublicSourceError("Original MaveDB variant identities duplicate after protein mapping")
    if strict_original and (all_rows!=EXPECTED_ORIGINAL_PAIRED_ROWS
            or len(usable)!=EXPECTED_PAIRED_MISSENSE):
        raise NUDT15PublicSourceError(
            "ORIGINAL_NUDT15_PAIRED_COHORT_NOT_QUALIFIED:"
            f" source={all_rows} paired_protein_missense={len(usable)}"
            f" unparsed={not_single_aa} missing_both={missing_scores}")
    if len(usable)<60 or len({r["position"] for r in usable})<20:
        raise NUDT15PublicSourceError("Insufficient source protein residues and measured outcomes")
    return usable, {
        "n_source_MaveDB_paired_score_set_rows":all_rows,
        "n_paired_missense_variant_ids":len(usable),
        "n_original_nonmissense_or_unparsed_variants":not_single_aa,
        "n_rows_with_missing_any_scores":missing_scores,
        "original_source_row_data_exported":False,
    }


def _percentile_interval(samples):
    values=sorted(samples)
    return [round(values[24],6),round(values[974],6)]


def benchmark_nudt15_original(rows):
    if len({r["position"] for r in rows})<20:
        raise NUDT15PublicSourceError("Insufficient independent heldout protein positions")
    if len({r["variant"] for r in rows})!=len(rows):
        raise NUDT15PublicSourceError("Duplicate source variant labels")
    methods={
        "training_mean":lambda r:[],
        "simple_AA_properties_only":_features,
        "measured_protein_stability":lambda r:[r["abundance"]],
        "AA_properties_plus_measured_stability":lambda r:[*_features(r),r["abundance"]],
    }
    observed=[r["thiopurine_tolerance"] for r in rows]
    predictions={label:[None]*len(rows) for label in methods}
    folds=[]
    for fold in range(5):
        train=[r for r in rows if ((r["position"]*2654435761)&0xFFFFFFFF)%5!=fold]
        test=[(i,r) for i,r in enumerate(rows) if
              ((r["position"]*2654435761)&0xFFFFFFFF)%5==fold]
        if min(len(train),len(test))<10:
            raise NUDT15PublicSourceError("Position-disjoint train/test underpowered")
        if {r["position"] for r in train}&{r["position"] for _,r in test}:
            raise NUDT15PublicSourceError("Variant position data leakage")
        folds.append(len(set(r["position"] for _,r in test)))
        for label,feature in methods.items():
            estimator=_train_ridge(
                [{"activity":r["thiopurine_tolerance"],
                  "abundance":r["abundance"],
                  "start":r["start"],"end":r["end"]}
                 for r in train],
                feature,
                alpha=30.0)
            for idx,record in test:
                predictions[label][idx]=estimator(record)
    if any(not all(isinstance(q,(float,int)) and isfinite(q) for q in v)
           for v in predictions.values()):
        raise NUDT15PublicSourceError("Missing numerical original heldout tolerance score")
    # This is an illustrative predeclared source-phenotype binary split,
    # not a clinical toxicity threshold, nor a CPIC recommendation.
    binary=[int(y>=0.8) for y in observed]
    metrics={}
    for label,pred in predictions.items():
        d=[p-y for p,y in zip(pred,observed)]
        metrics[label]={
            "heldout_MAE_thiopurine_tolerance_score":round(mean(abs(x) for x in d),6),
            "heldout_RMSE_thiopurine_tolerance_score":round(sqrt(mean(x*x for x in d)),6),
            "heldout_Spearman":(round(_correlation(_ranks(pred),_ranks(observed)),6)
                                if _correlation(_ranks(pred),_ranks(observed)) is not None
                                else None),
            "illustrative_AUROC_tolerance_ge_0_8_not_clinical_toxicity":round(
                _auroc(pred,binary),6),
        }
    groups=defaultdict(list)
    for i,r in enumerate(rows):
        groups[r["position"]].append(i)
    positions=sorted(groups)
    rand=Random(20261105)
    paired={}
    for test,base in (
        ("measured_protein_stability","training_mean"),
        ("AA_properties_plus_measured_stability","measured_protein_stability"),
    ):
        a,b=predictions[test],predictions[base]
        delta=mean(abs(a[i]-observed[i])-abs(b[i]-observed[i])
                   for i in range(len(rows)))
        bootstrap=[]
        for _ in range(999):
            sample=[rand.choice(positions) for _ in positions]
            idx=[i for residue in sample for i in groups[residue]]
            bootstrap.append(mean(abs(a[i]-observed[i])-abs(b[i]-observed[i])
                                  for i in idx))
        paired[test+"_minus_"+base]={
            "observed_delta_MAE_negative_is_better":round(delta,6),
            "descriptive_position_cluster_95pct_bootstrap_no_refit":_percentile_interval(bootstrap),
        }
    # NUDT15 score of >0.8 means resistance/tolerance, <0.3 high
    # thiopurine sensitivity. These are arbitrary research descriptive
    # thresholds, not CPIC or fitted clinical cutoffs.
    stable_sensitive=sum(r["abundance"]>=0.8 and r["thiopurine_tolerance"]<0.3
                         for r in rows)
    unstable_tolerant=sum(r["abundance"]<0.3 and r["thiopurine_tolerance"]>=0.8
                          for r in rows)
    return {
        "status":"PUBLIC_NUDT15_VARIANT_MATCHED_CELL_LINE_PD_NOT_HUMAN_PATIENT_TOXICITY",
        "original_MaveDB_score_set":URN,
        "source_reference_DOI":SOURCE_DOI,
        "n_paired_missense_variants":len(rows),
        "n_independent_aa_positions":len(positions),
        "n_independent_aa_positions_per_heldout_fold":folds,
        "held_position_out_methods":metrics,
        "paired_position_bootstrap":paired,
        "n_stable_but_thiopurine_sensitive_research_thresholds_0_8_0_3":stable_sensitive,
        "n_unstable_but_thiopurine_tolerant_research_thresholds_0_3_0_8":unstable_tolerant,
        "exposure_outcome_direction":"higher_score_means_tolerant_to_thiopurine",
        "functional_measure_not_same_as_protein_stability":True,
        "no_pretreatment_mRNA_or_protein_decay_rate_assayed":True,
        "not_direct_NUDT15_enzymatic_substrate_metabolite_measure":True,
        "not_original_human_PK_or_clinical_toxicity":True,
        "original_patient_drug_dose_or_MPA_tacrolimus_not_predicted":True,
        "original_donor_rows_or_variant_scores_not_exported":True,
        "CV_grouped_by_protein_position_not_random_variant":True,
        "no_model_uses_PD_score_as_omic_predictor":True,
    }
