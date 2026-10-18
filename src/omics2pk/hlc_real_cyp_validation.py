"""Reproducible same-original-liver CYP RNA and genetic-source validation.

User-authorized original Sage Human Liver Cohort archive bytes must remain
private and never enter GitHub, CI, git-lfs or an issue attachment.

Original historical study: Schadt et al. PLoS Biology 2008 and Yang et al.
Genome Research 2010. Input: curatedExpressionLiver.zip and
curatedPhenotypesynapse.zip. Optional: curatedGenotype.zip. These contain
measured donor activity assays, NOT in-vivo plasma AUC/CLint.

The public code has NO participant rows, genotype values, expression data or
source subject identifiers. All model transforms are fitted on training folds.
"""
from __future__ import annotations

from collections import Counter
from io import TextIOWrapper
from itertools import chain
from math import isfinite
from pathlib import Path
import csv
import re
from zipfile import ZipFile

import numpy as np

SEED=20261018
ACTIVITIES=(
    "CYP1A2_activity","CYP2A6_activity","CYP2B6_activity",
    "CYP2C8_activity","CYP2C9_activity","CYP2C19_activity",
    "CYP2D6_activity","CYP2E1_activity","CYP3A4M_activity",
    "CYP3A4T_activity",
)
GENES=("CYP1A2","CYP2A6","CYP2B6","CYP2C8","CYP2C9","CYP2C19",
       "CYP2D6","CYP2E1","CYP3A4","CYP3A5","POR","ABCB1","SLC22A2",
       "SLC47A1","SLC47A2")
PHENOTYPE_HISTORICAL_FIELDS=("AGE_(YRS)","GENDER",*ACTIVITIES)
ASSAY_N=466
PHENOTYPE_SOURCE_N=836


class OriginalHLCQualificationError(ValueError):
    pass


def _read_archived_table(filename,member,**kwargs):
    import pandas as pd
    with ZipFile(filename) as z:
        if member not in z.namelist():
            raise OriginalHLCQualificationError("Expected original source archive member absent")
        with z.open(member) as handle:
            return pd.read_csv(handle,sep="\t",**kwargs)


def load_original_hlc(expression_zip,phenotype_zip):
    """Only in-memory source; no originals or human sample IDs emitted."""
    feature=_read_archived_table(expression_zip,"features.txt",dtype=str)
    expr=_read_archived_table(expression_zip,"expression.txt",
                              index_col=0,dtype={"feature_id":str})
    ph=_read_archived_table(phenotype_zip,"phenotype.txt",index_col=0,dtype=str)
    if not {"feature_id","genesymbol"}.issubset(feature.columns):
        raise OriginalHLCQualificationError("Original array annotation missing")
    if (len(feature)!=40638 or expr.shape!=(40638,ASSAY_N)
        or ph.shape[1]!=PHENOTYPE_SOURCE_N or
        not np.array_equal(expr.index.to_numpy(dtype=str),
                           feature["feature_id"].to_numpy(dtype=str))):
        raise OriginalHLCQualificationError("Source array or phenotype dimensions/probe IDs changed")
    ids=expr.columns.astype(str)
    if (len(set(ids))!=ASSAY_N or not set(ids).issubset(set(ph.columns.astype(str)))
        or any(not re.fullmatch(r"\d{7}",id) for id in ids)):
        raise OriginalHLCQualificationError("Original clinical source individual key crosswalk invalid")
    if any(name not in ph.index for name in PHENOTYPE_HISTORICAL_FIELDS):
        raise OriginalHLCQualificationError("Original CYP measured activity row missing")
    from pandas import to_numeric
    ph=ph.loc[:,ids]
    age=to_numeric(ph.loc["AGE_(YRS)"],errors="coerce").to_numpy(float)
    sex=ph.loc["GENDER"].map({"Male":1.,"Female":0.}).to_numpy(float)
    clin=np.column_stack([age,sex])
    symbols=feature["genesymbol"].fillna("").to_numpy(dtype=str)
    positions={g:np.where(symbols==g)[0] for g in GENES}
    if any(len(positions[g])==0 for g in GENES[:9]):
        raise OriginalHLCQualificationError("Required CYP probe annotation missing")
    exp=np.asarray(expr.to_numpy(dtype=np.float32).T)
    exp[~np.isfinite(exp)]=np.nan
    yy={key:to_numeric(ph.loc[key],errors="coerce").to_numpy(float)
        for key in ACTIVITIES}
    return {"sample_ids":ids.tolist(),"clinical":clin,
            "rna":exp,"gene_rows":positions,"outcomes":yy}


def _prediction_errors(x,y,folds):
    from sklearn.impute import SimpleImputer
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import Ridge
    from sklearn.pipeline import make_pipeline
    p=np.empty(len(y))
    for train,test in folds:
        m=make_pipeline(SimpleImputer(strategy="median"),StandardScaler(),
                        Ridge(alpha=1.0))
        m.fit(x[train],np.log1p(y[train]))
        p[test]=np.maximum(0.,np.expm1(np.clip(m.predict(x[test]),-20.,20.)))
    return np.abs(p-y)


def evaluate_fold_matched_cyp(clin,target,panel,measured,*,seed=SEED,
                              permutations=199):
    from sklearn.model_selection import KFold
    clinical=np.asarray(clin,float)
    rna=np.asarray(target,float)
    multi=np.asarray(panel,float)
    y=np.asarray(measured,float)
    if (clinical.ndim!=2 or clinical.shape[1]!=2 or
        rna.ndim!=2 or multi.ndim!=2 or
        len(y)<30 or not (len(y)==len(clinical)==len(rna)==len(multi))):
        raise OriginalHLCQualificationError("No valid matched same-human RNA and CYP target design")
    if (not np.isfinite(y).all() or np.any(y<0)
        or not np.isfinite(clinical).all()):
        raise OriginalHLCQualificationError("Missing or unquantifiable original assay")
    folds=list(KFold(5,shuffle=True,random_state=seed).split(y))
    designs={
        "age_sex":clinical,
        "target_CYP_RNA_only":rna,
        "age_sex_plus_target_CYP_RNA":np.column_stack((clinical,rna)),
        "age_sex_plus_targeted_CYP_POR_transporter_panel":
            np.column_stack((clinical,multi)),
    }
    errors={k:_prediction_errors(x,y,folds) for k,x in designs.items()}
    gain=errors["age_sex_plus_target_CYP_RNA"]-errors["age_sex"]
    rng=np.random.default_rng(seed+len(y))
    draws=np.mean(rng.choice(gain,size=(1999,len(gain)),replace=True),axis=1)
    interval=np.quantile(draws,[.0025,.9975])
    # Approximate nuisance-preserving RNA negative control. Original clinical
    # sex and age unchanged; shuffle the entire RNA probe vector within
    # sex×age quartile strata, maintaining donor-level multimodal structure.
    age_group=np.digitize(clinical[:,0],np.quantile(clinical[:,0],[.25,.50,.75]),right=True)
    strat=4*clinical[:,1]+age_group
    placebo=[]
    for _ in range(permutations):
        scrambled=np.arange(len(y))
        for group in np.unique(strat):
            i=np.flatnonzero(strat==group)
            if len(i)>1:scrambled[i]=rng.permutation(i)
        fake=np.column_stack((clinical,rna[scrambled]))
        placebo.append(float(np.mean(_prediction_errors(fake,y,folds)-errors["age_sex"])))
    observed=float(gain.mean())
    p=(1+sum(null<=observed for null in placebo))/(1+permutations)
    return {
        "n_original_human_livers":len(y),
        "n_targeted_CYP_probe_features":rna.shape[1],
        "n_targeted_multi_gene_probe_features":multi.shape[1],
        "MAE_by_model_original_source_units":
            {k:round(float(np.mean(v)),6) for k,v in errors.items()},
        "delta_MAE_target_RNA_minus_age_sex":round(observed,6),
        "descriptive_paired_99point5pct_loss_bootstrap_delta":
            [round(float(x),6) for x in interval],
        "sex_and_age_quartile_stratified_RNA_shuffles":permutations,
        "one_sided_empirical_permutation_p":round(p,6),
        "absolute_error_improved_original_humans":int((gain<0).sum()),
        "absolute_error_worsened_original_humans":int((gain>0).sum()),
        "no_individual_original_predictions_or_source_values_exported":True,
    }


def benchmark_original_hlc_expression(expression_zip,phenotype_zip,
                                      *,permutations=199):
    input=load_original_hlc(expression_zip,phenotype_zip)
    clin=input["clinical"];rna=input["rna"];pos=input["gene_rows"]
    multi=np.asarray(sorted(set(j for g in GENES for j in pos[g])),int)
    output={}
    for k,target in enumerate(ACTIVITIES):
        y=input["outcomes"][target]
        mask=np.isfinite(y)&(y>=0)&np.isfinite(clin).all(axis=1)
        if mask.sum()<300:
            raise OriginalHLCQualificationError("Insufficient original matched source CYP donors")
        gene="CYP3A4" if target.startswith("CYP3A4") else target.split("_")[0]
        outcome=evaluate_fold_matched_cyp(
            clin[mask],rna[mask][:,pos[gene]],rna[mask][:,multi],y[mask],
            seed=SEED,permutations=permutations)
        outcome["assayed_enzymatic_activity"]=target
        outcome["target_gene_from_original_platform_metadata"]=gene
        output[target]=outcome
    return {
        "status":"ORIGINAL_HLC_SAME_HUMAN_CYP_RNA_IN_VITRO_ACTIVITY_VALIDATION",
        "n_original_phenotype_source_person_keys":PHENOTYPE_SOURCE_N,
        "n_exactly_matched_expression_and_phenotype_original_person_keys":ASSAY_N,
        "n_original_microarray_features":40638,
        "n_independent_CYP_activity_criteria":len(output),
        "fivefold_source_donor_key_CV":True,
        "all_transforms_fold_local":True,
        "source_processed_normalization_may_precede_this_CV":True,
        "bonferroni_alpha_for_ten_criteria":0.005,
        "observed_metabolic_activity_is_not_in_vivo_drug_PK":True,
        "protein_quantification_is_not_present_in_original_archives":True,
        "no_source_person_identifiers_or_assay_observations_redistributed":True,
        "endpoints":output,
    }
