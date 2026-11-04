"""Source-verified Rotroff n32 clinical PK versus named/unknown metabolites/pathways.

No patient-level rows or predictions are returned. Original labels are bound by
the original publisher XLSX checksum and same-source matrix column index. Paper:
Rotroff et al 2016 155 chemically named/229 unknown of 384 GC-TOF features.
In vitro function, genotype, clinical PopPK are NOT observed in this experiment.
The published pathway hypotheses came from metformin exposure/PD associations;
they are NOT independent PK mechanism discovery.
"""
from __future__ import annotations
from collections import Counter
from random import Random
from statistics import mean
import re
import numpy as np
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from omics2pk.rotroff33_clinical_e3 import (
    actual_original_tables,measured_original_data,SEED,E3SourceError
)
from scripts.audit_rotroff33_metformin_original_pk_omics import _as_float
from scripts.qualify_rotroff33_original_metformin_patient_keys import SAMPLE
from scripts.audit_rotroff33_predose_feature_labels import UNRESOLVED

# Strict original exact chemical labels. These are partial *candidate* panels.
# No source-authenticated full graph or reaction enzyme proteins for this cohort.
HYPOTHESES={
    "urea_cycle_partial":("citrulline","ornithine","urea"),
    "purine_degradation_partial":("inosine","hypoxanthine"),
}
PLANS={
    "clinical_only":None,
    "clinical_plus_all_384":"all",
    "clinical_plus_named_155":"named",
    "clinical_plus_unknown_229":"unknown",
    "clinical_plus_partial_urea_cycle":"urea_cycle_partial",
    "clinical_plus_partial_purine":"purine_degradation_partial",
    "clinical_plus_two_candidate_pathways":"paired_paths",
}
EXPECTED_AUC_CLINICAL=2826.0472
EXPECTED_CLR_CLINICAL=115.1113
EXPECTED_AUC_ALL=2825.0240
EXPECTED_CLR_ALL=118.1449

def qualify_original_feature_labels(raw):
    observed,original_source=measured_original_data(raw)
    _,omics=actual_original_tables(raw)
    matrix=omics[0]
    predose={}
    for row in matrix[1:]:
        m=SAMPLE.fullmatch(row.get(0,"").strip())
        if m and m.group(3)=="A":
            if m.group(2) in predose:
                raise E3SourceError("Duplicate original human predose specimen")
            predose[m.group(2)]=row
    if len(predose)!=33:
        raise E3SourceError("Wrong original predose subject/assay count")
    cols=sorted(j for j in set().union(*(set(r) for r in predose.values())) if j>0
                and sum(_as_float(r.get(j)) is not None for r in predose.values())>=25)
    if len(cols)!=384:raise E3SourceError("Original predose 384 feature count changed")
    labels=[matrix[0].get(j,"").strip() for j in cols]
    if len(set(labels))!=384 or any(not x for x in labels):
        raise E3SourceError("Original feature annotation labels are missing or ambiguous")
    named=[i for i,x in enumerate(labels) if not UNRESOLVED.fullmatch(x)]
    unknown=[i for i,x in enumerate(labels) if UNRESOLVED.fullmatch(x)]
    if len(named)!=155 or len(unknown)!=229:
        raise E3SourceError("Original source 155/229 annotation gate mismatch")
    idx_by_name={x.lower():i for i,x in enumerate(labels)}
    paths={}
    for pathway,names in HYPOTHESES.items():
        if any(x.lower() not in idx_by_name for x in names):
            raise E3SourceError("Biochemical hypothesis is not fully present in source")
        paths[pathway]=[idx_by_name[x.lower()] for x in names]
    for row in observed:
        original=predose[row["internal_patient_id"]]
        if len(row["predose_omics"])!=384:
            raise E3SourceError("Original predose molecular profile size mismatches labels")
        for j,i in enumerate(cols):
            parsed=_as_float(original.get(i))
            seen=row["predose_omics"][j]
            if (parsed is None)!=(seen is None) or (parsed is not None and parsed!=seen):
                raise E3SourceError("Original molecular value/feature-label pairing failed")
    return observed,original_source,{
        "n_source_features":len(labels),
        "n_named_candidate_chemical_labels":len(named),
        "n_numeric_unknown_source_feature_identifiers":len(unknown),
        "n_original_subjects_with_matched_feature_headers":len(observed),
        "panel_source_label_counts":{k:len(v) for k,v in paths.items()},
        "named":named,"unknown":unknown,"paths":paths,
    }

def benchmark_one_target(records,target_key,annotation):
    selected=[row for row in records if row["targets"][target_key] is not None]
    if len(selected)!=32:raise E3SourceError("32 actual PK-linked original donors required")
    y=np.asarray([row["targets"][target_key] for row in selected],dtype=float)
    if np.any(~np.isfinite(y)) or np.min(y)<=0:raise E3SourceError("Invalid measured patient PK endpoint")
    logy=np.log1p(y)
    clinic=np.asarray([row["clinical"] for row in selected],dtype=float)
    molecular=np.asarray([[np.nan if x is None else x for x in row["predose_omics"]] for row in selected],dtype=float)
    if molecular.shape!=(32,384):raise E3SourceError("Original numerical predose matrix shape drift")
    names=annotation["named"];unknown=annotation["unknown"];paths=annotation["paths"]
    groups={
      "all":list(range(384)),
      "named":names,
      "unknown":unknown,
      **paths,
      "paired_paths":sorted(set(paths["urea_cycle_partial"]+paths["purine_degradation_partial"])),
    }
    pred={k:np.full(32,np.nan) for k in PLANS}
    for holdout in range(32):
        tr=np.arange(32)!=holdout
        xcl=StandardScaler().fit(clinic[tr])
        a=xcl.transform(clinic[tr]);b=xcl.transform(clinic[holdout:holdout+1])
        for plan,grp in PLANS.items():
            if grp is None:
                xtrain,xtest=a,b
            else:
                matrix=molecular[:,groups[grp]]
                imputer=SimpleImputer(strategy="median")
                train=imputer.fit_transform(matrix[tr])
                test=imputer.transform(matrix[holdout:holdout+1])
                std=StandardScaler().fit(train)
                ztrain=std.transform(train);ztest=std.transform(test)
                # Three PCs for full, named and unknown; one for each incomplete
                # biochemical hypothesis to reduce overfitting on n=32.
                components=1 if grp in ("urea_cycle_partial","purine_degradation_partial","paired_paths") else 3
                pca=PCA(n_components=components,svd_solver="randomized",random_state=SEED)
                mtrain=pca.fit_transform(ztrain);mtest=pca.transform(ztest)
                xtrain=np.column_stack((a,mtrain));xtest=np.column_stack((b,mtest))
            fitted=Ridge(alpha=1.0,fit_intercept=True).fit(xtrain,logy[tr])
            val=max(0.0,float(np.expm1(min(float(fitted.predict(xtest)[0]),30.0))))
            if not np.isfinite(val):raise E3SourceError("Failed finite original held-out prediction")
            pred[plan][holdout]=val
    loss={k:np.abs(p-y) for k,p in pred.items()}
    result={}
    for k,p in pred.items():
        diff=p-y
        result[k]={
          "MAE_original_table_numeric_units":round(float(np.mean(np.abs(diff))),6),
          "RMSE_original_table_numeric_units":round(float(np.sqrt(np.mean(diff**2))),6),
          "within_20pct_of_measured_PK_percent":round(float(np.mean(np.abs(diff)/y<=.20)*100),4)
        }
    original_a=EXPECTED_AUC_CLINICAL if "AUClast" in target_key else EXPECTED_CLR_CLINICAL
    original_all=EXPECTED_AUC_ALL if "AUClast" in target_key else EXPECTED_CLR_ALL
    if abs(result["clinical_only"]["MAE_original_table_numeric_units"]-original_a)>0.02 or abs(result["clinical_plus_all_384"]["MAE_original_table_numeric_units"]-original_all)>0.02:
        raise E3SourceError("Reference original human Rotroff validation failed replay parity")
    rng=Random(SEED)
    deltas={}
    for new in PLANS:
        if new=="clinical_only":continue
        d=loss[new]-loss["clinical_only"]
        reps=sorted(float(np.mean([d[rng.randrange(32)] for _ in range(32)])) for _ in range(1999))
        deltas[new]={
          "delta_MAE_positive_is_worse":round(float(np.mean(d)),6),
          "descriptive_paired_95pct_CI":[round(reps[49],6),round(reps[1949],6)],
          "n_subjects_lower_heldout_absolute_error":int(np.sum(d<0)),
        }
    return {"n_unique_original_PK_subjects":32,
            "evaluation":"internal leave-one-original-patient-out, no external validation",
            "performance":result,"descriptive_paired_contrasts":deltas}

def run_source_real_experiment(raw):
    rows,source,annotation=qualify_original_feature_labels(raw)
    return {
      "status":"SOURCE_VERIFIED_REAL_PREDOSE_HUMAN_METFORMIN_PK_PATHWAY_EXPLORATORY",
      "study_DOI":"10.3389/fphar.2016.00135",
      "source_omic_specimen":"human circulating plasma predose GC-TOF (not liver, no genotype)",
      "source_original_clinical_SHA256":source["original_clinical_workbook_SHA256"],
      "source_original_omics_SHA256":source["original_untargeted_metabolomic_workbook_SHA256"],
      "source_annotation_counts":{key:value for key,value in annotation.items() if key not in ("named","unknown","paths")},
      "endpoints":{k:benchmark_one_target(rows,k,annotation) for k in sorted(rows[0]["targets"])},
      "original_PK_numeric_units_unresolved":True,
      "prior_selected_pathway_hypotheses_not_independent":"urea/purine candidate neighborhoods were explored in published Rotroff metformin PD/exposure article; not a discovery-free clinical validation",
      "full_reaction_pathway_coverage_or_MSI_level1_confirmed":False,
      "complete_original_genotype_or_PopPK_baseline_available":False,
      "real_subject_level_molecular_or_PK_data_exported":False,
      "interpretation":"A negative/positive selected pathway CV is hypothesis generating, not a mechanistic demonstration; use matched clinical+genotype+PopPK plus external cohort for primary thesis goal."
    }
