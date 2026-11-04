#!/usr/bin/env python3
"""NCI60 real original three-layer RNA/protein/metabolomics + XGBoost Tree SHAP.

Important: NCI60 doxorubicin GI50 = CELL-LINE drug response (PD), NOT human PK.
Same 47 line identities but omics from DIFFERENT source experiments/batches.
Use 3x5 cell-line-disjoint CV, fixed XGB params, TRAIN-only unsupervised feature
filters (300 RNA, 300 protein, 120 metabolite), heldout exact Tree SHAP.
Raw molecular vectors, patient/cell-line labels and predictions never logged.
"""
from __future__ import annotations
import argparse
from collections import Counter
import json
from pathlib import Path
import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from xgboost import XGBRegressor, DMatrix

from benchmark_public_nci60_twoomics import SEED,load_dataset,fit_view
from benchmark_original_nci60_threeomics_swath import (
    load_official_proteome,distinct_map,PREFIX_TO_CLASS,CAPS,PIN,DOI,
)

# Frozen prior to independent results. No use of XGB feature ranking to select
# which endpoint or model to report.
PARAMS=dict(n_estimators=90,max_depth=2,learning_rate=.05,
    min_child_weight=4,subsample=.8,colsample_bytree=.8,reg_lambda=10.,
    reg_alpha=0.,n_jobs=1,random_state=SEED,tree_method="hist",
    objective="reg:squarederror",verbosity=0)
MODALITIES=("rna","prot","met")
COMBOS={"tumor_only":(),
       "tumor_rna":("rna",),"tumor_protein":("prot",),
       "tumor_metabolome":("met",),
       "tumor_rna_protein":("rna","prot"),
       "tumor_rna_metabolome":("rna","met"),
       "tumor_protein_metabolome":("prot","met"),
       "tumor_threeomics":MODALITIES}
# Mechanisms selected from established doxorubicin literature, NOT from
# this source GI50 or SHAP. Expression/protein abundance != actual function.
PREDECLARED_DRUG_GENES=("TOP2A","TOP2B","ABCB1","ABCC1","ABCG2","CBR1","CBR3")


def predict_xgb_shap(a,b,y):
    model=XGBRegressor(**PARAMS)
    model.fit(a,y)
    expected=np.asarray(model.predict(b),dtype=float).reshape(-1)
    shap=np.asarray(model.get_booster().predict(
        DMatrix(b),pred_contribs=True,approx_contribs=False),dtype=float)
    if shap.shape!=(len(b),a.shape[1]+1) or not np.allclose(
            expected,shap.sum(axis=1),rtol=1e-5,atol=1e-5):
        raise ValueError("Exact Tree SHAP additive consistency failed")
    return expected,shap[:,:-1]


def fit_targeted_whitebox(a,b,y):
    if a.shape[1]<2:
        raise ValueError("No genuine original biological target gene features")
    imp=SimpleImputer(strategy="median",keep_empty_features=True)
    a=imp.fit_transform(a);b=imp.transform(b)
    sc=StandardScaler().fit(a)
    predictor=Ridge(alpha=10.0).fit(sc.transform(a),y)
    pred=np.asarray(predictor.predict(sc.transform(b)),dtype=float)
    if not np.isfinite(pred).all():raise ValueError("White box model invalid")
    return pred


def evaluate_original(ids,classes,drugresponse,rna_names,met_names,
                      rna,met,protein_names,protein_classes,protein,
                      protein_gene_symbols,repeats=3):
    base=distinct_map(ids,"original drug-response RNA+metabolome")
    pros=distinct_map(protein_names,"source SWATH identities")
    shared=[key for key in base if key in pros]
    if len(shared)!=47:
        raise ValueError("Original 47 matched cross-study NCI60 lines required")
    row=np.asarray([base[k] for k in shared]);prow=np.asarray([pros[k] for k in shared])
    lineage=np.asarray(classes)[row]
    if not all(lineage[i]==protein_classes[prow[i]] for i in range(47)):
        raise ValueError("Cancer lineage source conflict")
    y=np.asarray(drugresponse)[row]
    arrays={"rna":np.asarray(rna[row]),"prot":np.asarray(protein[prow]),
            "met":np.asarray(met[row])}
    columns={"rna":tuple(rna_names),"prot":tuple(protein_gene_symbols),
             "met":tuple(met_names)}
    for mode in MODALITIES:
        if arrays[mode].shape!=(47,len(columns[mode])):
            raise ValueError("Original source layer feature labels misaligned")
    target_idx={
        "rna":[i for i,name in enumerate(columns["rna"]) if name in PREDECLARED_DRUG_GENES],
        "prot":[i for i,name in enumerate(columns["prot"]) if name in PREDECLARED_DRUG_GENES],
    }
    # 3 repetitions of EXACT same original 47 biologically independent line IDs
    # in all 8 XGB arms + two white-box arms.
    arms=tuple(COMBOS)+("targeted_RNA_Ridge","targeted_RNA_protein_Ridge")
    predicted={name:np.full((repeats,47),np.nan) for name in arms}
    abs_contrib={m:np.zeros(arrays[m].shape[1]) for m in MODALITIES}
    freq={m:Counter() for m in MODALITIES}
    f_used={m:Counter() for m in MODALITIES}
    n_explained=0
    n_folds=repeats*5
    for rep in range(repeats):
        splitter=KFold(n_splits=5,shuffle=True,random_state=SEED+rep)
        for train,test in splitter.split(np.arange(47)):
            enc=OneHotEncoder(sparse_output=False,handle_unknown="ignore")
            xcl=enc.fit_transform(lineage[train,None])
            xct=enc.transform(lineage[test,None])
            modes={}
            for key in MODALITIES:
                # Original source-specific train-only median, variable
                # feature ranking and standardization. No supervised selection.
                pc_a,pc_b,a,b,selected=fit_view(
                    arrays[key][train],arrays[key][test],CAPS[key])
                if (a.shape[1]!=len(selected) or not np.isfinite(a).all()
                        or not np.isfinite(b).all()):
                    raise ValueError("Ill-defined omic feature selection")
                modes[key]=(a,b,selected)
                for j in selected:f_used[key][int(j)]+=1
            for arm,views in COMBOS.items():
                train_mat=np.column_stack([xcl]+[modes[v][0] for v in views])
                test_mat=np.column_stack([xct]+[modes[v][1] for v in views])
                pred,shap=predict_xgb_shap(train_mat,test_mat,y[train])
                predicted[arm][rep,test]=pred
                if arm=="tumor_threeomics":
                    n_explained+=len(test)
                    off=xcl.shape[1]
                    for key in MODALITIES:
                        ix=modes[key][2]
                        contributions=np.abs(shap[:,off:off+len(ix)])
                        if contributions.shape[1]!=len(ix):
                            raise ValueError("Attribution/source view alignment mismatch")
                        for col,source_ix in enumerate(ix):
                            abs_contrib[key][int(source_ix)]+=float(contributions[:,col].sum())
                            if np.any(contributions[:,col]>1e-8):
                                freq[key][int(source_ix)]+=1
                        off+=len(ix)
            # Genuine predeclared RNA/protein doxorubicin targets, without
            # SHAP-based selection on the same validation results.
            for label,keys in (
                ("targeted_RNA_Ridge",("rna",)),
                ("targeted_RNA_protein_Ridge",("rna","prot"))):
                tr=[xcl];te=[xct]
                for key in keys:
                    ii=target_idx[key]
                    if ii:
                        tr.append(arrays[key][train][:,ii])
                        te.append(arrays[key][test][:,ii])
                predicted[label][rep,test]=fit_targeted_whitebox(
                    np.column_stack(tr),np.column_stack(te),y[train])
    if n_explained!=47*repeats:
        raise ValueError("Wrong heldout independent biological response coverage")
    if any(not np.isfinite(v).all() for v in predicted.values()):
        raise ValueError("Missing heldout true drug response predictions")
    loss={name:np.abs(v-y[None,:]).mean(axis=0) for name,v in predicted.items()}
    summaries={}
    for name,err in loss.items():
        summaries[name]={"heldout_MAE_original_GI50_score":round(float(err.mean()),6),
                         "delta_MAE_vs_tumor_baseline":round(float((err-loss["tumor_only"]).mean()),6)}
    contrasts={}
    rng=np.random.default_rng(SEED)
    boot=rng.integers(0,47,size=(1999,47))
    for name in summaries:
        if name=="tumor_only":continue
        d=loss[name]-loss["tumor_only"]
        low,high=np.quantile(d[boot].mean(axis=1),[.025,.975])
        contrasts[name]={"delta_MAE_negative_is_better":round(float(d.mean()),6),
                         "descriptive_paired_cell_line_95pct_interval_no_refit":[round(float(low),6),round(float(high),6)]}
    importance={}
    totals={key:float(abs_contrib[key].sum()) for key in MODALITIES}
    common_total=sum(totals.values())
    for key in MODALITIES:
        ranked=np.argsort(-abs_contrib[key],kind="stable")
        result=[]
        for j in ranked[:15]:
            if abs_contrib[key][j]<=0:continue
            # SwissProt source protein gene symbols may be non-unique.
            # This field is a GENERIC label, not a unique protein accession.
            result.append({"original_source_name":columns[key][int(j)],
                           "n_train_folds_where_selected_by_unsupervised_variance":int(f_used[key].get(int(j),0)),
                           "n_outer_folds_with_nonzero_test_SHAP":int(freq[key].get(int(j),0)),
                           "mean_abs_heldout_SHAP_GI50_score":round(float(abs_contrib[key][int(j)]/n_explained),7),
                           "identity_grade":("source_gene_symbol_not_isoform_specific" if key in ("rna","prot")
                                             else "source_metabolite_name_MSI_unverified"),
                           "biological_function_inferred":False})
        importance[key]={"top15_model_attributions":result,
                          "fraction_of_absolute_test_SHAP_across_all_three_omic_layers":round(float(totals[key]/common_total),5) if common_total>0 else None,
                          "n_original_source_features":len(columns[key]),
                          "n_max_features_per_fold":CAPS[key]}
    return {
        "status":"ORIGINAL_NCI60_MULTIVIEW_XGB_TREE_SHAP_REAL_CELLLINE_GI50_PD_NOT_PATIENT_PK",
        "n_independent_cell_lines":47,
        "n_repeat_correlated_folds":repeats,
        "n_same_subject_test_records":n_explained,
        "n_actual_source_RNA_protein_metabolite_layer_features":{
            key:len(columns[key]) for key in MODALITIES
        },
        "source_protein_original_SHA256":PIN,
        "original_protein_DOI":DOI,
        "source_omic_assays_different_studies_and_biological_batches":True,
        "actual_endpoint":"original measured in-vitro doxorubicin GI50 score, NOT human drug PK",
        "models_on_identical_original_cell_lines":summaries,
        "comparisons_vs_tumor_baseline":contrasts,
        "heldout_tree_SHAP_model_attribution":importance,
        "predeclared_white_box_gene_coverage":{
            key:[columns[key][i] for i in indices] for key,indices in target_idx.items()},
        "predeclared_doxorubicin_gene_candidates":PREDECLARED_DRUG_GENES,
        "model_XGB_parameters":PARAMS,
        "feature_selection":"median imputation, unsupervised variance top RNA300/Proteome300/Metabolome120, TRAIN-ONLY per fold; no supervised selection",
        "SHAP_studies_predictive_association_not_transport_or_doxorubicin_metabolism":True,
        "individual_cell_line_protein_RNA_metabolite_or_response_values_exported":False,
        "clinical_pk_claimed":False,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--data",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    ids,classes,y,rnames,mnames,rna,met=load_dataset(args.data)
    pn,pc,prot,sha,genes=load_official_proteome(return_gene_symbols=True)
    result=evaluate_original(ids,classes,y,rnames,mnames,rna,met,pn,pc,prot,genes)
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf8")
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=="__main__":main()
