#!/usr/bin/env python3
"""NCI60 source real RNA/protein/met metabolome: nested minimal modality budgets.

True 47 original cell lines, 3x5 outer source-disjoint folds, inner 4-fold.
Sources from distinct laboratories and batches; this is doxorubicin GI50
CELLLINE DRUG RESPONSE (PD), never an actual human PK benchmark.
All 8 candidate omics combos XGB-hist, train-fold-only feature filters.
The outer test cell line is NEVER consulted to choose 80/90% target arm.
"""
from __future__ import annotations
import argparse,json
from collections import Counter
from pathlib import Path
import numpy as np
from sklearn.model_selection import KFold
from sklearn.preprocessing import OneHotEncoder
from xgboost import XGBRegressor
from benchmark_public_nci60_twoomics import load_dataset,fit_view,SEED
from benchmark_original_nci60_threeomics_swath import (
    load_official_proteome,distinct_map,CAPS,PIN
)
from benchmark_original_nci60_threeomics_xgb_shap import PARAMS,COMBOS
from omics2pk.modalities_information_budget import choose_minimum_modalities

def _fit_all(assays,classes,y,train,test):
    enc=OneHotEncoder(handle_unknown="ignore",sparse_output=False)
    tr_class=enc.fit_transform(classes[train,None])
    te_class=enc.transform(classes[test,None])
    per_mode={}
    for mode,array in assays.items():
        _,_,tr_m,te_m,_=fit_view(array[train],array[test],CAPS[mode])
        per_mode[mode]=(tr_m,te_m)
    output={}
    for arm,subset in COMBOS.items():
        a=np.column_stack([tr_class]+[per_mode[v][0] for v in subset])
        b=np.column_stack([te_class]+[per_mode[v][1] for v in subset])
        clf=XGBRegressor(**PARAMS).fit(a,y[train])
        pred=np.asarray(clf.predict(b),dtype=float)
        if pred.shape!=(len(test),) or not np.isfinite(pred).all():
            raise ValueError("Incomplete original source cell-line heldout")
        output[arm]=pred
    return output

def _ci(a,b,y,draws):
    delta=(np.abs(a-y[None,:])-np.abs(b-y[None,:])).mean(axis=0)
    lo,hi=np.quantile(delta[draws].mean(axis=1),[.025,.975])
    return {"delta_MAE_negative_is_better":round(float(delta.mean()),6),
            "paired_descriptive_cell_line_95pct_no_refit":[round(float(lo),6),round(float(hi),6)],
            "n_original_independent_cell_lines":len(y)}

def run_original(ids,classes,y,rna,met,protein_names,protein_classes,protein):
    base=distinct_map(ids,"IntLIM doxorubicin/RNA/met original keys")
    prot=distinct_map(protein_names,"Guo2019 SWATH original keys")
    shared=[k for k in base if k in prot]
    if len(shared)!=47 or len(set(shared))!=47:
        raise ValueError("Same original 47 NCI60 exact biological identities required")
    a=np.asarray([base[k] for k in shared],dtype=int)
    b=np.asarray([prot[k] for k in shared],dtype=int)
    tissue=np.asarray(classes)[a]
    if any(tissue[i]!=protein_classes[b[i]] for i in range(47)):
        raise ValueError("Cross-source cancer lineage identity drift")
    resp=np.asarray(y)[a]
    original={"rna":np.asarray(rna[a]),"met":np.asarray(met[a]),
              "prot":np.asarray(protein[b])}
    if any(array.shape[0]!=47 for array in original.values()):
        raise ValueError("Original 47 multiomics arrays not qualified")
    arms=tuple(COMBOS)
    preds={arm:np.full((3,47),np.nan) for arm in arms}
    adaptive={q:np.full((3,47),np.nan) for q in (.8,.9)}
    choices={q:Counter() for q in (.8,.9)}
    observed_best_frac={q:[] for q in (.8,.9)}
    for repeat in range(3):
        outer=KFold(n_splits=5,shuffle=True,random_state=SEED+repeat)
        for fold,(tr,te) in enumerate(outer.split(np.arange(47))):
            # Subdivide ONLY 37/38 original line identities available to training.
            # Inner fold-specific imputation, variance selection and XGB.
            inner=KFold(n_splits=4,shuffle=True,random_state=SEED+100*repeat+fold)
            inner_p={arm:np.full(len(tr),np.nan) for arm in arms}
            for subtrain,subvalid in inner.split(tr):
                ix=tr[subtrain];iv=tr[subvalid]
                p=_fit_all(original,tissue,resp,ix,iv)
                for arm in arms:inner_p[arm][subvalid]=p[arm]
            if any(not np.isfinite(x).all() for x in inner_p.values()):
                raise ValueError("Missing internal CV source line outcome")
            inner_losses={arm:float(np.mean(np.abs(v-resp[tr])))
                          for arm,v in inner_p.items()}
            selected={q:choose_minimum_modalities(inner_losses,COMBOS,fraction=q,
                                                  reference="tumor_only")
                      for q in (.8,.9)}
            # Exactly one model fit on outer 37/38 train and same heldout
            # original 9/10 lines per fold. Selection is sealed already.
            p=_fit_all(original,tissue,resp,tr,te)
            for arm in arms:preds[arm][repeat,te]=p[arm]
            for fraction,choice in selected.items():
                adaptive[fraction][repeat,te]=p[choice["arm"]]
                choices[fraction][choice["arm"]]+=1
                if choice["best_inner_gain"]>0:
                    observed_best_frac[fraction].append(
                        choice["selected_arm_inner_gain"]/choice["best_inner_gain"])
    for vals in list(preds.values())+list(adaptive.values()):
        if vals.shape!=(3,47) or not np.isfinite(vals).all():
            raise ValueError("Missing outer-validated original doxorubicin GI50")
    error={arm:np.abs(pred-resp[None,:]).mean(axis=0)
           for arm,pred in preds.items()}
    stats={arm:{"mae_original_GI50_units":round(float(err.mean()),6),
                "n_modalities":len(COMBOS[arm]),
                "modalities":list(COMBOS[arm]),
                "selected_max_source_features_per_fold":sum(CAPS[m] for m in COMBOS[arm])}
           for arm,err in error.items()}
    boot=np.random.default_rng(SEED).integers(0,47,size=(1999,47))
    adaptive_report={}
    for fraction,values in adaptive.items():
        d=(np.abs(values-resp[None,:])).mean(axis=0)
        adaptive_report[str(int(fraction*100))]={
          "mae_original_GI50_units":round(float(d.mean()),6),
          "n_outer_folds":15,"n_biologically_independent_source_lines":47,
          "nested_gain_fraction_predeclared":fraction,
          "selected_arm_by_outer_fold":dict(choices[fraction]),
          "fraction_outer_folds_choosing_no_omics":round(choices[fraction].get("tumor_only",0)/15,5),
          "mean_n_omic_assays_selected":round(
            sum(len(COMBOS[k])*v for k,v in choices[fraction].items())/15,6),
          "relative_to_tumor_only":_ci(values,preds["tumor_only"],resp,boot),
          "mean_inner_gain_fraction_achieved_when_positive":
            round(float(np.mean(observed_best_frac[fraction])),5) if observed_best_frac[fraction] else None,
          "no_original_outer_source_cell_line_used_for_choice":True,
        }
    return {
      "status":"REAL_NCI60_SOURCE_47_CELLLINE_NESTED_80_90_MULTIMODAL_GAIN_PD_NOT_PK",
      "n_independent_original_lines":47,
      "n_correlated_repetition_test_records":141,
      "three_repetitions_five_outer_folds_four_inner_folds":True,
      "source_omics_from_noncontemporaneous_preparations":True,
      "source_proteome_original_SHA256":PIN,
      "response":"doxorubicin original GI50 in vitro drug sensitivity; NOT human plasma PK",
      "all_mono_pair_threeomic_fixed_heldout_XGB_models":stats,
      "adaptive_80_90pct_inner_heldout_gain_assay_budget":adaptive_report,
      "posthoc_descriptive_full_3omics_vs_each_source_comparator":{
          arm:_ci(preds["tumor_threeomics"],p,resp,boot)
          for arm,p in preds.items() if arm!="tumor_threeomics"
      },
      "no_patient_or_original_cell_line_level_values_exported":True,
      "external_independent_new_biology_validation":False,
      "clinical_PK_assay_evidence_claimed":False,
    }

def main():
    arg=argparse.ArgumentParser()
    arg.add_argument("--data",type=Path,required=True)
    arg.add_argument("--out",type=Path,required=True)
    opts=arg.parse_args()
    ids,classes,y,_,_,rna,met=load_dataset(opts.data)
    pn,pc,p,sha=load_official_proteome()
    result=run_original(ids,classes,y,rna,met,pn,pc,p)
    opts.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf8")
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=="__main__":main()
