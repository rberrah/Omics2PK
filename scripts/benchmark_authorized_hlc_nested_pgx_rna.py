"""Private authorized original HLC: nested CYP3A4/5 RNA incremental to true source SNP.

Requires the user's *authorized local* original Sage/Synapse archives; NEVER
downloads/uploads/prints source donor identities or individual values.
This is in-vitro HLC liver CYP activities, NOT tacrolimus/clinical PK.
Genomics source has rs776746 CYP3A5 or rs1799853 CYP2C9; no full PGx.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
import numpy as np
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

from omics2pk.hlc_real_cyp_validation import load_original_hlc,OriginalHLCQualificationError
from omics2pk.hlc_real_genomics_validation import _load_cyp_genotypes,TARGETS


class HLCNestedOriginalError(ValueError):
    pass


def _fold_fit(x,y,train,test,alpha):
    """Original donor-group fitted only on train. Log1p activity, MAE raw."""
    if set(train)&set(test):
        raise HLCNestedOriginalError("Donor leakage in source folds")
    if alpha not in (1.,20.,100.):
        raise HLCNestedOriginalError("Nonpreset regularization")
    xx=np.asarray(x,dtype=float)
    tr=xx[train]
    te=xx[test]
    # Do not let an all-missing training probe silently choose a test-derived fill.
    allmissing=np.isnan(tr).all(axis=0)
    if any(allmissing):
        raise HLCNestedOriginalError("No molecular feature measured in inner train")
    fill=np.nanmedian(tr,axis=0)
    tr=np.where(np.isnan(tr),fill,tr)
    te=np.where(np.isnan(te),fill,te)
    if not (np.isfinite(tr).all() and np.isfinite(te).all()):
        raise HLCNestedOriginalError("Invalid source predictors after fold imputation")
    sc=StandardScaler()
    zz=sc.fit_transform(tr)
    zt=sc.transform(te)
    model=Ridge(alpha=alpha).fit(zz,np.log1p(y[train]))
    return np.maximum(0.,np.expm1(np.clip(model.predict(zt),-20,20)))


def _source_repeated_nested_cv(y,designs,*,seeds=(20261106,20261107,20261108)):
    yy=np.asarray(y,dtype=float)
    n=len(yy)
    if n<50 or n>500 or not np.isfinite(yy).all() or np.any(yy<0):
        raise HLCNestedOriginalError("Unqualified original human liver activity target")
    if any(np.asarray(x).shape[0]!=n for x in designs.values()):
        raise HLCNestedOriginalError("Same donor intersection mismatch")
    model_loss={k:np.zeros((len(seeds),n)) for k in designs}
    tuned={k:Counter() for k in designs}
    for repeat,seed in enumerate(seeds):
        outer=KFold(5,shuffle=True,random_state=seed)
        for train,test in outer.split(yy):
            inner_splits=list(KFold(4,shuffle=True,random_state=seed+123).split(train))
            for name,x in designs.items():
                best=[]
                for a in (1.,20.,100.):
                    errors=[]
                    for it,iv in inner_splits:
                        actual_train=train[it]
                        actual_val=train[iv]
                        preds=_fold_fit(x,yy,actual_train,actual_val,a)
                        errors.extend(abs(preds-yy[actual_val]))
                    best.append((float(np.mean(errors)), -a, a))
                _,_,chosen=min(best)
                tuned[name][str(chosen)]+=1
                pred=_fold_fit(x,yy,train,test,chosen)
                model_loss[name][repeat,test]=abs(pred-yy[test])
    # Average repeated CV errors per *original biological donor*, not each
    # repeat as an independent new observation.
    e={k:v.mean(axis=0) for k,v in model_loss.items()}
    if "clinical_PGx" not in e:
        raise HLCNestedOriginalError("No original PGx baseline")
    baseline=e["clinical_PGx"]
    rng=np.random.default_rng(20261106+n)
    summary={}
    for k,loss in e.items():
        out={
            "observed_human_liver_CYP_activity_MAE":round(float(loss.mean()),6),
            "nested_4fold_alpha_selected_counts_across_3x5_outer_splits":
                dict(sorted(tuned[k].items())),
            "n_same_original_genotyped_livers":n,
            "no_original_donor_values_returned":True,
        }
        if k!="clinical_PGx":
            delta=loss-baseline
            b=np.mean(delta[rng.integers(0,n,size=(1999,n))],axis=1)
            out["delta_MAE_vs_clinical_PGx_negative_better"]=round(float(delta.mean()),6)
            out["descriptive_95pct_paired_human_bootstrap_no_nested_refits"]=[
                round(float(v),6) for v in np.quantile(b,[.025,.975])]
        summary[k]=out
    return summary


def original_source_nested_hlc(expression,phenotype,genotype):
    source=load_original_hlc(expression,phenotype)
    ids,snps=_load_cyp_genotypes(genotype)
    # Same qualified original study intersection already checked in main branch.
    original_ids={key:i for i,key in enumerate(source["sample_ids"])}
    overlap=sum(x in original_ids for x in ids)
    if overlap!=178:
        raise HLCNestedOriginalError("Original expression-genotype joined donor set changed")
    out={}
    for target,rs,genes in TARGETS:
        ysource=source["outcomes"][target]
        candidates=[
            ("CYP_target",genes),
            ("CYP_target_POR_ABCB1",tuple(dict.fromkeys((*genes,"POR","ABCB1")))),
        ]
        probes={}
        for name,symbols in candidates:
            pos=sorted(set(int(p) for symbol in symbols for p in source["gene_rows"][symbol]))
            if len(pos)<1 or len(pos)>100:
                raise HLCNestedOriginalError("Source genotype-target annotation invalid")
            probes[name]=pos
        keep=[]
        for i,key in enumerate(ids):
            if key not in original_ids:continue
            ix=original_ids[key]
            y=ysource[ix]
            clin=source["clinical"][ix]
            if not np.isfinite(y) or y<0 or not np.isfinite(clin).all():continue
            dosage={"CC":0., "CT":1., "TT":2.,"NA":np.nan}.get(snps[rs][i])
            if dosage is None:raise HLCNestedOriginalError("Unexpected original dbSNP call")
            # NOTE T dosage is not declared as expressor allele on unknown strand.
            keep.append((ix,clin,dosage,float(y)))
        n=len(keep)
        if n not in (171,172):
            raise HLCNestedOriginalError("Original HLC PGx+RNA donor intersection changed")
        original_indices=[r[0] for r in keep]
        clinical=np.asarray([r[1] for r in keep],dtype=float)
        dosage=np.asarray([r[2] for r in keep],dtype=float)
        y=np.asarray([r[3] for r in keep],dtype=float)
        base=np.column_stack([clinical,dosage])
        designs={
            "clinical_PGx":base,
            "clinical_PGx_plus_target_RNA":np.column_stack([
                base,source["rna"][original_indices][:,probes["CYP_target"]]]),
            "clinical_PGx_plus_target_RNA_POR_ABCB1":np.column_stack([
                base,source["rna"][original_indices][:,probes["CYP_target_POR_ABCB1"]]]),
        }
        # Augmented models choose feature *groups*, not individual
        # SNP/omics by heldout results. A future nested group-selection
        # model can be added independently after fixed-arm auditing.
        res=_source_repeated_nested_cv(y,designs)
        out[target]={
            "status":"ACTUAL_HLC_NESTED_GENOTYPE_VS_RNA_IF_ORIGINAL_AUTHORIZED_ARCHIVES",
            "n_original_paired_donors":n,
            "tested_source_rs":"rs"+rs,
            "source_call_is_genomic_T_dosage_not_full_CYP_star_phenotype":True,
            "n_original_target_CYP_RNA_probes":len(probes["CYP_target"]),
            "n_original_target_POR_ABCB1_probes":len(probes["CYP_target_POR_ABCB1"]),
            "measured_endpoint":"human liver original in-vitro CYP microsomal activity",
            "nested_ridge_models":res,
        }
    return {
        "status":"ORIGINAL_HLC_LOCAL_AUTHORIZED_ONLY_NESTED_E2",
        "genotype_expression_original_shared_donors":overlap,
        "n_activity_endpoints":len(out),
        "molecular_source":"human liver expression microarray, not PBMC",
        "genetic_source":"actual 1 original source SNP per corresponding CYP",
        "all_outcome_and_genotype_models_fitted_donor_nested_only":True,
        "no_patient_clinical_AUC_CL":True,
        "no_raw_donor_rows_or_predictions_exported":True,
        "results":out,
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--expression",required=True)
    p.add_argument("--phenotype",required=True)
    p.add_argument("--genotype",required=True)
    args=p.parse_args()
    try:
        print(json.dumps(original_source_nested_hlc(
             args.expression,args.phenotype,args.genotype),sort_keys=True,indent=2))
    except (HLCNestedOriginalError,OriginalHLCQualificationError,OSError,ValueError) as exc:
        print(json.dumps({"status":"ABSTAIN_ORIGINAL_HLC_SOURCE_NOT_AVAILABLE_OR_UNQUALIFIED",
                          "reason_class":type(exc).__name__}))
        raise SystemExit(2)
