"""Unsupervised full-array 8-PC supplement to original HLC targeted CYP RNA.

Uses all 40,638 original processed array probe rows, but PCA/median
imputation/standardization are fitted exclusively within training sample
folds. The raw array source data MUST be supplied locally under legitimate
Synapse permissions. These unsupervised PCs have NO validated biological
gene/pathway annotations and may capture technical/liver-composition effects.
"""
from __future__ import annotations

import numpy as np

from omics2pk.hlc_real_cyp_validation import (
    load_original_hlc,OriginalHLCQualificationError,SEED
)

EVALUATED=("CYP1A2_activity","CYP2C9_activity","CYP2D6_activity",
           "CYP3A4M_activity","CYP3A4T_activity")


def evaluate_original_array_PCA8(clin,all_RNA,target_indices,y,*,seed=SEED):
    from sklearn.decomposition import PCA
    from sklearn.impute import SimpleImputer
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import Ridge
    from sklearn.model_selection import KFold
    from sklearn.pipeline import make_pipeline
    c=np.asarray(clin,float);x=np.asarray(all_RNA,dtype=np.float32);yy=np.asarray(y,float)
    if (len(c)!=len(x) or len(c)!=len(yy) or len(yy)<30
        or x.ndim!=2 or x.shape[1]<20 or
        not np.isfinite(c).all() or not np.isfinite(yy).all() or
        np.any(yy<0)):
        raise OriginalHLCQualificationError("Invalid original liver full-array endpoint")
    result={name:np.empty(len(yy)) for name in (
        "age_sex","age_sex_plus_target_CYP_RNA",
        "age_sex_plus_fullarray_trainonly_PCA8",
        "age_sex_plus_target_CYP_RNA_plus_fullarray_PCA8")}
    for train,test in KFold(5,shuffle=True,random_state=seed).split(yy):
        imp=SimpleImputer(strategy="median",keep_empty_features=True)
        z=imp.fit_transform(x[train])
        q=imp.transform(x[test])
        sc=StandardScaler()
        z=sc.fit_transform(z)
        q=sc.transform(q)
        pca=PCA(n_components=8,svd_solver="randomized",
                iterated_power=2,random_state=seed)
        pcs=pca.fit_transform(z)
        tq=pca.transform(q)
        candidates={
            "age_sex":(c[train],c[test]),
            "age_sex_plus_target_CYP_RNA":(
                np.column_stack((c[train],z[:,target_indices])),
                np.column_stack((c[test],q[:,target_indices]))),
            "age_sex_plus_fullarray_trainonly_PCA8":(
                np.column_stack((c[train],pcs)),
                np.column_stack((c[test],tq))),
            "age_sex_plus_target_CYP_RNA_plus_fullarray_PCA8":(
                np.column_stack((c[train],z[:,target_indices],pcs)),
                np.column_stack((c[test],q[:,target_indices],tq))),
        }
        for name,(xx,vv) in candidates.items():
            model=make_pipeline(StandardScaler(),Ridge(alpha=1.0))
            model.fit(xx,np.log1p(yy[train]))
            result[name][test]=np.maximum(0.,np.expm1(
                np.clip(model.predict(vv),-20.,20.)))
    return {
        "n_original_heldout_biological_liver_samples":len(yy),
        "n_original_unannotated_microarray_probes":x.shape[1],
        "n_train_only_PCA_components":8,
        "MAE_by_model_in_original_assay_units":
            {k:round(float(np.abs(v-yy).mean()),6) for k,v in result.items()},
        "no_original_individual_samples_or_measured_values_redistributed":True,
    }


def benchmark_original_hlc_fullarray(expression_zip,phenotype_zip):
    source=load_original_hlc(expression_zip,phenotype_zip)
    outputs={}
    for name in EVALUATED:
        y=source["outcomes"][name]
        valid=np.isfinite(y)&(y>=0)&np.isfinite(source["clinical"]).all(axis=1)
        gene="CYP3A4" if name.startswith("CYP3A4") else name.split("_")[0]
        outputs[name]=evaluate_original_array_PCA8(
            source["clinical"][valid],source["rna"][valid],
            source["gene_rows"][gene],y[valid])
    return {
        "status":"REAL_HLC_EXPLORATORY_WHOLE_TRANSCRIPTOME_VS_TARGET_CYP_5_PRESPECIFIED_E2_ASSAYS",
        "n_studied_original_activity_endpoints":len(outputs),
        "method":"Ridge alpha1, 5 held sample folds, all preprocessing/PCA8 fitted training-only",
        "original_source_only_in_vitro_human_liver_enzyme_activity":True,
        "not_clinical_in_vivo_pharmacokinetic_AUC":True,
        "original_assay_person_identifiers_and_expression_rows_not_exported":True,
        "results":outputs,
    }
