"""Fixed XGBoost comparison on original HLC CYP metabolism (source private).

Same original n387–397 donor IDs, same KFold(5, seed=20261018) as Ridge,
same measured in-vitro CYP probe activity. XGB is intentionally FIXED,
not tuned on outer heldout outcomes. This does NOT constitute nested XGB
hyperparameter-optimized clinical model performance or external replication.
"""
from __future__ import annotations

import numpy as np

from omics2pk.hlc_real_cyp_validation import (
    load_original_hlc,SEED,ACTIVITIES,GENES,
    OriginalHLCQualificationError
)


def fit_fixed_source_xgb(design,y,*,seed=SEED):
    from sklearn.model_selection import KFold
    from xgboost import XGBRegressor
    x=np.asarray(design,float)
    yy=np.asarray(y,float)
    if (len(x)!=len(yy) or len(yy)<30
        or not np.isfinite(yy).all() or np.any(yy<0)):
        raise OriginalHLCQualificationError("Original CYP donor/outcome data invalid")
    pred=np.empty(len(yy))
    for train,test in KFold(5,shuffle=True,random_state=seed).split(yy):
        m=XGBRegressor(objective="reg:squarederror",
                       n_estimators=150,max_depth=2,learning_rate=.035,
                       min_child_weight=8,subsample=.85,colsample_bytree=.9,
                       reg_lambda=10.,tree_method="hist",
                       random_state=seed,n_jobs=1,verbosity=0)
        m.fit(x[train],np.log1p(yy[train]))
        pred[test]=np.maximum(0.,np.expm1(
            np.clip(m.predict(x[test]),-20.,20.)))
    return round(float(np.mean(np.abs(pred-yy))),6)


def benchmark_original_HLC_fixed_xgboost(expression_zip,phenotype_zip):
    source=load_original_hlc(expression_zip,phenotype_zip)
    c=source["clinical"];rna=source["rna"];genes=source["gene_rows"]
    panel=np.asarray(sorted({int(k) for g in GENES for k in genes[g]}))
    result={}
    for name in ACTIVITIES:
        y=source["outcomes"][name]
        good=np.isfinite(y)&(y>=0)&np.isfinite(c).all(axis=1)
        gene="CYP3A4" if name.startswith("CYP3A4") else name.split("_")[0]
        clinical=c[good]
        target=rna[good][:,genes[gene]]
        multi=rna[good][:,panel]
        candidates={
            "clinical_age_sex":clinical,
            "clinical_plus_target_CYP_expression":
                np.column_stack((clinical,target)),
            "clinical_plus_CYP_POR_transporter_panel":
                np.column_stack((clinical,multi)),
        }
        result[name]={
            "n_original_same_biological_donor_heldout_CYP_assays":int(good.sum()),
            "MAE_source_activity_units":
                {k:fit_fixed_source_xgb(v,y[good]) for k,v in candidates.items()},
            "original_source_individual_donor_data_exported":False,
        }
    return {
        "status":"ORIGINAL_HLC_MEASURED_CYP_IN_VITRO_FUNCTION_FIXED_NONLINEAR_XGB_EXPLORATORY",
        "fixed_XGB":"150 trees, depth2, lr.035, min_child_weight8, subsample.85, colsample.9, lambda10",
        "comparison_same_CV_as_Ridge":True,
        "source_drug_in_vitro_CYP_assay_not_human_plasma_PK":True,
        "source_whole_proteome_or_CYP_protein_phenotype_not_measured":True,
        "original_individuals_CYP_or_transcript_values_not_exported":False,
        "results":result,
    }
