"""Fold-safe mono-/multiomic ablations against *observed* pharmacokinetics.

Research-only. One genuine biological individual per row and per modality,
real same-person PK outcome, authentic PGx or explicit PGx missingness.
Every transform and PK estimator is fitted ONLY inside each outer training
partition. This is a deliberately conservative exploratory benchmark,
NOT a nonlinear mixed-effects PopPK model or proof of clinical superiority.
Do not print patient records, identifiers, or individual predictions.
"""
from __future__ import annotations

from datetime import datetime
from hashlib import sha256
from math import isfinite
from typing import Mapping

import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold, LeaveOneGroupOut
from sklearn.preprocessing import StandardScaler


class MultimodalPKBenchmarkError(ValueError):
    pass


def _assert(condition, error):
    if not condition:
        raise MultimodalPKBenchmarkError(error)


def _numeric(v, *, rows=None):
    try:
        x = np.asarray(v, dtype=float)
    except (ValueError, TypeError, OverflowError):
        raise MultimodalPKBenchmarkError("Invalid numerical source observations") from None
    _assert(x.ndim == 2 and np.isfinite(x).all() and
            (rows is None or x.shape[0] == rows),
            "Nonfinite or unmatched original measurements")
    return x


def _date(v):
    _assert(isinstance(v, str) and bool(v), "Absolute source timestamp missing")
    try:
        d = datetime.fromisoformat(v.replace("Z", "+00:00"))
    except ValueError:
        raise MultimodalPKBenchmarkError("Invalid ISO molecular/index-dose time") from None
    _assert(d.tzinfo is not None and d.utcoffset() is not None,
            "Original timestamp must include timezone, not just T0")
    return d


def _fit_ridge(Xtrain, Xtest, ytrain):
    scaler = StandardScaler().fit(Xtrain)
    fit = Ridge(alpha=1.0).fit(scaler.transform(Xtrain), np.log(ytrain))
    predict = np.exp(fit.predict(scaler.transform(Xtest)))
    _assert(np.isfinite(predict).all() and np.all(predict > 0),
            "Invalid original-PK research prediction")
    return predict


def _two_frozen_pcs(train, test):
    sc = StandardScaler().fit(train)
    ztr, zte = sc.transform(train), sc.transform(test)
    n = min(2, ztr.shape[0]-1, ztr.shape[1])
    _assert(n >= 1, "Not enough train people for molecular factors")
    p = PCA(n_components=n, svd_solver="randomized", random_state=20261019)
    return p.fit_transform(ztr), p.transform(zte)


def _source_scale(values, value_type):
    """Per-individual MS1 normalisation only, or source normalized omics.

    No information from *other* patients is used at this stage. Upstream
    batch correction is NOT certified by this operation.
    """
    typ = str(value_type).lower()
    if "raw" in typ and "count" in typ:
        raise MultimodalPKBenchmarkError("Unnormalized RNA counts cannot enter PK features")
    if any(v in typ for v in ("peak_area", "ms1_area", "ion_intensity")):
        _assert(np.all(values >= 0), "MS1 abundance cannot be negative")
        area = values.sum(axis=1, keepdims=True)
        _assert(np.all(area > 0), "Empty per-person MS1 peak-area spectrum")
        return np.log1p(values / area * 1_000_000), "per_person_TIC_log1p_research_only"
    _assert(any(k in typ for k in (
        "normalized", "vst", "rlog", "log2", "log-transformed",
        "protein_abundance", "tpm", "fpkm", "intensity_normalized"
    )), "Molecular input must declare an acceptable measurement scale")
    return values, "upstream_source_scale_not_independently_verified"


def compare_modalities_pk(
    *, subject_ids, modalities: Mapping[str, Mapping[str, object]], clinical,
    measured_pk, index_dose_at, pk_target_name, pk_target_unit,
    source_kind, source_sha256, genotype_dosages=None, fixed_scores=None,
    n_outer_splits=5, max_lookback_days=30,
    include_xgboost=False, validation_mode="participant_kfold", site_ids=None,
):
    """Compare single-omic, combined and optionally fixed biological scores.

    modalities[omic] requires keys:
      feature_ids: unique source IDs (3..10000)
      values: person-by-feature original numeric observations, exact subject order
      collected_at: timezone-aware per-person original sampling ISO timestamps
      tissue: nonempty original specimen tissue/compartment
      value_type: explicit raw MS1 peak area or normalized molecular scale

    fixed_scores maps omic -> original feature_id -> preregistered SIGNED
    weight; if these were selected using the same cohort PK, they are NOT
    qualified as preregistered biological weights.

    Return aggregate only, no participant IDs or predictions. Models are
    descriptive research (log outcome Ridge alpha 1). Baseline genotype
    comparator uses observed 0/1/2 allele dosages if and ONLY if available.
    """
    ids=tuple(subject_ids)
    n=len(ids)
    _assert(12 <= n <= 400 and len(set(ids)) == n and all(
        isinstance(x,str) and x.strip() for x in ids),
        "Need 12-400 unique genuine biological subject identifiers")
    _assert(source_kind in ("synthetic","original_source_declared"),
            "Actual-versus-fictional data classification required")
    _assert(isinstance(source_sha256,str) and len(source_sha256)==64
            and all(k in "0123456789abcdef" for k in source_sha256),
            "Exact original source SHA-256 declaration required")
    _assert(isinstance(pk_target_name,str) and pk_target_name.strip()
            and isinstance(pk_target_unit,str) and pk_target_unit.strip(),
            "Observed PK endpoint name and physical units required")
    _assert(isinstance(max_lookback_days,(float,int)) and 0 < max_lookback_days <= 365,
            "Invalid prespecified sampling lookback")
    _assert(2 <= n_outer_splits <= min(10,n), "Invalid fixed outer patient folds")
    _assert(validation_mode in ("participant_kfold", "leave_one_site_out"),
            "Explicit participant-vs-site validation mode required")
    if validation_mode == "participant_kfold":
        _assert(site_ids is None,
                "Site labels supplied: use leave_one_site_out; no silent center mixing")
    else:
        _assert(site_ids is not None and len(site_ids) == n,
                "Site-disjoint validation requires a known site for every original person")
        _assert(all(isinstance(site,str) and site.strip() for site in site_ids),
                "Missing or empty original site provenance is not qualified")
    _assert(isinstance(modalities,Mapping) and 1 <= len(modalities) <= 3
            and set(modalities) <= {"transcriptomics","proteomics","metabolomics"},
            "One to three recognized individual omic modalities required")
    _assert(len(index_dose_at)==n, "PK index dose time must match every person")

    C=_numeric(clinical, rows=n)
    _assert(1 <= C.shape[1] <= 50,"Unmatched or excessive clinical covariates")
    y=np.asarray(measured_pk, dtype=float)
    _assert(y.shape==(n,) and np.isfinite(y).all() and np.all(y>0),
            "Must have genuine positive measured individual PK outcomes")
    G=None
    if genotype_dosages is not None:
        G=_numeric(genotype_dosages,rows=n)
        _assert(1 <= G.shape[1] <= 20 and np.all(G == np.floor(G))
                and np.all((G>=0)&(G<=2)),
                "Authentic 0/1/2 genotype dosages required; no imputation")
    baseline="clinical_pgx" if G is not None else "clinical"
    B=C if G is None else np.column_stack((C,G))

    # Stable whole-person folds independent of original per-assay row order.
    canonical=np.argsort(np.asarray(ids),kind="stable")
    C=C[canonical]; B=B[canonical]; y=y[canonical]
    groups=None
    outer_folds=n_outer_splits
    if validation_mode == "leave_one_site_out":
        groups=np.asarray(site_ids,dtype=str)[canonical]
        original_sites=set(groups.tolist())
        outer_folds=len(original_sites)
        _assert(2 <= outer_folds <= 10,
                "Need 2-10 distinct original source sites for whole-center validation")
        counts={site:int(np.count_nonzero(groups == site)) for site in original_sites}
        _assert(n-max(counts.values())>=8,
                "Each held-out site must leave >=8 actual training people")
    dose_dates=[_date(index_dose_at[i]) for i in canonical]
    samples={}
    score_weights={}
    source_scales={}
    fixed_scores = {} if fixed_scores is None else fixed_scores
    _assert(isinstance(fixed_scores,Mapping) and set(fixed_scores)<=set(modalities),
            "Only externally prespecified modality score definitions accepted")
    for omic, spec in sorted(modalities.items()):
        _assert(isinstance(spec,Mapping), "Molecular assay record is not a mapping")
        features=tuple(spec.get("feature_ids",()))
        _assert(3<=len(features)<=10000 and len(set(features))==len(features)
                and all(isinstance(f,str) and f.strip() for f in features),
                "Need unique assay source feature IDs (3..10000)")
        _assert(isinstance(spec.get("tissue"),str) and spec["tissue"].strip(),
                "Molecular tissue/compartment cannot be inferred")
        times=spec.get("collected_at",())
        _assert(len(times)==n, "Every omic must be matched to original subject time")
        for i in range(n):
            dt=_date(times[canonical[i]])
            delta=(dose_dates[i]-dt).total_seconds()/86400
            _assert(0<=delta<=max_lookback_days,
                    "Post-landmark or out-of-window omic specimen forbidden")
        X=_numeric(spec.get("values"),rows=n)
        _assert(X.shape[1]==len(features),
                "Molecular matrix cannot be matched by position without source IDs")
        scaled, scale_note=_source_scale(X[canonical],spec.get("value_type",""))
        samples[omic]=(scaled,features)
        source_scales[omic]=scale_note
        if omic in fixed_scores:
            weights=fixed_scores[omic]
            _assert(isinstance(weights,Mapping) and 0<len(weights)<=len(features)
                    and set(weights)<=set(features),
                    "Fixed biological weights must map to original measured feature IDs")
            try:
                w=np.array([float(weights.get(g,0)) for g in features])
            except (TypeError,ValueError,OverflowError):
                raise MultimodalPKBenchmarkError("Invalid predefined biological weights") from None
            _assert(np.isfinite(w).all() and np.sum(np.abs(w))>0,
                    "Weights must be finite, nonzero and externally fixed")
            score_weights[omic]=w/np.sum(np.abs(w))
    order=sorted(samples)
    names=["clinical"]
    if G is not None:
        names+=["clinical_pgx"]
    for omic in order:
        names.extend((f"{omic}_pca2",f"{omic}_direct_ridge"))
        if omic in score_weights:
            names.append(f"{omic}_fixed_score")
    if len(order)>1:
        names+=["early_fusion_pca2","late_fusion_fixed_average"]
    if include_xgboost:
        _assert(n>=80, "XGBoost with molecular factors requires >=80 people")
        try:
            from xgboost import XGBRegressor
        except ImportError:
            raise MultimodalPKBenchmarkError("Install optional benchmark xgboost dependency") from None
        names += ["reference_xgboost","omic_pca_xgboost"]
    predictions={name:np.full(n,np.nan,dtype=float) for name in names}
    # KFold is deterministic on *sorted unique biological IDs*, no assay rows.
    folds=(KFold(n_splits=n_outer_splits,shuffle=True,
                 random_state=20261019).split(np.arange(n))
           if validation_mode=="participant_kfold"
           else LeaveOneGroupOut().split(np.arange(n),groups=groups))
    for train,test in folds:
        pred={"clinical":_fit_ridge(C[train],C[test],y[train])}
        if G is not None:
            pred["clinical_pgx"]=_fit_ridge(B[train],B[test],y[train])
        features_train=[]
        features_test=[]
        single_predictions=[]
        for omic in order:
            x,_=samples[omic]
            a,b=_two_frozen_pcs(x[train],x[test])
            features_train.append(a); features_test.append(b)
            pred[f"{omic}_pca2"]=_fit_ridge(
                np.column_stack((B[train],a)),np.column_stack((B[test],b)),y[train])
            single_predictions.append(pred[f"{omic}_pca2"])
            z=StandardScaler().fit(x[train])
            pred[f"{omic}_direct_ridge"]=_fit_ridge(
                np.column_stack((B[train],z.transform(x[train]))),
                np.column_stack((B[test],z.transform(x[test]))),y[train])
            if omic in score_weights:
                score_train=z.transform(x[train])@score_weights[omic]
                score_test=z.transform(x[test])@score_weights[omic]
                pred[f"{omic}_fixed_score"]=_fit_ridge(
                    np.column_stack((B[train],score_train)),
                    np.column_stack((B[test],score_test)),y[train])
        joint_train=np.column_stack((B[train],*features_train))
        joint_test=np.column_stack((B[test],*features_test))
        if len(order)>1:
            pred["early_fusion_pca2"]=_fit_ridge(joint_train,joint_test,y[train])
            pred["late_fusion_fixed_average"]=np.mean(single_predictions,axis=0)
        if include_xgboost:
            def xgb_predict(a,b):
                z=StandardScaler().fit(a)
                model=XGBRegressor(
                    n_estimators=100,max_depth=2,learning_rate=0.04,
                    subsample=0.9,colsample_bytree=0.9,
                    reg_lambda=10,random_state=20261019,
                    n_jobs=1,verbosity=0)
                model.fit(z.transform(a),np.log(y[train]))
                return np.exp(model.predict(z.transform(b)))
            pred["reference_xgboost"]=xgb_predict(B[train],B[test])
            pred["omic_pca_xgboost"]=xgb_predict(joint_train,joint_test)
        for name in names:
            predictions[name][test]=pred[name]

    _assert(all(np.isfinite(x).all() and np.all(x>0) for x in predictions.values()),
            "Incomplete held-person PK predictions")
    errors={name:np.abs(p-y) for name,p in predictions.items()}
    base_error=errors[baseline]
    rng=np.random.default_rng(20261019)
    resamples=rng.integers(0,n,size=(999,n))
    results={}
    for name in names:
        pred=predictions[name];err=errors[name];diff=err-base_error
        distribution=np.mean(diff[resamples],axis=1)
        results[name]={
            "mae":round(float(err.mean()),7),
            "rmse":round(float(np.sqrt(np.mean((pred-y)**2))),7),
            "bias":round(float(np.mean(pred-y)),7),
            "p20":round(float(np.mean(err/y<=0.2)),6),
            "delta_mae_vs_reference":round(float(diff.mean()),7),
            "delta_mae_bootstrap_95pct_no_refit":[round(float(v),7) for v in
                                                 np.percentile(distribution,[2.5,97.5])],
            "n_heldout_biological_people":n,
        }
    return {
        "status":("SYNTHETIC_METHOD_CONTRACT_NO_HUMAN_PK_EVIDENCE" if source_kind=="synthetic"
                  else "ORIGINAL_SOURCE_DECLARED_EXPLORATORY_NOT_EXTERNAL_VALIDATION"),
        "method":("source_site_disjoint_heldout_fold_fixed_arms"
                  if validation_mode=="leave_one_site_out"
                  else "canonical_person_heldout_kfold_ridge_log_pk_fixed_arms"),
        "validation_mode":validation_mode,
        "n_attested_unique_source_sites":(outer_folds if groups is not None else None),
        "site_identity_caller_attested_not_independently_verified":groups is not None,
        "site_disjoint_holdout_is_not_independent_external_cohort":groups is not None,
        "n_independent_people":n,
        "n_outer_splits":outer_folds,
        "pk_target_name":pk_target_name,
        "pk_target_unit_as_declared":pk_target_unit,
        "source_sha256":source_sha256,
        "genotype_available":G is not None,
        "primary_reference":baseline,
        "modalities":order,
        "n_features_by_modality":{omic:len(samples[omic][1]) for omic in order},
        "original_tissue_and_time_checked_per_modality":True,
        "source_scale_notes":source_scales,
        "fixed_scores_source_author_attested_not_independently_verified":bool(score_weights),
        "early_fusion_only_if_two_or_more_omics":len(order)>1,
        "late_fusion_fixed_no_model_selection":len(order)>1,
        "upstream_batch_preprocessing_provenance_not_verified":True,
        "external_validation":False,
        "physical_units_and_source_originality_declarative_only":True,
        "models":results,
        "no_person_identifiers_values_or_predictions_exported":True,
    }
