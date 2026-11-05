"""Endpoint-aware heldout prediction diagnostics, with biological unit as resample unit.

Used for NCI-60 original processed doxorubicin GI50 (cell-line PD, not PK).
Never log/antilog unknown source-response units or reinterpret GI50 as AUC/CL.
For actual positive PK endpoints use an independently specified PK assay/CoU
contract; do not treat transformed, censored or unavailable values as positive.
"""
from __future__ import annotations
import numpy as np


class DiagnosticsError(ValueError):
    pass


def _check(y, predictions):
    obs=np.asarray(y,dtype=float)
    pred=np.asarray(predictions,dtype=float)
    if (obs.ndim!=1 or len(obs)<20 or pred.ndim!=2
            or pred.shape[1]!=len(obs) or not 1<=pred.shape[0]<=10
            or not np.isfinite(obs).all() or not np.isfinite(pred).all()
            or np.std(obs)<1e-12):
        raise DiagnosticsError("Distinct measured outcomes and complete heldout repeated predictions required")
    return obs,pred


def _average_tie_ranks(x):
    x=np.asarray(x,dtype=float)
    order=np.argsort(x,kind="stable")
    rank=np.empty(len(x),dtype=float)
    start=0
    while start<len(x):
        stop=start+1
        while stop<len(x) and x[order[start]]==x[order[stop]]:
            stop+=1
        rank[order[start:stop]]=(start+stop+1)/2.0
        start=stop
    return rank


def _corr(a,b):
    if np.std(a)<1e-12 or np.std(b)<1e-12:
        return None
    return float(np.corrcoef(a,b)[0,1])


def _round(x):
    return None if x is None or not np.isfinite(x) else round(float(x),6)


def summarize_heldout(y, predictions, *, endpoint_kind="original_NCI60_processed_doxorubicin_GI50"):
    """Descriptive metrics on original biological units, no line-wise results.

    MAE/RMSE/bias include repeated heldouts; ranks/calibration use the
    *mean of repeated predictions per original line* (not pseudo-replicates).
    This is a pooled OOF diagnostic; no external calibration or CV within
    which calibration has been re-estimated.
    """
    if endpoint_kind!="original_NCI60_processed_doxorubicin_GI50":
        raise DiagnosticsError("No unverified endpoint rescaling or PK metrics")
    y,p=_check(y,predictions)
    e=p-y[None,:]
    abs_e=np.abs(e)
    per_line=abs_e.mean(axis=0)
    pm=p.mean(axis=0)
    pearson=_corr(y,pm)
    spearman=_corr(_average_tie_ranks(y),_average_tie_ranks(pm))
    denom=float(np.sum((pm-pm.mean())**2))
    slope=None if denom<1e-12 else float(
        np.sum((pm-pm.mean())*(y-y.mean()))/denom)
    ccc_den=float(np.var(y)+np.var(pm)+(np.mean(y)-np.mean(pm))**2)
    ccc=None if ccc_den<1e-12 else 2*float(np.cov(y,pm,bias=True)[0,1])/ccc_den
    # All error magnitudes are in the ORIGINAL, processed GI50 source units.
    return {
        "endpoint_kind":endpoint_kind,
        "n_independent_source_cell_lines":len(y),
        "n_repeated_predictions_not_independent":int(p.size),
        "MAE_original_source_scale":_round(np.mean(abs_e)),
        "RMSE_original_source_scale":_round(np.sqrt(np.mean(e*e))),
        "median_per_cell_line_mean_absolute_error":_round(np.median(per_line)),
        "p90_per_cell_line_mean_absolute_error":_round(np.quantile(per_line,.9)),
        "mean_signed_error_prediction_minus_observed":_round(np.mean(e)),
        "mean_absolute_error_by_repeat":[_round(x) for x in np.mean(abs_e,axis=1)],
        "out_of_fold_R2_original_source_scale":_round(
            1-np.mean(e*e)/np.mean((y-y.mean())**2)),
        "Spearman_per_unique_cell_line_repeat_mean_prediction":_round(spearman),
        "Pearson_per_unique_cell_line_repeat_mean_prediction":_round(pearson),
        "concordance_correlation_per_unique_cell_line":_round(ccc),
        "pooled_calibration_slope_observed_on_predicted_descriptive":_round(slope),
        "pooled_calibration_intercept_descriptive":_round(
            y.mean()-slope*pm.mean() if slope is not None else None),
        "no_new_holdout_calibration_performed":True,
        "source_has_not_proven_physical_GI50_unit_or_log_transform":True,
        "symmetric_absolute_errors_may_mask_GI50_response_direction":True,
        "not_independent_external_validation":True,
    }


def paired_heldout_contrast(y, candidate, reference, *, n_draws=1999, seed=20261105):
    """Paired uncertainty from n original cell-line loss units, exploratory only.

    These are percentiles of losses *without fitting* the models again;
    no external replication, model-selection correction or CI proof.
    """
    yy,a=_check(y,candidate)
    _,b=_check(y,reference)
    if a.shape!=b.shape or not isinstance(n_draws,int) or not 99<=n_draws<=9999:
        raise DiagnosticsError("Same independent units and repeated folds required")
    delta=np.abs(a-yy[None,:]).mean(axis=0)-np.abs(b-yy[None,:]).mean(axis=0)
    sampled=np.random.default_rng(seed).integers(0,len(yy),size=(n_draws,len(yy)))
    estimates=delta[sampled].mean(axis=1)
    lo,hi=np.quantile(estimates,[.025,.975])
    return {
        "delta_MAE_candidate_minus_reference_negative_is_better":_round(delta.mean()),
        "descriptive_paired_cell_line_bootstrap_95pct":[_round(lo),_round(hi)],
        "n_independent_line_units_resampled":len(yy),
        "n_bootstrap_draws":n_draws,
        "repeat_errors_averaged_within_cell_line_before_resampling":True,
        "not_refitted_in_bootstrap_or_multiplicity_adjusted":True,
        "not_an_external_significance_test":True,
    }


def multiplicative_pk_diagnostics(y, predictions, *, original_units_verified=False,
                                  positive_endpoint=None, clinically_relevant_ratio=1.2):
    """PK-only metric contract: explicit endpoint, positive ORIGINAL-scale values.

    Never apply to processed GI50 / in-vitro probe activity, negative values,
    logarithms, censoring, or ratios formed from an unverified PK matrix.
    """
    if positive_endpoint not in ("AUC","CL","C0","Cmax"):
        raise DiagnosticsError("Unqualified pharmacokinetic endpoint")
    if original_units_verified is not True:
        raise DiagnosticsError("Unverified original-scale positive PK values")
    yy,p=_check(y,predictions)
    if np.any(yy<=0) or np.any(p<=0):
        raise DiagnosticsError("Multiplicative PK errors require strictly positive original values")
    if not 1<clinically_relevant_ratio<=3:
        raise DiagnosticsError("Prospective fold-accuracy threshold required")
    lr=np.log(p/yy[None,:])
    abs_lr=np.abs(lr)
    ratio=np.exp(lr)
    return {
        "endpoint":positive_endpoint,
        "n_independent_subjects":len(yy),
        "geometric_mean_fold_error_pred_over_observed":_round(np.exp(np.mean(lr))),
        "geometric_mean_absolute_fold_error":_round(np.exp(np.mean(abs_lr))),
        "mean_absolute_log_prediction_error":_round(np.mean(abs_lr)),
        "RMSE_log_ratio":_round(np.sqrt(np.mean(lr**2))),
        "mean_relative_prediction_error_percent":_round(
            100*np.mean(ratio-1)),
        "within_prespecified_fold_ratio_fraction":_round(
            np.mean(abs_lr<=np.log(clinically_relevant_ratio))),
        "prespecified_fold_threshold":clinically_relevant_ratio,
        "metric_derived_from_verified_positive_original_PK":True,
        "patient_level_holdout_repetitions_not_independent":True,
        "bias_and_precision_are_separate_from_accuracy":True,
    }
