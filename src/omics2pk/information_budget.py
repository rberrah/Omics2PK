"""Omics2PK: choose minimal omic feature budget from *validated gain*, not SHAP mass.

Selection is allowed only in INNER cross-validation. A trained model's SHAP
absolute-value share is a descriptive explanation, NOT a portion of its gain.
No positive gain over clinical reference => choose 0 omic features.

This file has no dataset, clinical records, or model fitting dependency.
"""
from __future__ import annotations
import math

DEFAULT_BUDGETS = (0, 1, 2, 3, 5, 8, 12, 20, 40, 80, 160, 384)
TARGETS = (0.8, 0.9)


class InformationBudgetError(ValueError):
    pass


def _finite_number(value):
    if isinstance(value, bool):
        raise InformationBudgetError("Non-numeric loss")
    try:
        f = float(value)
    except (ValueError, TypeError, OverflowError):
        raise InformationBudgetError("Non-numeric loss") from None
    if not math.isfinite(f) or f < 0:
        raise InformationBudgetError("Loss must be finite and nonnegative")
    return f


def validate_budgets(values, n_features):
    if isinstance(n_features, bool) or not isinstance(n_features, int) or n_features <= 0:
        raise InformationBudgetError("Source number of measured features required")
    if not isinstance(values, (tuple, list)) or not values:
        raise InformationBudgetError("Budget schedule required")
    if any(isinstance(k, bool) or not isinstance(k, int) for k in values):
        raise InformationBudgetError("Budgets must be integer feature counts")
    result=tuple(sorted(set(values)))
    if result[0] != 0 or result[-1] != n_features or any(k < 0 or k > n_features for k in result):
        raise InformationBudgetError("Budget schedule must include 0 and full original feature count")
    return result


def select_gain_retention_budget(loss_by_k, fraction, *, require_positive_gain=True):
    """Pick smallest k preserving fraction of best positive inner-validated gain.

    Input loss_by_k is produced from ONE outer TRAINING set by independent
    inner-fold heldout predictions. The outer subject must not enter ranking,
    loss estimates, choice of k, or estimate of best gain.
    """
    if isinstance(fraction, bool) or not isinstance(fraction, (int,float)) or not 0 < fraction <= 1:
        raise InformationBudgetError("Target must be a fraction from 0 to 1")
    if not isinstance(loss_by_k, dict) or 0 not in loss_by_k:
        raise InformationBudgetError("True clinical-only k=0 reference is mandatory")
    values={}
    for k,loss in loss_by_k.items():
        if isinstance(k,bool) or not isinstance(k,int) or k<0:
            raise InformationBudgetError("Invalid omic budget")
        values[k]=_finite_number(loss)
    ref=values[0]
    improvement={k:ref-v for k,v in values.items()}
    best=max(0.0,max(improvement.values()))
    if best<=0:
        return {"k":0,"status":"NO_POSITIVE_INNER_VALIDATED_GAIN",
                "reference_inner_loss":ref,"best_positive_inner_gain":0.0,
                "gain_retention_fraction":float(fraction),
                "selected_inner_gain":0.0,
                "min_features_selected":True}
    threshold=best*float(fraction)
    candidate=[k for k in sorted(values) if k>0 and improvement[k]>=threshold-1e-12]
    if not candidate:
        raise InformationBudgetError("No candidate preserves validated positive gain")
    k=candidate[0]
    return {"k":k,"status":"INNER_VALIDATED_GAIN_RETENTION_CANDIDATE",
            "reference_inner_loss":ref,"best_positive_inner_gain":best,
            "gain_retention_fraction":float(fraction),
            "selected_inner_gain":improvement[k],
            "min_features_selected":True}


def cumulative_shap_coverage(absolute_shap_importance, selected_indices):
    """Descriptive sum(|TreeSHAP|) in a *training-fold* full omic model."""
    vals=[_finite_number(x) for x in absolute_shap_importance]
    if not vals:
        raise InformationBudgetError("Missing measured feature importances")
    selected=list(selected_indices)
    if len(selected)!=len(set(selected)) or any(
            isinstance(j,bool) or not isinstance(j,int) or j<0 or j>=len(vals)
            for j in selected):
        raise InformationBudgetError("Invalid or repeated feature selection")
    total=sum(vals)
    if total==0:
        return {"status":"NO_MOLECULAR_SHAP_ATTRIBUTION",
                "fraction_of_absolute_SHAP":None,"fraction_of_predictive_gain":None}
    return {"status":"DESCRIPTIVE_MODEL_ATTRIBUTION_NOT_PREDICTIVE_GAIN",
            "fraction_of_absolute_SHAP":sum(vals[j] for j in selected)/total,
            "fraction_of_predictive_gain":None}


def pareto_frontier(model_points):
    """Find nondominated (n_omic_features, heldout_loss) points.

    Model 0 may dominate all omic models: a scientifically valid outcome.
    Uses descriptive heldout loss; NEVER for hyperparameter/threshold selection.
    """
    if not isinstance(model_points, dict) or not model_points:
        raise InformationBudgetError("Feature budget model points missing")
    points=[]
    for label,item in model_points.items():
        if not isinstance(label,str) or not label or not isinstance(item,dict):
            raise InformationBudgetError("Each candidate needs a named model and its loss")
        k=item.get("k")
        if isinstance(k,bool) or not isinstance(k,int) or k<0:
            raise InformationBudgetError("Invalid feature count")
        loss=_finite_number(item.get("loss"))
        points.append({"model":label,"k":k,"loss":loss})
    points.sort(key=lambda p:(p["k"],p["loss"],p["model"]))
    frontier=[]
    for p in points:
        if not any(q["k"]<=p["k"] and q["loss"]<=p["loss"] and
                   (q["k"]<p["k"] or q["loss"]<p["loss"])
                   for q in points):
            frontier.append(p)
    return {"criterion":"minimize (number of measured omic features, heldout prediction loss)",
            "non_dominated_models":frontier,
            "do_not_choose_using_outer_test_losses":True}
