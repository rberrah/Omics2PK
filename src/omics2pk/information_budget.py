"""Omics2PK: choose minimal omic feature budget from *validated gain*, not SHAP mass.

Selection is allowed only in INNER cross-validation. A trained model's SHAP
absolute-value share is a descriptive explanation, NOT a portion of its gain.
No positive gain over clinical reference => choose 0 omic features.

This file has no dataset, clinical records, or model fitting dependency.
"""
from __future__ import annotations
import math
from collections import Counter
from itertools import combinations

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


def summarize_train_panel_stability(panels, n_features, *, max_reported=20):
    """Descriptive feature-panel overlap across *training-derived* model fits.

    All candidate selection and test performance evaluation MUST be done
    separately, within nested CV. In LOOCV, training sets overlap heavily:
    high Jaccard is not independent replication or biological confirmation.
    Report empty panels separately so they cannot spuriously inflate overlap.
    """
    if isinstance(n_features, bool) or not isinstance(n_features, int) or n_features <= 0:
        raise InformationBudgetError("Need positive original measured feature count")
    if isinstance(max_reported, bool) or not isinstance(max_reported, int) or max_reported < 1:
        raise InformationBudgetError("Need positive reporting limit")
    if not isinstance(panels, (list, tuple)) or len(panels) < 2:
        raise InformationBudgetError("At least two train-derived panels required")
    normalized = []
    for panel in panels:
        if not isinstance(panel, (list, tuple)):
            raise InformationBudgetError("Panel must be an ordered source-index collection")
        indices = list(panel)
        if (any(isinstance(i, bool) or not isinstance(i, int) or i < 0 or i >= n_features
                for i in indices) or len(set(indices)) != len(indices)):
            raise InformationBudgetError("Duplicate or non-source molecular feature")
        normalized.append(set(indices))
    nonempty = [panel for panel in normalized if panel]
    intersections, jaccards, random_intersections = [], [], []
    for left, right in combinations(nonempty, 2):
        common = len(left & right)
        intersections.append(common)
        jaccards.append(common / len(left | right))
        # Uniform random sets with these same panel sizes: descriptive
        # expectation only. Correlated omics / overlapping donor folds are
        # NOT covered by this null, so do not calculate a p-value.
        random_intersections.append(len(left) * len(right) / n_features)
    frequencies = Counter(i for panel in nonempty for i in panel)
    top = sorted(frequencies.items(), key=lambda item: (-item[1], item[0]))[:max_reported]
    return {
        "status": "TRAIN_ONLY_PANEL_STABILITY_DESCRIPTIVE_NOT_EXTERNAL_VALIDATION",
        "n_model_fits": len(normalized),
        "n_nonempty_panels": len(nonempty),
        "n_zero_omic_panels": len(normalized) - len(nonempty),
        "n_pairs_both_nonempty": len(jaccards),
        "mean_pairwise_Jaccard_conditional_on_both_nonempty":
            round(sum(jaccards) / len(jaccards), 6) if jaccards else None,
        "mean_pairwise_intersection_conditional_on_both_nonempty":
            round(sum(intersections) / len(intersections), 6) if intersections else None,
        "uniform_random_sets_expected_intersection_matched_panel_sizes":
            round(sum(random_intersections) / len(random_intersections), 6)
            if random_intersections else None,
        "source_feature_selection_frequency": [
            {"feature_index": idx, "n_training_panels": count,
             "fraction_all_training_panels": round(count / len(normalized), 6),
             "fraction_nonempty_training_panels": round(count / len(nonempty), 6)}
            for idx, count in top
        ],
        "overlapping_train_folds_invalidate_independence_or_significance_claims": True,
        "used_to_choose_panel_or_evaluate_predictive_gain": False,
    }
