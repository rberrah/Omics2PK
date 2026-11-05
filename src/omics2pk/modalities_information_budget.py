"""Minimal multi-omics modality budget based on inner-heldout gain.

For a panel with RNA/protein/metabolites, choose the cheapest combination
preserving 80/90% of best positive inner-CV improvement vs no-omic baseline.
No model is selected using outer heldout biological outcomes.
"""
from __future__ import annotations
import math

class ModalityBudgetError(ValueError):pass

def choose_minimum_modalities(inner_loss_by_arm, omic_modes_by_arm,
                             *,fraction,reference="tumor_only",
                             source_assay_cost_weights=None):
    if not 0<fraction<=1 or isinstance(fraction,bool):
        raise ModalityBudgetError("Invalid gain retention fraction")
    if not isinstance(inner_loss_by_arm,dict) or not isinstance(omic_modes_by_arm,dict):
        raise ModalityBudgetError("Need original model losses and source assays")
    if set(inner_loss_by_arm)!=set(omic_modes_by_arm) or reference not in inner_loss_by_arm:
        raise ModalityBudgetError("Reference and each source modality arm mandatory")
    if omic_modes_by_arm[reference]!=():
        raise ModalityBudgetError("Reference must have no molecular assay")
    values={}
    for name,loss in inner_loss_by_arm.items():
        modes=omic_modes_by_arm[name]
        if not isinstance(modes,tuple) or len(set(modes))!=len(modes) or not all(
                x in ("rna","prot","met") for x in modes):
            raise ModalityBudgetError("Unknown or duplicated original source assay")
        if isinstance(loss,bool) or not isinstance(loss,(int,float)) or not math.isfinite(loss) or loss<0:
            raise ModalityBudgetError("Model loss not source-qualified")
        values[name]=float(loss)
    weights=source_assay_cost_weights or {"rna":1.0,"prot":1.0,"met":1.0}
    if set(weights)!={"rna","prot","met"} or any(
            not isinstance(x,(int,float)) or isinstance(x,bool)
            or not math.isfinite(x) or x<=0 for x in weights.values()):
        raise ModalityBudgetError("Positive per-omic cost weights must cover exactly the measured modalities")
    base=values[reference]
    best=max(0.0,max(base-loss for loss in values.values()))
    if best<=0:
        return {"arm":reference,"n_omic_modalities":0,"status":"NO_POSITIVE_INNER_CV_GAIN",
                "best_inner_gain":0.0,"target_fraction":float(fraction),
                "selected_arm_inner_gain":0.0}
    eligible=[arm for arm,loss in values.items() if arm!=reference
              and base-loss>=fraction*best-1e-12]
    arm=min(eligible,key=lambda a:(sum(weights[m] for m in omic_modes_by_arm[a]),
                                    len(omic_modes_by_arm[a]),values[a],a))
    return {"arm":arm,"n_omic_modalities":len(omic_modes_by_arm[arm]),
            "status":"INNER_CV_MINIMUM_OMICS_COST_CANDIDATE",
            "best_inner_gain":best,"target_fraction":float(fraction),
            "selected_arm_inner_gain":base-values[arm],
            "selected_weighted_assay_cost":sum(weights[m] for m in omic_modes_by_arm[arm]),
            "source_assay_cost_weights":weights,
            "selected_only_using_inner_training_data":True}
