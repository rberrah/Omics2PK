"""Same-donor hepatic CYP proteomics vs functional enzyme-phenotype fusion.

Tucatinib and imatinib are metabolized in human liver microsomes, with
original-source published donor-matched CYP2C8 and CYP3A abundance, enzyme
probe activities and observed apparent in-vitro clearance.

This is NOT a transcriptome assay and 'proteins+activity' is NOT multi-omics.
No original donor data are bundled. Never pool the same 21 liver donors as
42 biologically independent study participants across two drugs.
"""
from __future__ import annotations

from math import isfinite
from random import Random
from statistics import mean

from omics2pk.epa_liver_cyp_validation import _ridge as _any_dim_ridge
from omics2pk.matched_tki_metabolism import MatchedMetabolismError, _lodo, qualify_pair

SEED=20261018
SOURCE_LABELS=(
    (1, "tucatinib_app_HLM_clearance_uL_min_mg"),
    (2, "imatinib_app_HLM_clearance_quantifiable_only_uL_min_mg"),
    (3, "imatinib_N_desmethyl_formation_peak_area_ratio_not_CL"),
)
ARMS=("cohort_mean","CYP_proteins","CYP_probe_activities","CYP_proteins_plus_activities")


def _prediction(records, held, arm, idx):
    if arm=="cohort_mean":
        return mean(r[idx] for r in records)
    if arm=="CYP_proteins":
        vec=lambda r:r[4]
    elif arm=="CYP_probe_activities":
        vec=lambda r:r[5]
    elif arm=="CYP_proteins_plus_activities":
        vec=lambda r:(*r[4],*r[5])
    else:
        raise MatchedMetabolismError("Unknown original-source modality")
    ytrain=[r[idx] for r in records]
    try:
        pred=_any_dim_ridge([vec(r) for r in records],ytrain,vec(held))
    except (TypeError,ArithmeticError,ValueError) as e:
        raise MatchedMetabolismError("Invalid train-only ridge design") from e
    if not isfinite(pred):
        raise MatchedMetabolismError("Nonfinite source-trained prediction")
    return max(0.0,pred)


def _losses(records, idx):
    if len(records)<12 or len({r[0] for r in records})!=len(records):
        raise MatchedMetabolismError("Too few unique independent original donors")
    scores={k:[] for k in ARMS}
    for i,held in enumerate(records):
        training=[r for j,r in enumerate(records) if j!=i]
        for name in ARMS:
            p=_prediction(training,held,name,idx)
            scores[name].append(abs(p-held[idx]))
    return scores


def _quantile(values,p):
    values=sorted(values)
    z=(len(values)-1)*p
    low=int(z)
    up=min(low+1,len(values)-1)
    return values[low]+(values[up]-values[low])*(z-low)


def run_fusion_on_verified_same_donors(paired):
    """Caller already source-qualified pairs; no original donor data output."""
    if len(paired)!=21 or len({p[0] for p in paired})!=21:
        raise MatchedMetabolismError("Expected exactly 21 original HLM biological donors")
    if sum(p[2] is None for p in paired)!=4:
        raise MatchedMetabolismError("Four original imatinib nonquantifiable clearances must remain missing")
    results={}
    for idx,label in SOURCE_LABELS:
        eligible=[r for r in paired if r[idx] is not None]
        errs=_losses(eligible,idx)
        # Verify established original-source 2-variable models are numerically
        # identical to newly constructed common 4-variable code.
        baseline=_lodo(paired,idx)["heldout_donor_MAE"]
        expected={
            "cohort_mean":"baseline",
            "CYP_proteins":"protein",
            "CYP_probe_activities":"activity"
        }
        for arm,old in expected.items():
            if abs(mean(errs[arm])-baseline[old])>0.000002:
                raise MatchedMetabolismError("Established original single-modality donor-heldout parity failed")
        delta=[
            x-y for x,y in zip(errs["CYP_proteins_plus_activities"],errs["CYP_proteins"])
        ]
        # Descriptive paired error resampling across *real original biological
        # donor units*. No refitting; not an external uncertainty interval.
        rng=Random(SEED+idx)
        bootstrap=sorted(mean(delta[rng.randrange(len(delta))] for _ in delta)
                         for _ in range(1999))
        results[label]={
            "n_original_donors_same_for_all_arms":len(eligible),
            "observed_in_vitro_metabolic_measurement":label,
            "source_endpoint_missingness_imputed":False,
            "MAE_by_fixed_modality_arm_original_units":{
                k:round(mean(v),6) for k,v in errs.items()
            },
            "delta_MAE_fusion_minus_proteins":round(mean(delta),6),
            "exploratory_95pct_paired_loss_bootstrap_delta":[
                round(_quantile(bootstrap,.025),6),
                round(_quantile(bootstrap,.975),6),
            ],
            "n_donors_fusion_lower_absolute_error_than_protein":sum(x<0 for x in delta),
            "n_donors_fusion_higher_absolute_error_than_protein":sum(x>0 for x in delta),
            "original_donor_level_losses_and_predictions_not_emitted":True,
        }
    return {
        "status":"ORIGINAL_21_SHARED_DONOR_HEPATIC_PROTEIN_ACTIVITY_FUSION_E2_IN_VITRO",
        "original_biological_donors_total":21,
        "n_actual_source_modalities":2,
        "precise_modalities":["targeted CYP2C8/CYP3A proteomic protein abundance",
                               "functional CYP2C8/CYP3A probe enzyme-activity phenotyping"],
        "transcriptomic_RNA_data_for_these_21_donors_present":False,
        "metabolomics_data_for_these_21_donors_present":False,
        "can_test_RNA_only_vs_RNA_plus_proteins_on_same_donors":False,
        "not_a_multiomic_comparison":True,
        "frozen_models":"Ridge(alpha=1), 2 or 4 features, training-fold centering/scaling, donor-level LODO",
        "independent_external_derived_drug_clinical_PK":False,
        "same_donors_across_compounds_not_independent_external_replication":True,
        "n_original_donor_rows_exported":0,
        "outcomes":results,
    }


def from_original_source_rows(tucatinib_rows,imatinib_rows):
    paired,checks=qualify_pair(tucatinib_rows,imatinib_rows,verify_reference=True)
    return {
        **run_fusion_on_verified_same_donors(paired),
        "original_imatinib_source_rank_correlations_qualified":checks,
        "original_source_DOIs":[
            "10.1016/j.dmd.2025.100061",
            "10.1016/j.dmd.2025.100196",
        ]
    }
