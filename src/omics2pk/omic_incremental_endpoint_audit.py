"""Evaluate *actual* held-out added-omic prediction vs genetic baseline.

This module reports matched out-of-fold prediction *descriptives* and a
predeclared metric. It never infers a genetic effect from observed omics,
never converts enzyme activity into in-vivo human PK, never replaces a
patient genotype panel with source-supplied one-SNP pseudo-diplotypes,
and never claims clinical impact from CV alone.

PK positive endpoint units are original/source-provided: MAE/RMSE share
those units. Mean absolute log ratio is scale-robust but presupposes
strictly positive observations and predictions. P20 is prediction
accuracy, NOT clinical target attainment. 95% intervals are subject-
level bootstrap of heldout errors, without refitting/selection adjustment.
"""
from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping,Sequence
from math import isfinite,log,sqrt
from random import Random
from statistics import mean

class OmicIncrementalEvidenceError(ValueError):
    pass

TARGETS={
    "human_plasma_AUC":"CLINICAL_PK",
    "human_plasma_C0":"CLINICAL_PK",
    "human_plasma_CL":"CLINICAL_PK",
    "human_plasma_Cmax":"CLINICAL_PK",
    "human_renal_clearance":"CLINICAL_PK",
    "human_hepatic_enzyme_activity":"BIOLOGICAL_FUNCTION_NOT_CLINICAL_PK",
    "PHH_imatinib_invitro_clearance":"IN_VITRO_DRUG_METABOLISM_NOT_CLINICAL_PK",
}
ERRORS=("mean_absolute_log_error","MAE","RMSE","P20")
GENO={"source_complete_diplotype","verified_partial_PGx","one_measured_SNP","not_assessed"}
REGIMENS={"true_out_of_fold_per_biological_subject","independent_external_cohort"}
STRATA={"known_no_function_variant","known_reduced_function_variant",
        "known_function_variant","unresolved_missense",
        "no_reported_variant_in_assayed_sites","untyped"}


def _valid_float(value,*,key):
    if isinstance(value,bool):
        raise OmicIncrementalEvidenceError(key+" cannot be bool")
    try:out=float(value)
    except (ValueError,TypeError,OverflowError) as exc:
        raise OmicIncrementalEvidenceError(key+" must be a real positive number") from exc
    if not isfinite(out) or out<=0:
        raise OmicIncrementalEvidenceError(key+" must be positive and finite")
    return out


def _metrics(y,p):
    if len(y)!=len(p) or len(y)<5:
        raise OmicIncrementalEvidenceError("Unmatched source PK outcome/prediction lengths")
    diff=[v-t for v,t in zip(p,y)]
    ratio=[v/t for v,t in zip(p,y)]
    ll=[log(r) for r in ratio]
    return {
        "MAE":mean(abs(v) for v in diff),
        "RMSE":sqrt(mean(v*v for v in diff)),
        "mean_absolute_log_error":mean(abs(v) for v in ll),
        "root_mean_squared_log_error":sqrt(mean(v*v for v in ll)),
        "signed_log_bias":mean(ll),
        "P20":mean(.8<=r<=1.2 for r in ratio),
        "within_factor_two":mean(.5<=r<=2.0 for r in ratio),
    }


def _losses(y,p,primary):
    if primary=="MAE":
        return [abs(a-b) for a,b in zip(p,y)]
    if primary=="RMSE":
        # Paired sqrt(MSE) is nonlinear; compare *MSE*, then express
        # the observed difference as RMSE and bootstrap paired RMSE.
        return [(a-b)**2 for a,b in zip(p,y)]
    if primary=="mean_absolute_log_error":
        return [abs(log(a/b)) for a,b in zip(p,y)]
    if primary=="P20":
        return [int(not (.8<=a/b<=1.2)) for a,b in zip(p,y)]
    raise OmicIncrementalEvidenceError("Unqualified primary PK metric")


def _contrast(y,ref,arm,primary,boots=499):
    def stat(indices,values):
        return (sqrt(mean(values[i] for i in indices))
                if primary=="RMSE" else mean(values[i] for i in indices))
    x=_losses(y,arm,primary)
    b=_losses(y,ref,primary)
    ids=tuple(range(len(y)))
    observed=stat(ids,x)-stat(ids,b)
    draws=[]
    random=Random(20261105)
    for _ in range(boots):
        chosen=[random.randrange(len(y)) for _ in y]
        draws.append(stat(chosen,x)-stat(chosen,b))
    draws.sort()
    lo=draws[int(.025*(boots-1))]
    hi=draws[int(.975*(boots-1))]
    return {
        "delta_primary_loss_negative_favors_added_omics":round(observed,6),
        "descriptive_subject_bootstrap_95pct_interval":[round(lo,6),round(hi,6)],
        "interval_includes_zero":lo<=0<=hi,
        "n_bootstrap":boots,
        "bootstrap_refitted_any_model":False,
        "bootstrap_accounts_for_omic_selection":False,
    }


def audit_omic_incremental_PK_prediction(
    *,records,endpoint,original_units,genotype_qualifier,
    evaluation_regime,predeclared_primary_metric,
    source_omic_sample_prelandmark_verified,
    source_prediction_training_is_subject_disjoint,
    original_source_omics_and_PK_same_subject_verified,
    qualified_original_specimen_site_evidence=None,
):
    """Research-only, source-attested *matched* original held-out predictions.

    records: one per biological subject with keys
      subject_id, fold_id, observed, genotype_stratum,
      predictions {clinical_plus_genotype: number,
                   <each prespecified omic arm>: number}.
    Each source subject appears once, all model arms same heldout fold,
    and all predicted numbers must be >0 (strictly positive PK).
    """
    if endpoint not in TARGETS:
        raise OmicIncrementalEvidenceError("Endpoint not qualified as positive PK/function outcome")
    if not isinstance(original_units,str) or not original_units.strip():
        raise OmicIncrementalEvidenceError("Actual original endpoint units required")
    if genotype_qualifier not in GENO or evaluation_regime not in REGIMENS:
        raise OmicIncrementalEvidenceError("Original genotype and evaluation scope must be explicit")
    if genotype_qualifier=="not_assessed":
        raise OmicIncrementalEvidenceError(
            "Cannot compare additional omics against genetic baseline if no original genotype was measured")
    if predeclared_primary_metric not in ERRORS:
        raise OmicIncrementalEvidenceError("Outcome-appropriate primary metric must be prespecified")
    flags=(
        source_omic_sample_prelandmark_verified,
        source_prediction_training_is_subject_disjoint,
        original_source_omics_and_PK_same_subject_verified,
    )
    if any(flag is not True for flag in flags):
        raise OmicIncrementalEvidenceError(
            "No original predose/fit/PK-subject join verification: ABSTAIN")
    if not isinstance(records,(list,tuple)) or not 10<=len(records)<=10000:
        raise OmicIncrementalEvidenceError("Need >=10 independent source biological subjects")
    names=None;seen=set();folds=set();rows=[]
    for raw in records:
        if not isinstance(raw,Mapping):
            raise OmicIncrementalEvidenceError("Each original subject must have an independent source record")
        sid=raw.get("subject_id")
        fold=raw.get("fold_id")
        if not isinstance(sid,str) or not sid.strip() or sid in seen:
            raise OmicIncrementalEvidenceError("Duplicate or missing original biological subject")
        if not isinstance(fold,(str,int)) or isinstance(fold,bool) or str(fold)=="":
            raise OmicIncrementalEvidenceError("Original heldout fold identity required")
        seen.add(sid);folds.add(str(fold))
        y=_valid_float(raw.get("observed"),key="measured original positive outcome")
        variant_stratum=raw.get("genotype_stratum","untyped")
        if variant_stratum not in STRATA:
            raise OmicIncrementalEvidenceError("No auto-inference of known mutation functional consequences")
        px=raw.get("predictions")
        if not isinstance(px,Mapping) or "clinical_plus_genotype" not in px:
            raise OmicIncrementalEvidenceError("Source clinical+actual genotype reference must be included")
        if names is None:
            names=tuple(sorted(px))
            if len(names)<2 or len(names)>20:
                raise OmicIncrementalEvidenceError("One predeclared omic predictor arm minimum")
            if any(not isinstance(a,str) or not a for a in names):
                raise OmicIncrementalEvidenceError("Invalid arm name")
        if tuple(sorted(px))!=names:
            raise OmicIncrementalEvidenceError("Model arms must predict exactly the same subject outcome")
        preds={arm:_valid_float(px[arm],key="original heldout "+arm) for arm in names}
        rows.append((sid,str(fold),variant_stratum,y,preds))
    if len(folds)<3 and evaluation_regime=="true_out_of_fold_per_biological_subject":
        raise OmicIncrementalEvidenceError("At least three original subject-heldout partitions required")
    # No nested-CV verification from a table of predictions alone: this
    # evaluates PRESPECIFIED model arms, not optimized choice of best arm.
    y=[r[3] for r in rows]
    allp={arm:[r[4][arm] for r in rows] for arm in names}
    compared={}
    for arm in names:
        vals=_metrics(y,allp[arm])
        compared[arm]={k:round(v,6) for k,v in vals.items()}
    # Numeric original held-out predictions alone cannot prove where
    # each omic was physically sampled. The independent source/assay
    # catalogue must qualify the *same exact model arm names*, not a
    # guessed organ based on gene symbol or metabolite name.
    specimen_verified=False
    if qualified_original_specimen_site_evidence is not None:
        e=qualified_original_specimen_site_evidence
        if (not isinstance(e,Mapping) or
                e.get("status")!="QUALIFIED_ASSAY_TISSUE_SITE_AND_ORIGINAL_STUDY_FOR_PRESPECIFIED_ARMS"
                or not isinstance(e.get("arm_original_collection_sites"),Mapping)
                or set(e["arm_original_collection_sites"])!=set(names)):
            raise OmicIncrementalEvidenceError(
                "Missing or inconsistent actual collected specimen site for a compared original omic arm")
        if e.get("endpoint_class")!=TARGETS[endpoint]:
            raise OmicIncrementalEvidenceError(
                "Sampled ex-vivo enzyme and actual clinical patient PK endpoint cannot be exchanged")
        if not e.get("source_site_does_not_prove_actual_original_subject_time_join"):
            raise OmicIncrementalEvidenceError(
                "Specimen catalogue cannot assert it independently verified biological donor identity")
        specimen_verified=True
    baseline=allp["clinical_plus_genotype"]
    contrasts={
        arm:_contrast(y,baseline,allp[arm],predeclared_primary_metric)
        for arm in names if arm!="clinical_plus_genotype"
    }
    # Never claim subgroup-specific modality superiority when few
    # original independent carriers exist. Original assayed negative
    # panel is not a confirmed wildtype functional enzyme.
    strata=defaultdict(list)
    for i,r in enumerate(rows):strata[r[2]].append(i)
    per_mutation={}
    for group,ids in sorted(strata.items()):
        if len(ids)<20:
            per_mutation[group]={
                "n_observed_source_subjects":len(ids),
                "status":"ABSTAIN_TOO_FEW_INDEPENDENT_VARIANT_STRATUM_SUBJECTS",
                "no_omic_winner_inferred":True,
            }
            continue
        sy=[y[i] for i in ids]
        sp={arm:[allp[arm][i] for i in ids] for arm in names}
        per_mutation[group]={
            "n_observed_source_subjects":len(ids),
            "status":"DESCRIPTIVE_WITHIN_KNOWN_VARIANT_STRATUM_NOT_CAUSAL",
            "source_genotype_phasing_and_function_not_verified_by_this_module":True,
            "compared_prespecified_arm_errors":{
                arm:{key:round(_metrics(sy,sp[arm])[key],6)
                     for key in ("MAE","mean_absolute_log_error","P20")}
                for arm in names
            },
            "no_omic_winner_inferred":True,
        }
    return {
        "status":("EXTERNAL_SOURCE_HELDOUT_DESCRIPTIVES_NOT_DOSE_BENEFIT"
                  if evaluation_regime=="independent_external_cohort"
                  else "PRESPECIFIED_OUT_OF_FOLD_DESCRIPTIVES_NOT_EXTERNAL_VALIDATION"),
        "n_original_unique_biological_subjects":len(rows),
        "n_original_heldout_groups":len(folds),
        "endpoint":endpoint,
        "endpoint_class":TARGETS[endpoint],
        "original_outcome_units":original_units,
        "original_PGx_coverage":genotype_qualifier,
        "geno_reference_was_complete_diplotype":genotype_qualifier=="source_complete_diplotype",
        "predeclared_primary_loss":predeclared_primary_metric,
        "same_source_subjects_same_prediction_fold_all_arms":True,
        "baseline_clinical_plus_genotype":True,
        "original_sampling_site_and_assay_method_qualified":specimen_verified,
        "original_sampling_site_missing_requires_ABSTAIN_mechanistic_claim":not specimen_verified,
        "qualified_sampling_site_is_not_real_patient_PK_validation":True,
        "measured_prelandmark_molecular_values_required":True,
        "heldout_descriptives_by_predeclared_model_arm":compared,
        "added_omic_minus_genotype_primary_loss_paired_bootstrap":contrasts,
        "within_variant_function_groups":per_mutation,
        "selected_best_omic_claimed":False,
        "cannot_claim_incremental_benefit_without_external_replication_or_selection_adjustment":True,
        "nested_hyperparameter_or_modality_selection_independently_verified":False,
        "clinical_dose_recommendation_or_improved_patient_PK_proven":False,
        "original_subject_ids_outcomes_predictions_or_genetic_rows_emitted":False,
    }
