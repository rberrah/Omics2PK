"""Independent E2 donor-label negative control on ORIGINAL EPA100 CYP data.

Hypothesis: matched source CYP qPCR adds biologically relevant predictive
information beyond chronological age and sex, rather than simply more
degrees of freedom.

Within each CYP, permute original RNA values across ORIGINAL unique liver
donors separately within recorded sex strata. Randomization never sees the
enzyme phenotype and happens before leave-one-DONOR-out fitting. All scaling
is refit on training donors separately for every outer fold. This removes
RNA↔functional activity biological matching while preserving marginal sex-
specific RNA and all age/sex/outcome profiles. Observed and placebo models
use IDENTICAL 100 source subjects and held-donor folds.

Source-specific/finite donor perm tests are not external validation or
clinical PK, cannot prove mechanistic causality. Results are aggregate ONLY.
Original individual donor keys/molecular/activity values never exported.
"""
from __future__ import annotations
from math import expm1, isfinite, log1p
from random import Random
from statistics import mean

from omics2pk.epa_liver_cyp_validation import (
    EPALinkageError, EXPECTED_GENES, _ridge,
)
from omics2pk.epa_eight_cyp_sensitivity import _point_loss_contrast

SEED = 20261018
PERMUTATIONS = 199
N_SOURCE = 100
BONFERRONI_ALPHA = 0.05 / len(EXPECTED_GENES)


def _check_records(records):
    if not (60 <= len(records) <= N_SOURCE):
        raise EPALinkageError("Source original complete biological donor count out of range")
    if len({r[0] for r in records}) != len(records):
        raise EPALinkageError("Duplicate original biological donor key")
    if not all(len(r) == 5 and r[1] > 0 and r[2] >= 0 and r[4] in (0.0,1.0)
               and all(isfinite(v) for v in r[1:])
               for r in records):
        raise EPALinkageError("Invalid original CYP matched RNA/activity/demographics rows")
    for sex in (0.0,1.0):
        if sum(row[4] == sex for row in records) < 3:
            raise EPALinkageError("Insufficient independent donors per sex stratum")


def _permute_rna_within_sex(records, rng):
    """Never shuffle enzyme outcomes, subject keys, sex or age."""
    result = list(records)
    for sex in (0.0,1.0):
        ix=[i for i,r in enumerate(records) if r[4]==sex]
        values=[records[i][1] for i in ix]
        rng.shuffle(values)
        for i, surrogate in zip(ix,values):
            original=result[i]
            result[i]=(original[0], surrogate, *original[2:])
    return result


def _loo_rna_adjusted_mae(rows):
    """Exact same fixed ridge model as source-verified comparator; fold-local."""
    errors=[]
    for h,held in enumerate(rows):
        train=[row for k,row in enumerate(rows) if k!=h]
        x_train=[[r[3],r[4],log1p(r[1])] for r in train]
        train_y=[log1p(r[2]) for r in train]
        x_test=[held[3],held[4],log1p(held[1])]
        pred=_ridge(x_train,train_y,x_test)
        yhat=max(0.,expm1(min(30.,pred)))
        if not isfinite(yhat):
            raise EPALinkageError("Invalid heldout source estimate")
        errors.append(abs(held[2]-yhat))
    return mean(errors)


def evaluate_epa_permutation_negative_control(original, *, n_permutations=PERMUTATIONS,
                                              seed=SEED):
    if set(original)!=set(EXPECTED_GENES):
        raise EPALinkageError("Expected exactly eight prespecified original EPA CYP assays")
    if isinstance(n_permutations,bool) or not 19<=n_permutations<=399:
        raise EPALinkageError("Need 19..399 independently prespecified permuted samples")
    source_ids=set().union(*(set(row[0] for row in original[g]) for g in EXPECTED_GENES))
    if len(source_ids)!=N_SOURCE:
        raise EPALinkageError("Not exactly 100 distinct original human donors")
    summaries={}
    for i,gene in enumerate(EXPECTED_GENES):
        rows=sorted(original[gene],key=lambda r:r[0])
        _check_records(rows)
        observed_losses=_point_loss_contrast(rows)
        ref_mae=mean(p[0] for p in observed_losses.values())
        observed_mae=mean(p[1] for p in observed_losses.values())
        # Verify placebo machinery reproduces actual, unshuffled original
        # predictor within real numerical precision (same fold and solver).
        direct_mae=_loo_rna_adjusted_mae(rows)
        if abs(direct_mae-observed_mae)>1e-7:
            raise EPALinkageError("Source baseline vs negative control machinery disagreement")
        observed_delta=observed_mae-ref_mae
        rng=Random(seed+7919*i)
        permuted_delta=[]
        for _ in range(n_permutations):
            placebo=_permute_rna_within_sex(rows,rng)
            maep=_loo_rna_adjusted_mae(placebo)
            permuted_delta.append(maep-ref_mae)
        # One-sided empirical test: donor-matched RNA improves prediction
        # MORE than as many deliberately mislabeled RNA profiles.
        extreme=sum(value<=observed_delta for value in permuted_delta)
        p=(extreme+1)/(n_permutations+1)
        summaries[gene]={
            "n_unique_original_human_liver_donors":len(rows),
            "baseline_age_sex_mae_original_units":round(ref_mae,6),
            "real_matched_RNA_mae_original_units":round(observed_mae,6),
            "real_matched_RNA_delta_mae":round(observed_delta,6),
            "negative_control_sex_stratified_permuted_RNA_median_delta_mae":round(
                sorted(permuted_delta)[n_permutations//2],6),
            "negative_control_sex_stratified_permuted_RNA_mean_delta_mae":round(
                mean(permuted_delta),6),
            "null_permutations_beating_or_tying_actual_matching":extreme,
            "one_sided_empirical_permutation_p":round(p,6),
            "Bonferroni_8_alpha":BONFERRONI_ALPHA,
            "Bonferroni_8_empirical_permutation_p_below_threshold":p<=BONFERRONI_ALPHA,
            "source_rows_or_individual_predictions_exported":False,
        }
    return {
        "status":"SOURCE_VERIFIED_EPA100_8_CYP_SEX_STRATIFIED_RNA_NEGATIVE_CONTROL_E2_ONLY",
        "n_unique_original_biological_donors":N_SOURCE,
        "n_prespecified_enzyme_endpoints":len(EXPECTED_GENES),
        "number_of_rna_mismatch_permutations_per_endpoint":n_permutations,
        "permutation_seed":seed,
        "fold_design":"identical leave-one-original-donor-out folds, one enzyme at a time",
        "negative_control":"shuffle RNA across donor keys within biological sex before CV; keep target, age, sex and unique source IDs fixed",
        "null":"RNA-to-enzyme pairing contains no donor-specific information beyond what survives within-sex shuffling",
        "one_sided_empirical_p_formula":"(1 + # placebo delta <= observed delta)/(1 + B)",
        "n_bonferroni_adjusted_positive_RNA_matching_associations":sum(
            row["real_matched_RNA_delta_mae"]<0 and row["Bonferroni_8_empirical_permutation_p_below_threshold"]
            for row in summaries.values()),
        "endpoints":summaries,
        "not_patient_PK_or_PopPK":True,
        "not_external_source_replication":True,
        "not_biological_causality":True,
        "original_individual_donor_rows_ever_exported":False,
        "limitations":[
            "Eight assay endpoints share up to 100 original biological donors, not 800 independent patients",
            "Within-sex shuffle preserves marginal sex-specific RNA but not the RNA-age relationship",
            "The permutation uses RNA measurements from all study donors without consulting phenotypes; train fold scaler and regression never use heldout phenotype",
            "Observed and permuted models use fixed ridge(alpha=1), no training-tuned feature selection or SWEET",
            "In-vitro microsomal measured probe-enzyme activity, NOT clinical drug concentration, AUC, systemic PK or PopPK",
            "The matched evidence comes from one public human liver donor panel; not an independent cohort",
        ],
    }
