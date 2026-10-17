"""Original 100-donor EPA eight-CYP research sensitivity: shared-donor multiplicity.

This is reanalysis of 8 dependent *in-vitro* probe-enzyme endpoints, not
patient PK. The original two-ID linked donor observations are NEVER exported.
Eight point-estimate comparisons and paired resampling are performed across
the SAME 100 source people; endpoint missingness stays explicit.

All model settings fixed before inspecting bootstrap outcomes:
  8 original CYPs, alpha=1, log1p original activity and source qPCR,
  leave-one-BIOLOGICAL-donor-out, train-only feature scaling,
  5000 shared-donor resamples (seed 20261017);
  Bonferroni percentile 99.375% per endpoint (0.05/8).
"""
from __future__ import annotations

from math import expm1, isfinite, log1p
from random import Random
from statistics import mean

from omics2pk.epa_liver_cyp_validation import (
    EPALinkageError, EXPECTED_GENES, _ridge, validate_one_gene,
)

SEED=20261017
RESAMPLES=5000
ALPHA=0.05
PER_CYP_ALPHA=ALPHA/len(EXPECTED_GENES)


def _point_loss_contrast(records):
    """No subject values returned to report: use ONLY transient donor-keyed map."""
    if not 60<=len(records)<=100:
        raise EPALinkageError("Need >=60 original biological donors")
    keys=[r[0] for r in records]
    if len(set(keys))!=len(keys):
        raise EPALinkageError("Duplicate original biological identifier")
    losses={}
    for held, row in enumerate(records):
        train=[r for i,r in enumerate(records) if i!=held]
        ylog=[log1p(r[2]) for r in train]
        clinical=[row[3],row[4]]
        mrna=log1p(row[1])
        p_ref=_ridge([[r[3],r[4]] for r in train],ylog,clinical)
        p_full=_ridge([[r[3],r[4],log1p(r[1])] for r in train],
                      ylog,[*clinical,mrna])
        reff=max(0.,expm1(min(30,p_ref)))
        full=max(0.,expm1(min(30,p_full)))
        if not all(isfinite(x) for x in (reff,full)):
            raise EPALinkageError("Nonfinite out-of-fold predictions")
        losses[row[0]]=(abs(row[2]-reff),abs(row[2]-full))
    return losses


def _quantile(values,q):
    vals=sorted(values)
    if not vals:
        raise EPALinkageError("No bootstrap distribution")
    # Fixed linear percentile interpolation, deterministic and no scipy.
    idx=(len(vals)-1)*q
    low=int(idx)
    high=min(low+1,len(vals)-1)
    return vals[low]+(idx-low)*(vals[high]-vals[low])


def calculate_eight_gene_sensitivity(original_records, *, n_resamples=RESAMPLES,
                                     seed=SEED):
    if set(original_records)!=set(EXPECTED_GENES):
        raise EPALinkageError("Eight prespecified CYP endpoints exactly required")
    if n_resamples<100 or isinstance(n_resamples,bool):
        raise EPALinkageError("Need >=100 prespecified source resamples")
    by_gene={}
    for gene in EXPECTED_GENES:
        records=original_records[gene]
        loss=_point_loss_contrast(records)
        original=validate_one_gene(records)
        baseline=round(mean(r[0] for r in loss.values()),6)
        expanded=round(mean(r[1] for r in loss.values()),6)
        expected=original["mae_original_units"]
        if (abs(baseline-expected["age_sex"])>0.000001
                or abs(expanded-expected["age_sex_plus_rna"])>0.000001):
            raise EPALinkageError(
                "Source original EPA eight-gene LOOCV parity failed; abstain"
            )
        by_gene[gene]=loss
    source_keys=set().union(*(set(v) for v in by_gene.values()))
    if len(source_keys)!=100:
        raise EPALinkageError("Not exactly 100 total original liver donors")
    for gene, losses in by_gene.items():
        if not set(losses).issubset(source_keys):
            raise EPALinkageError("Invalid source donor key")

    # Shared cluster index across all eight CYPs, not eight independent cohorts.
    ordered_keys=sorted(source_keys)
    rng=Random(seed)
    samples={gene:[] for gene in EXPECTED_GENES}
    for i in range(n_resamples):
        draw=[ordered_keys[rng.randrange(len(ordered_keys))]
              for _ in ordered_keys]
        for gene, source in by_gene.items():
            delta=[source[key][1]-source[key][0] for key in draw if key in source]
            if len(delta)<60:
                raise EPALinkageError("Too few complete source donors in resample")
            samples[gene].append(mean(delta))

    results={}
    count_nominal=count_bonferroni=0
    for gene, source in by_gene.items():
        delta=[row[1]-row[0] for row in source.values()]
        estimate=mean(delta)
        lower95=_quantile(samples[gene],0.025)
        upper95=_quantile(samples[gene],0.975)
        lb=_quantile(samples[gene],PER_CYP_ALPHA/2)
        ub=_quantile(samples[gene],1-PER_CYP_ALPHA/2)
        sig95=(lower95>0 or upper95<0)
        sig_multi=(lb>0 or ub<0)
        count_nominal+=int(sig95 and estimate<0)
        count_bonferroni+=int(sig_multi and estimate<0)
        results[gene]={
            "n_original_donors":len(source),
            "mae_delta_plus_rna_minus_age_sex":round(estimate,6),
            "paired_percentile_CI_95_unadjusted":[round(lower95,6),
                                                 round(upper95,6)],
            "paired_percentile_CI_99_375_Bonferroni_8":[round(lb,6),round(ub,6)],
            "nominal_CI_excludes_zero":sig95,
            "Bonferroni_CI_excludes_zero":sig_multi,
            "RNA_MAE_point_gain":estimate<0,
            "n_donors_RNA_lower_absolute_error":sum(x<0 for x in delta),
            "worst_single_loss_difference":round(max(delta),6),
            "best_single_loss_difference":round(min(delta),6),
        }
    return {
        "status":"ORIGINAL_EPA_100_LIVER_CYP_8_MULTIPLICITY_SENSITIVITY",
        "original_biological_donors":100,
        "n_related_CYP_endpoints":len(EXPECTED_GENES),
        "n_prespecified_shared_donor_bootstrap_draws":n_resamples,
        "random_seed":seed,
        "familywise_method":"Eight Bonferroni-adjusted per-phenotype 99.375% percentile intervals",
        "bootstrap_unit":"resample the SAME 100 biological donor IDs for all eight CYPs",
        "n_negative_delta_MAE_point_estimates":sum(
            d["RNA_MAE_point_gain"] for d in results.values()),
        "n_unadjusted_excluding_zero_improvements":count_nominal,
        "n_eight_endpoints_adjusted_excluding_zero_improvements":count_bonferroni,
        "results":results,
        "limitations":[
            "Percentile paired bootstrap of previously fitted LODO losses, not full refitted nested bootstrap",
            "Bonferroni controls nominal familywise error only under interval coverage assumptions",
            "Clustered eight endpoints are correlated within identical donors, not eight independent biological cohorts",
            "Small E2 in-vitro qPCR->enzyme activity; not clinical PK, omics->PopPK, genotype-adjusted or causal proof",
            "Selected qPCR genes had preexisting biologically specified enzyme matches; no global feature search",
        ],
        "original_individual_donor_records_exported":False,
    }
