# CYP3A5 functional activity: split-conformal uncertainty and target leakage (2026-11-03)

## Scope

Research-only extension of `benchmark_cyp3a5_activity`. The response is an independently **measured, validated CYP3A5-specific activity** endpoint in explicit source units. This is not a Bayesian model, a mechanistic clearance estimate, or a tacrolimus dosing algorithm. Source provenance and assay validation are still *declared*, not independently verified by this software.

## Same-assay leakage gate

Each input modality must declare `assay_relation_to_outcome` as either `independent_assay` or `independent_specimen`. `same_activity_assay_derived` and `unknown` are refused. In particular, a metabolite directly generated or quantified as part of the activity assay **cannot be recycled as an independent predictor of that activity**. Real source qualification must authenticate the declaration with original protocols, biospecimen IDs, processing history and specimen dates. Even if a metabolite is measured in an independent assay, overlapping biological pathways do not establish CYP3A5 specificity.

Conventional pre-dose PK prediction has an even stricter time criterion: no postdose pharmacometabolomic value or downstream concentration-derived feature may be used to predict a future PK outcome. A contemporaneous E2 activity-phenotype association is not automatically a prospective E3 PK predictor.

## Optional uncertainty audit

Set `interval_alpha=0.10` to request nominal 90% intervals. For each of the same donor-heldout outer folds and each comparator arm:

1. Independently split only the outer-training donors into a **proper-fit** set and a disjoint **calibration** set; use the *same* randomized split for all arms. All input transformations (Ridge mean/SD and regression coefficients) are fit exclusively on the proper-fit donors.
2. Fit fixed Ridge(alpha=1) on log measured activity; calculate out-of-proper-fit absolute log-residuals on calibration donors.
3. Use conformal order statistic rank `ceil((n_cal + 1) × (1 − alpha))`; **abstain** rather than pretend to have a finite 90%/99% interval if insufficient calibration donors remain while keeping at least eight proper training subjects.
4. For the held-out donors, predict with the **same proper-fit model**, convert the log-scale radius back to an interval on the original activity scale, and report only *aggregate* coverage, mean interval width and MAE. The original point-MAE comparators remain fitted on full outer-training folds and must not be misidentified as the centers of these intervals.

Split-conformal prediction has a marginal finite-sample coverage rationale **only if donors are exchangeable and the proper-fit/calibration/test separation is respected**. Repeated outer-fold audits do **not** establish externally validated, conditional/genotype-specific, center-robust or longitudinal uncertainty. These are predictive intervals for measurements, not biological posterior credible intervals on latent CYP3A5 function.

## Synthetic software audit only

A deterministic synthetic 60-donor dataset with prespecified signal was used solely to check code execution, input order invariance, abstention, and leakage prevention. For requested 90% marginal intervals the observed held-out coverage values were:

| Arm | Synthetic empirical coverage | Synthetic mean interval width (arbitrary simulated activity units) |
|---|---:|---:|
| Genotype only | 0.867 | 2.060 |
| Genotype + transcriptomics | 0.917 | 0.500 |
| Genotype + proteomics | 0.900 | 0.581 |
| Genotype + metabolomics | 0.900 | 0.481 |
| Genotype + all omics | 0.933 | 0.443 |
| Genotype + knowledge-selected omics | 0.900 | 0.406 |

These are deliberately simulated values, **not** human enzyme or PK validation. The metabolomic synthetic signal was built to be more closely associated with the simulated activity, so these numbers are not empirical evidence for metabolomic superiority. For nominal 99% intervals in a small synthetic n=30 cohort, the module correctly reports `ABSTAIN_INSUFFICIENT_PROPER_TRAIN_OR_CALIBRATION_DONORS` rather than fabricated finite intervals.

## Necessary next evidence

Before a genuine omics+knowledge E2 study, obtain original same-donor genotype (including *6/*7 coverage), isoform- or peptide-resolved measurements, independently assayed CYP3A5-specific activity and a documented metabolite/enzymatic reaction map. CYP3A shared probe activities and CYP3A5-specific outcomes must not be treated as identical. For the E3 stage, use original *pre-landmark* omics, medication histories, time-resolved tacrolimus PK and a true clinical+PGx comparator on identical held-out patients. Both stages require source-validated IDs and missingness, original LOQ treatment and examination of genotype subgroup performance.
