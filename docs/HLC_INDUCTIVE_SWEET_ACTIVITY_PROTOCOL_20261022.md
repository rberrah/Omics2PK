# HLC — inductive SWEET v1 to measured CYP activity benchmark

## Purpose

The Human Liver Cohort (HLC; Synapse project `syn4499`) is the strongest
currently identified public-source candidate for a **whole-transcriptome
functional benchmark** of the experimental Omics2PK SWEET v1 extension.

This is an **E2 functional metabolism benchmark**, not clinical PK.

Primary sources:
- Yang et al. 2010, *Genome Research*, DOI `10.1101/gr.103341.109`;
- GEO `GSE9588`;
- Synapse HLC project `syn4499`;
- current project child entities include `syn88644`
  (`curatedExpressionLiver`), `syn89614` (`curatedGenotype`) and
  `syn4629` (`curatedPhenotype`).

The source article reports microsomal activities for nine CYP enzymes in
398 HLC liver samples, with two CYP3A4 probe substrates (midazolam and
testosterone), and source-normalized liver expression from the same parent
cohort. The cohort sizes at the project/entity level are **not** an automatic
same-donor crosswalk.

## Why this benchmark is preferable to forcing SWEET onto EPA100

EPA100 is a valuable same-donor RNA-to-activity positive/negative E2 benchmark,
but its RNA layer is targeted qPCR. SWEET is a network method and should be
tested on a genuine broader transcriptome rather than on a tiny already
biologically selected qPCR panel.

HLC provides the appropriate whole-transcriptome setting. It also provides a
harder transportability test because the original paper reports study-site /
processing effects as important covariates.

## Source gate before any numerical claim

A real run is permitted only after an authenticated local join establishes:

1. one unique biological donor per expression column;
2. an original HLC donor key linking expression and measured CYP activity;
3. an explicit original study/tissue-collection site for every donor;
4. age and sex from the original source;
5. one source-annotated direct CYP transcript/probe for the endpoint;
6. a physical activity unit and named probe substrate;
7. no row-order matching, figure digitisation or cross-study fusion.

The public repository must never contain the original donor crosswalk,
individual expression values or individual activity measurements.

## Prespecified primary validation

**Outer validation: leave one original HLC site out.**

Site is used only to define the external-style outer fold and is **not** passed
as a predictor. This deliberately asks whether a model fitted on the remaining
sites transfers to a biologically independent source collection / processing
context.

Within every outer training set:

- source-normalized expression is passed through unchanged;
- optional edge genes are selected by deterministic variance using only outer
  training RNA;
- SWEET v1 training features are subject-disjoint cross-fitted;
- the held-out site is transformed/networked only from the frozen outer
  training reference;
- model standardization is fitted on outer training donors only;
- no CYP activity value enters preprocessing, edge selection, graph
  construction or SWEET sample-weight calibration.

## Prespecified model arms

All arms use identical eligible donors and the same held-site folds:

1. `age_sex`;
2. `age_sex_plus_direct_cyp_rna`;
3. `age_sex_plus_train_selected_rna`;
4. `age_sex_plus_train_selected_pca4`;
5. `age_sex_plus_inductive_sweet_v1`.

The direct RNA arm is essential: SWEET must beat or complement the actual CYP
transcript, not merely an age/sex baseline.

The selected-RNA and PCA arms are equally important because a network should
not be credited for signal that a simpler leakage-controlled expression
representation already captures.

Primary SWEET contrasts:

- `MAE(SWEET) - MAE(direct CYP RNA)`;
- `MAE(SWEET) - MAE(selected RNA Ridge)`;
- `MAE(SWEET) - MAE(PCA4)`.

Lower is better. Results must include held-site SWEET clamp and similarity-shift
QC.

## Interpretation

A positive result would establish only that an inductive network representation
contains held-site information about **measured microsomal enzyme function**.
It would not establish:

- in-vivo systemic clearance;
- tacrolimus AUC or dose requirement;
- superiority over authentic CYP3A genotype in patients;
- clinical benefit;
- external validation in an independent study.

A null or worse SWEET result is scientifically useful. It would support using
direct RNA / low-dimensional molecular scores rather than a network for the
downstream Omics2PK clinical benchmark.

## Current execution status — 22 October 2026

The software benchmark is implemented in
`src/omics2pk/hlc_sweet_activity_benchmark.py` with synthetic leakage tests.

The owner-provided HLC archives have previously been obtained, but the current
automation environment cannot materialize their raw bytes from the project
Library. Therefore **no new HLC numerical SWEET-vs-RNA result is claimed by
this PR**. The next authorized local run should consume those original files
without committing them and retain only aggregate metrics and source hashes.
