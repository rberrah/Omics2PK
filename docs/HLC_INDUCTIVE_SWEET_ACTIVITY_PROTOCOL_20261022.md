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
cohort. Current Synapse metadata reports 465 samples for `syn88644` and 836 for
`syn4629`; these entity-level counts are **not** an automatic complete-case
same-donor count.


### Preferred source join: Synapse-to-Synapse, not GEO-to-Synapse

The primary local acquisition route should use the owner-obtained original
`curatedExpressionLiver` (`syn88644`) and `curatedPhenotype` (`syn4629`)
assets from the same HLC project. If their source donor keys are demonstrably
compatible, this is preferable to forcing a `GSE9588 GSM -> Synapse donor`
crosswalk. GSE9588 remains an independently public expression reference and a
useful numerical stress test, but it is **not required** for the HLC activity
benchmark when the original Synapse expression entity is available.

No row-order join is allowed: the shared biological donor key must be explicit
in the source assets. The historical 466-study / 465-current-expression /
427-GEO / 398-activity counts must remain distinct.

A **previous owner-authorized local audit** of the obtained HLC archives
(reported only as aggregate QC and not rerun in this PR) found 466 expression
identifiers with an exact 466-ID match to the phenotype asset, while
endpoint-specific CYP activities were available for approximately 387–398
donors. This makes the current Synapse metadata value of 465 a provenance
discrepancy to record, **not** a reason to discard an otherwise authenticated
466-row source archive. The source gate therefore accepts only 465 or 466
unique expression donor keys, records which documented count is observed, and
still requires explicit donor-key matching for every modeled activity.

## Why this benchmark is preferable to forcing SWEET onto EPA100

EPA100 is a valuable same-donor RNA-to-activity positive/negative E2 benchmark,
but its RNA layer is targeted qPCR. SWEET is a network method and should be
tested on a genuine broader transcriptome rather than on a tiny already
biologically selected qPCR panel.

HLC provides the appropriate whole-transcriptome setting. It also provides a
harder transportability test because the original paper reports study-site /
processing effects as important covariates.

**Important limitation:** the available curated HLC expression is already a
source-processed representation. Published HLC methods report adjustment for
age, sex and centre of origin before downstream association analyses. Thus a
leave-one-centre-out Omics2PK analysis is a stronger stress test than random
patient folds, but it is **not a pristine raw-array external validation**:
source-level normalization/covariate processing can precede the Omics2PK split.
This limitation must accompany any numerical result.

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
- SWEET v1 training features are subject-disjoint cross-fitted; grouped
  cross-fitting reuses one immutable frozen reference per training fold without
  changing the v1 formula;
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

The owner-provided HLC archives have previously been obtained and a prior
authorized local audit established the aggregate expression/phenotype identity
structure described above. The current automation environment cannot
materialize their raw bytes from the project Library, so **no new HLC
SWEET-vs-RNA numerical result is claimed by this PR**. The next authorized
local run should consume those same original files without committing them and
retain only aggregate metrics and source hashes.
