# HLC / GSE9588 site-disjoint SWEET functional benchmark protocol

Updated: **22 October 2026**

## Scientific question

On the same human liver specimens, does an inductive **train-frozen SWEET-inspired**
RNA representation improve prediction of measured CYP metabolic activity beyond:

1. a training-outcome mean;
2. the prespecified corresponding CYP transcript;
3. a low-dimensional broad-RNA PCA comparator?

This is an **E2 functional-metabolism benchmark**. It is not patient plasma PK,
PopPK, MIPD or a clinical validation claim.

## Why HLC is now the preferred SWEET functional target

The Human Liver Cohort underlying GSE9588 contains global liver expression and
published microsomal CYP activity measurements. The literature also reports
substantial study-site effects in this cohort. Therefore a random donor split is
not the primary validation design: it could reward site/batch recognition.

The primary design is **leave-one-source-site-out**.

## Frozen protocol

For each held-out study site:

- all donors from that site are excluded before any fitted molecular operation;
- RNA imputation / scaling are fit on training sites only;
- PCA basis is fit on training sites only;
- SWEET reference and sample-weight calibration are fit on training sites only;
- SWEET features for training donors are themselves leave-one-training-donor-out
  cross-fitted before the activity regression;
- the held-out site is predicted with the frozen training artifact;
- fixed Ridge alpha=1 is used for every learned outcome arm;
- activity is log1p transformed for fitting and evaluated on its original scale;
- no endpoint-driven gene selection, network tuning or hyperparameter search is
  permitted.

The SWEET arm is the Omics2PK **inductive query-anchored v1**, not the published
transductive SWEET algorithm.

## Fixed arms

| Arm | Representation |
|---|---|
| training_logmean | no molecular predictor |
| target_RNA | prespecified matching CYP transcript |
| PCA | broad RNA, PCA fit on outer-training sites |
| SWEET_v1 | train-frozen inductive SWEET gene-strength vector |

Primary comparison:

`MAE(SWEET_v1) - MAE(target_RNA)`

Secondary:

`MAE(PCA) - MAE(target_RNA)`

A negative delta favors the richer representation.

## Source gate

A numerical run requires an authenticated same-donor extract containing:

- unique biological donor identity;
- source study site / collection center;
- global liver expression;
- measured CYP activity and physical assay units;
- explicit mapping from expression sample to the same biological liver;
- a prespecified CYP target and a prespecified edge panel or a separately
  qualified training-only selector.

No row-order join, figure digitization, reconstructed correlations or
cross-study fusion is allowed.

The user's previously obtained HLC Synapse archives are appropriate candidate
sources, but they must be materialized and structurally audited before fitting.
The public repository never stores donor-level source rows.

## Current status

- train-frozen SWEET software: **qualified at the software-contract level**;
- real GSE9588 expression stress test: **qualified**;
- site-disjoint HLC benchmark implementation: **implemented with synthetic CI**;
- authenticated HLC expression + activity + site rows in the current execution:
  **not materialized**;
- real HLC SWEET-vs-target-RNA result: **ABSTAIN pending source materialization**.

A positive result would support a transcriptomic network representation as a
functional CYP phenotype proxy. It would still not establish improved clinical
tacrolimus PK prediction; that remains an E3/E4 question requiring paired
patient omics and drug exposure.
