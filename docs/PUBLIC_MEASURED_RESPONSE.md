# Public measured-response validation — NCI-60 multidrug

The **mixOmics `multidrug`** dataset combines **60 NCI-60 tumor cell
lines**, expression of **48 ABC transporter genes measured by RT-qPCR**,
and experimentally measured activity of **1,429 compounds** (GI50).

Sources: https://github.com/mixOmics-org/mixOmics
and Szakács et al. (2004), Cancer Cell.

**Endpoint is in-vitro drug sensitivity, NOT clinical PK, clearance or AUC.**
This provides an initial real-observation pharmacological benchmark, not
evidence that transcriptomics predicts individual PopPK.

Planned comparisons with subject-disjoint nested folds: intercept-only,
RNA expression + Ridge, PCA on RNA + Ridge, exploratory SWEET-inspired
network + Ridge, and RNA + network + Ridge. Any optional nonlinear models
must have identical data access and tuning folds. Include shuffled-subject
network negative controls and report uncertainty across splits.

The first gate only audits upstream RData and lists matching compounds.
No compound is selected based on model performance. A named compound must
be selected/frozen after the audit and **before** predictions are run.
Never select the best-scoring compound out of 1,429 post hoc.
