# Same-donor hepatic RNA + protein + measured CYP activity — source acquisition gate

**Date:** 22 October 2026  
**Scope:** research-only in-vitro hepatic metabolism; not clinical plasma PK.

## Scientific question

On the **same original human liver donors**, does quantitative CYP protein add
held-donor-out information about measured drug-metabolism activity beyond the
corresponding RNA measurement?

The benchmark must compare, on identical donor folds:

1. training mean;
2. RNA only;
3. protein only;
4. RNA + protein.

A positive association in a paper, a scatterplot, group summary, or a
cross-study RNA/protein merge cannot answer this question.

## Priority 1 — Zhang et al. 2016, n=100

**DOI:** 10.1124/jpet.116.233635  
**PMID:** 27189963

The published study reports **100 human normal liver samples** with measurements
covering 10 CYP isoforms, including CYP protein content, mRNA levels and
measured microsomal probe-substrate intrinsic activities. The methods include
midazolam 1′-hydroxylation among the measured CYP activities.

This is currently the strongest acquisition target because it combines:
- a materially larger interindividual cohort than Ohtsuki n=17;
- targeted RNA and protein measurements;
- measured drug-substrate metabolism;
- multiple CYP axes in one source study.

### Important CYP3A interpretation

Midazolam metabolism in human liver microsomes reflects the CYP3A axis and must
not be relabelled as a perfectly isoform-specific CYP3A4 activity measurement
when CYP3A4 and CYP3A5 both contribute. Any eventual analysis must preserve the
source assay wording and report exactly which molecular abundance(s) are used.

### Current source status

**PUBLISHED_MATCHED / ABSTAIN_SOURCE_ROWS_NOT_ACQUIRED.**

The article-level evidence is strong enough to prioritize acquisition, but an
authenticated original donor-level table containing donor IDs + RNA + protein +
activity has **not** been retrieved into the Omics2PK source workflow. No model
result may be claimed from the published correlations or figures.

## Priority 2 — Ohtsuki et al. 2012, n=17

**DOI:** 10.1124/dmd.111.042259  
**PMID:** 21994437

The study reports 17 human liver biopsies with absolute protein quantification
of multiple CYPs and comparison to mRNA expression and functional enzyme
activities. The publisher page explicitly states that supplemental material is
available. The published methods further state that **individual-donor protein
quantification is reported in Supplemental Tables S5 and S6**, and that values
can be under the limit of quantification when fewer than three MRM transitions
are obtained. This is why the Omics2PK source schema now models BLQ explicitly
rather than requiring 17 positive values.

This is a valuable **independent secondary source**, especially if individual
RNA, protein and activity rows can be recovered with stable donor identity. Its
small n makes it less suitable than Zhang n=100 as the primary proof.

Current status: **PUBLISHED_MATCHED / ABSTAIN_SOURCE_ROWS_NOT_ACQUIRED.**

A dedicated public source audit now probes the historical ASPET DOI supplement
landing route
`https://dmd.aspetjournals.org/lookup/suppl/doi:10.1124/dmd.111.042259/-/DC1`
and the current DMD article route. On 22 October 2026 **both returned HTTP 403**
from GitHub Actions. No file, donor row or supplement manifest was recovered.
This result is reproducible with
`scripts/audit_ohtsuki2012_official_supplement.py` and does not imply the
supplement is absent; it establishes that the tested official programmatic
publisher routes are currently inaccessible to this workflow.

## Other candidates

- **Jamwal 2018 n=74:** CYP3A4 mRNA, targeted protein and measured midazolam
  1′-hydroxylation in the same study; authenticated individual matrix still
  unavailable in the current workflow.
- **EPA n≈100:** authentic RNA + measured CYP activity, but no matched protein.
- **Tucatinib/imatinib HLM:** authentic protein + activity/CLint, but no matched
  RNA.
- **Wegler n=15:** protein/proteomics + measured CLint, no qualified matched RNA.
- **GSE9588:** authentic expression matrix and published CYP activities, but no
  qualified RNA↔activity donor crosswalk and no matched protein layer.

None of these may be combined across studies to manufacture a multi-view
cohort.

## Required original-source schema

The local validator
`src/omics2pk/multicyp_rna_protein_activity_source_gate.py` accepts an
authorized UTF-8 CSV with **one row per donor × included CYP**:

`donor_id, source_doi, row_kind, cyp_isoform, probe_substrate, rna_value, rna_unit, rna_status, protein_value, protein_unit, protein_status, activity_value, activity_unit, activity_status`

Allowed measurement statuses are:
- `observed`: numeric value required;
- `blq`: value field must be blank;
- `missing`: value field must be blank;
- `not_measured`: value field must be blank.

Rules:

- `row_kind=individual` only;
- one DOI per file;
- stable biological donor IDs, never row-order matching;
- every source donor remains represented, including explicit BLQ or missing
  measurements; BLQ/missing values are **never coerced to zero**;
- observed measurements must be finite and strictly positive;
- consistent physical/relative units within each CYP assay;
- one prespecified probe activity per CYP benchmark;
- no pooled rows, simulated PBPK rows, figure-digitized pseudo-data or
  reconstructed correlations;
- every included CYP must retain one row for the complete published source
  cohort:
  - Zhang 2016: **100 source donors represented**;
  - Ohtsuki 2012: **17 source donors represented**.

A file may contain a subset of CYP endpoints from the publication, but it may
not contain a post-hoc subset of donors for an included CYP. Numerical fitting
uses only the donors with RNA, protein and activity all observed for that CYP,
on the **same complete-case donor set across every compared arm**. No
imputation is performed. Fewer than 8 complete donors returns
`ABSTAIN_INSUFFICIENT_COMPLETE_DONORS`.

## Preregistered numerical analysis

For each included CYP independently:

- preserve the full source cohort and report RNA/protein/activity missingness;
- define one complete-case donor set requiring all three measurements observed;
- no BLQ substitution, no zero filling and no statistical imputation;
- donor leave-one-out on that fixed complete-case set;
- fixed Ridge alpha=1;
- log1p RNA, protein and measured activity;
- fold-local predictor scaling through the existing verified native Ridge;
- same held donor for RNA, protein and RNA+protein;
- primary contrast:
  **ΔMAE = MAE(RNA+protein) − MAE(RNA)**;
- secondary contrast:
  **MAE(protein) − MAE(RNA)**;
- descriptive paired donor-loss bootstrap;
- aggregate outputs only.

No hyperparameter search, CYP selection by observed gain, or post-hoc choice of
probe substrate is allowed.

## Go / no-go

**GO** only after the original publisher/author file is legitimately obtained,
hashed, and transformed into the exact local schema with source units and donor
IDs preserved.

Until then the correct result is:

`ABSTAIN_SOURCE_ROWS_NOT_ACQUIRED`

The software being ready is not evidence that protein improves RNA, and this
in-vitro benchmark is not evidence of clinical Omics→PK superiority.
