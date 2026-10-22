# Ohtsuki 2012 — official supplemental-material acquisition audit

**Date:** 22 October 2026  
**Scope:** source acquisition for the same-donor RNA / protein / CYP-activity
comparison in issue #63. This is not a fitted model.

## Why this source matters

Ohtsuki et al., *Drug Metabolism and Disposition* 2012, DOI
https://doi.org/10.1124/dmd.111.042259, studied **17 human liver biopsies** and
compared hepatic mRNA expression, absolute protein abundance and functional
enzyme activities.

The current publisher article text establishes several unusually useful source
facts:

- individual-donor targeted protein measurements are reported in
  **Supplemental Tables S5 and S6**;
- protein measurements with fewer than three usable MRM transitions are defined
  as **under the limit of quantification**, rather than zero;
- microsomal P450 assay conditions are in **Supplemental Table S1**;
- nine P450 enzyme activities were measured with specific probe substrates;
- the paper explicitly directs readers to the supplemental data for the
  **absolute activity data**;
- CYP3A4 activity was assayed with testosterone, midazolam and nifedipine.

These facts justify Ohtsuki as a genuine publication-level same-liver
multi-view candidate and justify the explicit BLQ/missingness policy in
`multicyp_rna_protein_activity_source_gate.py`.

## Exact official access test

The following publisher routes were preregistered and probed from GitHub
Actions with bounded HTTP requests:

1. historical ASPET DOI supplemental landing route  
   `https://dmd.aspetjournals.org/lookup/suppl/doi:10.1124/dmd.111.042259/-/DC1`
2. current DMD article route  
   `https://dmd.aspetjournals.org/article/S0090-9556%2824%2922586-1/abstract`

Result on 22 October 2026:

- historical supplemental route: **HTTP 403**;
- current article route: **HTTP 403**;
- discovered official supplemental file links: **0**;
- original supplemental bytes acquired: **NO**;
- donor rows parsed: **NO**;
- RNA/protein/activity crosswalk qualified: **NO**;
- model fitted: **NO**.

Reproducible source-only workflow:
GitHub Actions run **39792724502**.

Audit implementation:
`scripts/audit_ohtsuki2012_official_supplement.py`.

The workflow stores only metadata such as response state, URL shape and file
type if available. It never uploads publisher supplement bytes or donor-level
measurements.

## Scientific interpretation

This is an **access blocker, not a biological-data absence finding**.

The publication establishes that individual protein data and absolute enzyme
activity data exist in the supplement, but the tested official publisher routes
do not currently permit automated retrieval. Therefore the correct Omics2PK
state remains:

`PUBLISHED_MATCHED / ABSTAIN_SOURCE_ROWS_NOT_ACQUIRED`

The supplement must not be reconstructed from scatterplots, reported
correlations or group statistics.

## Next legitimate routes

1. authenticated institutional/publisher access to the original supplement;
2. an author- or institution-hosted original supplemental file with provenance;
3. direct author request for the original donor-level tables.

After acquisition, the source must be hashed and checked for stable donor IDs,
RNA units, protein units/status, probe-specific activity units/status and the
same donor cohort before the preregistered LODO benchmark can run.

Ohtsuki remains the strongest **independent replication** candidate. Zhang 2016
n=100 remains the preferred primary acquisition target because of its larger
sample size.
