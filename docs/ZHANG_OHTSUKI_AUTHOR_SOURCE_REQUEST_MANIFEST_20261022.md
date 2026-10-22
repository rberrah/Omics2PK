# Same-donor CYP RNA + protein + activity — author/source acquisition manifest

**Date:** 22 October 2026  
**Purpose:** obtain the minimum original donor-level source needed to answer the
Omics2PK E2 question without reconstructing data from figures or mixing cohorts.

## Scientific question

On the **same human liver donors**, does quantitative CYP protein improve
held-donor-out prediction of measured CYP probe activity beyond matched CYP RNA?

This is an in-vitro functional-metabolism benchmark. It is not clinical plasma
PK, PopPK or dose prediction.

## Primary request — Zhang et al. 2016

Primary study:

- Zhang et al., *J Pharmacol Exp Ther* 2016
- DOI: `10.1124/jpet.116.233635`
- PMID: `27189963`
- published cohort: **100 human normal liver microsomes**

The primary paper reports mRNA, absolute/quantitative protein content and
intrinsic activities for 10 CYP isoforms in the 100-liver cohort.

A separate open-access paper from the same group, Zhang et al. 2016
(`10.1124/dmd.116.069849`, PMID `27271371`), provides an important
provenance cross-check: it explicitly states that mRNA, protein and activity
levels of the 10 CYPs were **simultaneously quantified in the same set of
100 HLMs**, and describes a database relating POR expression to those CYP
measurements. This strengthens the evidence that the desired donor-level
three-layer matrix existed. It is **not an independent replication cohort** and
its aggregate correlations must not be substituted for the original rows.

### Minimum useful source package

Preferred: the **original publisher supplement / source table exactly as used by
the authors**, including its native donor identifiers and metadata.

If the original table cannot be shared in its native format, a de-identified
extract is sufficient provided one stable pseudonymous donor key is retained
across all three measurement layers.

For every included CYP, request:

1. stable donor ID for **all 100 source donors**;
2. CYP isoform name;
3. probe substrate / activity assay name;
4. original CYP mRNA value and unit or normalization definition;
5. original CYP protein abundance/content and physical unit;
6. original measured microsomal activity or intrinsic-activity value and unit;
7. analytical status for every unavailable value: BLQ, missing or not measured;
8. if available, assay LLOQ/LOQ definition and replicate-summary rule;
9. enough metadata to distinguish microsomal activity per mg protein from
   isoform-normalized activity per pmol CYP.

No direct identifiers, dates of birth, medical record numbers, addresses or
other clinical identity fields are required.

### Preferred scope

The best package contains all ten published CYP axes. If only a subset can be
shared, the donor cohort must still retain **all 100 source donors** for every
included CYP, with unavailable measurements explicitly marked rather than
dropping rows.

For the CYP3A axis, preserve the source wording for the midazolam assay and do
not relabel it as perfectly CYP3A4-specific when CYP3A5 may contribute.

## Secondary request — Ohtsuki et al. 2012

- Ohtsuki et al., *Drug Metab Dispos* 2012
- DOI: `10.1124/dmd.111.042259`
- PMID: `21994437`
- published cohort: **17 human liver biopsies**

The publisher record states that supplemental material exists. The paper
documents individual-donor protein values in Supplemental Tables S5/S6 and
activity-assay information in Supplemental Table S1.

Preferred acquisition is therefore the **complete original supplemental
bundle**, rather than a manually reconstructed extract. If only a reformatted
table is available, use the same schema as below and retain all 17 donors,
including BLQ/missing statuses.

Ohtsuki is an independent secondary cohort and should remain separate from
Zhang. Never pool their donor rows.

## Exact Omics2PK import schema

The public validator
`src/omics2pk/multicyp_rna_protein_activity_source_gate.py` expects UTF-8
long-form CSV with exactly:

`donor_id,source_doi,row_kind,cyp_isoform,probe_substrate,rna_value,rna_unit,rna_status,protein_value,protein_unit,protein_status,activity_value,activity_unit,activity_status`

Allowed statuses:

- `observed` — numeric value required;
- `blq` — value field blank;
- `missing` — value field blank;
- `not_measured` — value field blank.

A header-only template is stored at
`docs/templates/multicyp_same_donor_source_template.csv`.

## Source-integrity rules

Before any model:

- hash the received original file;
- record the DOI and source provenance;
- preserve the complete published donor cohort;
- require one stable donor key across RNA, protein and activity;
- never join by row order;
- never replace BLQ/missing with zero;
- never digitize figures to manufacture donor rows;
- never combine RNA from one cohort with protein/activity from another;
- never export donor-level source rows or individual predictions to the public
  repository.

The source may remain entirely private. Only source hashes, aggregate
missingness and aggregate held-donor-out metrics need to leave the authorized
workspace.

## Preregistered analysis after acquisition

For each eligible CYP, use one identical complete-case donor set for:

1. training mean;
2. RNA only;
3. protein only;
4. RNA + protein.

Fixed donor LODO, Ridge alpha=1, train-only scaling, no imputation.

Primary contrast:

`MAE(RNA + protein) - MAE(RNA)`

Secondary contrast:

`MAE(protein) - MAE(RNA)`

The result remains E2 human-liver functional metabolism and is not evidence of
clinical Omics-to-PK superiority.

## Stop condition

If neither the original supplement nor an author-provided stable-donor extract
can be obtained legitimately, retain:

`ABSTAIN_SOURCE_ROWS_NOT_ACQUIRED`

Do not spend further model-development effort trying to infer the missing
three-layer donor matrix from published summary statistics.
