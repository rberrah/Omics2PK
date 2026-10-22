# Zhang 2016 / Ohtsuki 2012 same-donor RNA + protein + activity acquisition gate

Updated: **22 October 2026**

## Question

Can Omics2PK compare, on the **same original human liver donors**, whether
quantitative CYP protein improves prediction of measured drug-metabolism
activity beyond matched CYP RNA?

This document records what is genuinely source-supported and what is still
missing. It is an acquisition gate, not a fitted result.

## A. Zhang et al. 2016 — primary n=100 target

**DOI:** 10.1124/jpet.116.233635  
**PMID:** 27189963

The original article reports **100 human normal liver samples** with large
interindividual variation in:
- CYP mRNA;
- CYP protein content;
- intrinsic CYP activities;
- selected CYP genotypes.

The current publisher/ScienceDirect record also states that supplemental
material exists.

### What is NOT yet qualified

No current public route inspected on 22 October 2026 yielded an authenticated,
downloadable donor-level supplement containing stable donor IDs plus matched
RNA + protein + activity values.

A dedicated bounded source probe subsequently tested both the historical ASPET
JPET supplement route
`https://jpet.aspetjournals.org/lookup/suppl/doi:10.1124/jpet.116.233635/-/DC1`
and the current publisher article route. From GitHub Actions, no legitimate
supplement file or manifest was acquired. This is an **access result**, not
evidence that the source data do not exist.

Therefore:

`ABSTAIN_SOURCE_ROWS_NOT_ACQUIRED`

Do not reconstruct the matrix from figures, correlation coefficients,
aggregate ranges, row order, related publications or later reuse of the same
liver bank.

## B. Ohtsuki et al. 2012 — secondary n=17 target

**DOI:** 10.1124/dmd.111.042259  
**PMID:** 21994437

The publisher full text provides unusually strong evidence that the desired
biology exists in the original experiment:

- **17 human liver biopsies**;
- absolute CYP/UGT protein quantification;
- matched mRNA expression;
- functional enzyme activities;
- protein values below LOQ are possible;
- **individual-donor protein data are explicitly stated to be presented in
  Supplemental Tables S5 and S6**;
- enzyme-activity assay details are stated to be in **Supplemental Table S1**.

The methods also state that RNA assays were run as technical triplicates.

### Why this is still not enough

The live publisher text proves study design and individual-donor protein
availability, but the complete supplemental asset itself has not been
materialized in the current workflow. A bounded source probe tested the exact
historical ASPET supplement endpoint
`https://dmd.aspetjournals.org/lookup/suppl/doi:10.1124/dmd.111.042259/-/DC1`
and the current DMD article route; **both returned HTTP 403** from GitHub
Actions. No supplement bytes, donor rows or file manifest were obtained. This
is an automated-access blocker, not evidence of biological-data absence.

The article text alone does not prove that one downloadable table contains:
- stable donor keys;
- individual RNA values;
- individual protein values;
- individual activity values;
- units/status fields sufficient for the preregistered LODO comparison.

Therefore Ohtsuki remains:

`ABSTAIN_SOURCE_ROWS_NOT_ACQUIRED`

It is nevertheless a **high-priority acquisition target** because individual
protein rows are explicitly documented and the cohort is small enough for a
manual source-integrity audit once the original supplement is legitimately
obtained.

## C. Required source contract before any fit

The existing Omics2PK gate requires long-form donor × CYP rows with:

`donor_id, source_doi, row_kind, cyp_isoform, probe_substrate, rna_value, rna_unit, rna_status, protein_value, protein_unit, protein_status, activity_value, activity_unit, activity_status`

Allowed measurement statuses:
- `observed`
- `blq`
- `missing`
- `not_measured`

Hard rules:
1. every published source donor remains represented;
2. the same donor-ID cohort is retained across included CYP endpoints;
3. BLQ/missing values are explicit and never coerced to zero;
4. no figure digitisation or row-order matching;
5. no cross-study fusion;
6. numerical comparison uses one fixed complete-case donor set per CYP;
7. fewer than 8 complete donors returns
   `ABSTAIN_INSUFFICIENT_COMPLETE_DONORS`.

## D. Preregistered comparison

For each eligible CYP, on identical donor folds:

1. training mean;
2. RNA only;
3. protein only;
4. RNA + protein.

Primary contrast:

`ΔMAE = MAE(RNA+protein) − MAE(RNA)`

Secondary contrast:

`MAE(protein) − MAE(RNA)`

Fixed donor-LODO, Ridge alpha=1, train-only scaling, aggregate outputs only.

## E. Current conclusion

The software and statistical gate are ready. Public literature now confirms
that both Zhang and Ohtsuki contain the relevant biological measurements, and
Ohtsuki explicitly documents individual-donor protein supplements.

What is still missing is **the authenticated donor-level numerical source
matrix itself**. The next legitimate route is authenticated institutional /
publisher access or an author-provided original file, followed by hashing and
the existing strict donor-level source gate.

Until it is obtained legitimately, RNA vs protein vs RNA+protein superiority
remains **unestablished**.
