# Wegler 2021 — source structure verified, donor-CLint still unverified

Date: **2026-10-16**. Public GitHub Actions: https://github.com/rberrah/Omics2PK/actions/runs/38922776553

## Source provenance
Wegler C et al., *Molecular Pharmaceutics* 2021, DOI [10.1021/acs.molpharmaceut.1c00053](https://doi.org/10.1021/acs.molpharmaceut.1c00053), PMCID PMC8041379. The paper studies donor-matched HH and HLM in-vitro intrinsic clearance (CLint) and proteomics in 15 human liver donors for five CYP probe substrates. This is **NOT** a clinical PopPK dataset of dose and patient plasma concentrations.

All three supporting files were downloaded temporarily through official PMC AWS (version PMC8041379.1), checked against fixed SHA-256 hashes and inspected **without storing donor IDs, individual protein values or CLint outcomes** in GitHub CI artifacts.

## Observed source structure

| Public support | Structural result | Scientific interpretation |
|---|---|---|
| Supplemental PDF | 21 pages, all with native text; 4 distinct table-number references, 10 distinct figure-number references | Contents searchable, but **no donor-level CLint table confirmed** |
| Data S1 — Hepatocytes | 4,785 XML rows; 34 cell elements in each of first four rows | Large global proteomics matrix; rows/cells are **not independent donors** |
| Data S1 — Microsomes | 5,173 XML rows; 34 cell elements in each of first four rows | Separate compartment-specific global proteome |
| Data S1 — Description | 8 XML rows | Metadata, not a measured PK outcome table |
| Data S2 — bioinformatics | Original contrast worksheets and fixed hash verified | Cannot be presumed to contain donor-level PK outcomes |

A PDF keyword match, clearance plot, workbook sheet, or correlation does not prove that original numerical CLint observations per donor and compound can be matched to source proteome columns. The source audit explicitly leaves donor-level CLint table verification and numerical outcome extraction at FALSE.

## Blocking scientific work, issue #11

1. Independently review the original PDF: which tables include exact **donor × five drugs × HH/HLM** numerical CLint values (not merely figure axes or group summaries)? Record original assay units, unbound/total basis, missingness/censoring and technical repeats.
2. Identify provenance-supported donor IDs in the S1 Hepatocytes and Microsomes matrices and prove the correspondence to CLint measurements. 34 populated cells in an XML row do **not** mean 34 donors.
3. Only then build a donor/compound/system crosswalk; if exact per-donor source CLint is not tabulated, request original outcome tables from study authors. Image/plot digitization, if used, must be a separately labelled **low-fidelity exploratory** sensitivity, never a "measured original data" benchmark.
4. Once eligible, run prespecified donor-disjoint leave-one-donor-out drug-specific reference vs measured-CYP protein vs carefully defined HLM→HH transfer models with clustered uncertainty and all null findings.

## Thesis fit
- Chapter 4: molecular functional phenotype and IVIVE.
- Chapter 10: molecular protein–in-vitro CLint proof; **not** yet RNA→clinical PK.
- Chapter 13: provenance, mediation and abstention when source data cannot be transported.
- Clinical R2 and future dosing require real linked individual dose/time/concentration records, not this liver-protein bench.
- Negative NCI-60 GI50 results remain PD-related, not PK.
