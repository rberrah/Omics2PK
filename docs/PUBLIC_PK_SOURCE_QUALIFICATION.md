# Public PK/proteomics qualification: Wegler 2021 and Yin/Battail 2025

**Status (16 October 2026): source-format inventory, not a donor-level PK fit or clinical validation.**

## 1. Priority: measured in-vitro hepatic clearance and proteomics (Wegler 2021)

Wegler C, Matsson P, Krogstad V, et al. "Influence of Proteome Profiles and Intracellular Drug Exposure on Differences in CYP Activity in Donor-Matched Human Liver Microsomes and Hepatocytes", Molecular Pharmaceutics 2021, DOI [10.1021/acs.molpharmaceut.1c00053](https://doi.org/10.1021/acs.molpharmaceut.1c00053), PMID 33739838, PMCID PMC8041379.

The study describes 15 original human liver donors, liver microsomes (HLM), hepatocytes (HH) and proteomics, with **measured intrinsic hepatic clearance (CLint)** of five probe compounds.

| Probe compound | CYP marker (not necessarily exclusive enzyme) |
|---|---|
| Bupropion | CYP2B6 |
| Diclofenac | CYP2C9 |
| Omeprazole | CYP2C19 |
| Bufuralol | CYP2D6 |
| Midazolam | CYP3A4, with possible other CYP3A contribution |

Official PMC supporting files:
- mp1c00053_si_001.pdf — supporting experimental methods, kinetics and donor tables;
- mp1c00053_si_002.xlsx — Data S1, global proteomics in hepatocytes/microsomes/tissue/pellet;
- mp1c00053_si_003.xlsx — Data S2, proteomics bioinformatics;
- PRIDE PXD014131 — raw LC-MS/MS, which is not needed for a first model.

**Important 2026 retrieval change:** In August 2026, NLM/PMC retired its legacy article-dataset FTP and cloud package distribution. The supported public channel is now the PMC AWS Cloud Service (bucket pmc-oa-opendata), whose article metadata includes per-object supplementary URLs and license data. The audit therefore probes the official PMC8041379.1 metadata object as well as historical article bin URLs; individual older downloads can legitimately fail even though the archived article and new cloud objects remain available. Source metadata access by itself proves no quantitative donor-level values. See https://pmc.ncbi.nlm.nih.gov/tools/pmcaws/ and https://pmc-oa-opendata.s3.amazonaws.com/README.txt .

The automated inventory probes format accessibility only. Neither published proteomic protein matrices nor the article figures automatically establish a verified **donor ID ↔ compound ↔ system ↔ CLint ↔ protein abundance** join.

### Prespecified first modeling protocol (strictly conditional)

**Biological unit:** 15 independent liver donors, NOT 75 unrelated individuals (five correlated compounds measured within donors).

**Endpoint candidate:** measured donor-compound CLint in **hepatocytes**, log-transformed with original assay units and censoring rules recorded, **not** patient CL, CL/F, AUC or MIPD dose. Keep the HLM and HH systems separate.

**Prediction-time design:** distinguish proteomics-only prediction from HLM-to-HH assay transfer (where donor HLM CLint might be a legitimate predictor). Do not mix them or present an HLM measurement as a molecular-only predictor.

**Outer split:** leave-one-DONOR-out, holding all five compounds and tissue fractions for that donor out of all training-stage preprocessing, feature selection, and model tuning.

**Baselines:** fixed drug-specific training-only mean/log-median and a small predeclared ridge baseline with directly related CYP proteins. Benchmark more complex approaches only after n=15 sample constraints are explicitly analyzed. XGBoost with full proteome is not justified as a primary baseline.

**Experimental dependencies:** retain original normalization, CYP assay fraction, units, protein/cell factors and measurement status. Kp_uu is measured experimental intracellular accessibility, not always a prospectively available predictor; evaluate separately as a mechanistic sensitivity.

**Metrics:** donor-clustered MAE/log fold-error, bias, by-compound errors and descriptive donor bootstrap intervals. No pseudo-75-patient inference, feature network fitted with test donors, or outcome-selected hyperparameters.

### Mandatory go/no-go

1. Official supplement contents can actually be retrieved, with immutable hashes.
2. Original donor-level CLint and proteomics have verified matching identifiers and compatible sample fractions.
3. Clear extraction of original assay units, missingness/censoring/replicate rules and labels.
4. If donors cannot be joined without digitizing figures or guessing aliases, **do not fit a model**. Obtain tabular primary CLint data from the authors or use another source.

The first source availability job does not complete any of these target-pairing gates merely because an XLSX can be downloaded.

## 2. Independent empirical SWEET validation (Yin/Battail 2025)

Yin L et al. "Sample-specific network analysis identifies gene co-expression patterns of immunotherapy response in clear cell renal cell carcinoma". iScience 2025, DOI [10.1016/j.isci.2025.113061](https://doi.org/10.1016/j.isci.2025.113061), PMCID PMC12309972.

The study describes 309 advanced ccRCC patients from CheckMate 009, 010 and 025, separated by nivolumab/everolimus and primary tumor/metastasis. The authors identify [Zenodo 15723817](https://doi.org/10.5281/zenodo.15723817) as the data/code archive.

This is real treatment response (PD-related), NOT drug exposure, individual CLint, plasma PK, or dosing. The source job only verifies Zenodo record metadata and inventories file names.

Before replication: verify exact archive/version/license, available outcome labels, patient identity and treatment/biopsy partitions; reproduce original SWEET descriptively without mislabeling the **transductive** calculation as out-of-sample. Compare locked expression-only and strictly inductive SWEET networks using held-out patient/biopsy-group splits, training-only preprocessing, no outcome-derived test-network references; handle survival censoring rigorously. Preserve separate experimental upstream z-score discrepancy tracked in issue #4.

## Audit software interpretation

The GitHub Actions workflow:
- Runs **offline mock tests** for valid XLSX workbook names, PDF signature, Zenodo ID, bad upstream HTML responses and false "clinical data validated" assertions.
- Uses bounded explicit source URLs, never GitHub raw data or user/patient content. Outputs only format-status metadata.
- On inaccessible sources, records UNAVAILABLE_OR_INVALID; successful CI must never be interpreted as confirmed accessible files unless JSON source states are checked.
- Never uploads raw biological/proteomic records or supplementary file bytes to any artifact.

**Next engineering target:** strict local-only donor/drug/system/CLint-unit/protein-matrix crosswalk manifest with synthetic anti-leakage tests. Actual numerical CLint models are blocked until a genuine donor-outcome linkage is demonstrated.
