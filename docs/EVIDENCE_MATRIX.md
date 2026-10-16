# Qualification evidence matrix — Omics2PK (2026-10-16)

All quantitative performance claims require a recorded observation origin,
sample/assay identity, prospective-vs-transductive method status and independent
validation. Do not combine distinct evidence levels.

| Dataset | Biological unit | Molecular inputs | Observed pharmacology | Scientifically valid conclusion |
|---|---|---|---|---|
| Official [SWEET example](https://github.com/SysMednet/SWEET) | 80 preprocessed expression samples, TCGA | RNA expression, 100 genes | **None** | Numerical reproduction, stability/inductive procedure only |
| [mixOmics multidrug](https://github.com/mixOmics-org/mixOmics) | 60 NCI-60 cancer cell lines | 48 ABC transporter RT-qPCR | Measured *in vitro* drug GI50 (1,429 compounds) | RNA-derived network feature versus simple baseline for **cell-line drug response** |
| [Yin, Battail et al., iScience 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12309972/) | Patients / tumors, subject to verified sample-level labels and approvals | Tumor expression, sample-specific networks | Reported immunotherapy response | Independent treatment response benchmark if downloaded samples match outcomes; **not drug PK** |
| [Wegler et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8041379/) | 15 human liver donors, documented paired measurements | Liver proteomics (not transcriptomics) | Measured in vitro hepatic intrinsic clearance | Protein → CLint study; no SWEET RNA inference without measured RNA |
| Local BIOMARGIN / ABIS–ISBA | Individual transplant recipients only if matched and authorized | Cohort omics + possible clinical/PGx | May include dose/concentration information in separate records | True PopPK only after patient/time reconciliation and ethics/governance checks |

**Specific paired-data pitfalls:**
- A real response to anticancer treatment (PD) is **not** tacrolimus clearance, bioavailability, or clinical AUC.
- QPCR transporter expression is a molecular state measurement, **not a DNA genotype**, and cannot be rebranded as expression for other metabolic enzymes.
- In liver-donor proteomics, 15 donors are too few for high-dimensional SWEET and any generalized M4 graph claims.
- Clinical trough-only tacrolimus PK typically identifies apparent CL/F rather than separating CL and oral bioavailability F without external information.
- Repeated measurements and technical replicates must not be treated as independent patients.
- Omics measurements and PK observations may not be contemporaneous. Record tissue, treatment start, sampling time, co-medication and disease state.
- Molecular state can reflect a pharmacological response; post-dose measurements cannot be presented as pre-dose predictors without temporal qualification.

## Ordered milestones

**Gate D0 – Public-code reproducibility.** Omics2PK imported MultiOmics V1 source producer on synthetic three-omics; official SWEET weights/raw-edge reproduction. Source example z-score archive mismatch recorded in issue #4. PASS for source compatibility only.

**Gate D1 – Paired measured PD.** Prespecify a drug, use exact independent cell lines and outer-fold-only network references. Compare intercept/raw expression/PCA/network/M1/M1bis/M2 and shuffled graphs using common splits. Publish all results, including no improvements. Doxorubicin is frozen before outcome-model execution, and other compounds must be prespecified independently before evaluation.

**Gate D2 – Multi-omics.** True matched samples of RNA+protein+metabolites and response; train-only MOFA2/other factor fits; control missingness and additional sample attrition. NCI-60 IntLIM expression+metabolites has *no drug outcome in its sample metadata* by default; linking GI50 requires a documented line-name crosswalk, not mere sample-order matching.

**Gate D3 – Mechanistic score M3.** Predeclare the drug/tissue evidence priors from Pillar A; model scores as pharmacokinetic covariates in actual NLME likelihood. Compare genotype alone, phenotype molecular, multiomics, and with/without Pillar A. Include mechanism-irrelevant network controls and a prior-only module.

**Gate D4 – Graph-PopPK M4.** Implement a parameter-level graph response (e.g. CLint mechanisms) in an identifiable NLME model, not merely a graph embedding put into an arbitrary regressor. Report identifiability, reference/scale stability and cost.

**Gate D5 – External clinical validation.** Real observations, independent locked cohort, dose histories, formulations, concomitant covariates; robust sparse-TDM posterior predictive checks (C0 and C0+C1), exposure metrics, bias, uncertainty and subgroup shifts. No validated claim until achieved.

## Intellectual property / cloud

Original Omics2PK code is licensed as described at repository root.
Third-party mixOmics datasets and any downloaded components retain their own
original license and attribution (mixOmics software uses GPL >=2).
No source file, data or license is overwritten here.
Only publicly available aggregates are uploaded by public CI. Controlled
omics, patient-level exposure data and full local R/Julia research analyses
remain local/institutionally authorized.
