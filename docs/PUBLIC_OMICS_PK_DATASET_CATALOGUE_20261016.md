# Inventory of public **omics × real PK/metabolism** datasets (verified 2026-10-16)

**Project objective:** compare a reproducible molecular phenotype (RNA alone; protein alone; metabolite alone; qualified combinations) against measured drug exposure, intrinsic clearance or enzyme function. Preserve **exact same-donor or same-patient identifiers, source units and chronology**. Cytotoxicity (GI50/IC50), transplant outcomes and treatment-response classifications are NOT PK endpoints.

**Owner-confirmed constraint:** BIOMARGIN has **no PK data**. A theoretical ABIS pairing is only a permissioned *future feasibility question*, not a public benchmark or an established crosswalk. See issue #14.

## A. Numerically published, real in-vitro metabolism plus molecular phenotype

### A1. Tucatinib — Fallon/Jackson et al., 2025 — first actionable proof-of-concept candidate
- Original paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC13095627/ ; DOI [10.1016/j.dmd.2025.100061](https://doi.org/10.1016/j.dmd.2025.100061).
- Proteomics: https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD057282 .
- **21 original single-donor human liver microsomal (HLM) preparations and ONE 150-donor pooled reference**, not 22 independent human donors. Table 1 *of the actual public article* provides individual microsomal tucatinib clearance (uL/min/mg protein), CYP2C8/CYP3A4/CYP3A5 protein concentrations (pmol/mg protein), midazolam CYP3A activity and amodiaquine CYP2C8 activity in donor-aligned rows.
- **Direct data availability:** donor-level numbers in the published Table 1; raw targeted quantitative proteomics underlying measurements in PRIDE PXD057282. Source rights/re-use and definitive correspondence of raw files to donor labels must still be checked before redistribution.
- **Scope:** targeted measured protein → measured in-vitro clearance of ONE molecule; no transcriptome, no clinical C(t) or plasma CL/F. Because published analysis already fitted CYP models, reproduce with **pre-specified leave-one-donor-out** and honest performance, not claim discovery or external clinical validation.
- **Scientific gates:** distinguish CYP3A total vs CYP3A4+CYP3A5, activity vs protein, pooled reference vs individual donors; compare training-only mean, fixed targeted-CYP ridge and activity-only+protein sensitivity on identical donor folds; n=21, no high-dimensional models.

### A2. Wegler et al., 2021 — cross-substrate proteomics/CLint, source verified
- Paper and three supplements: https://pmc.ncbi.nlm.nih.gov/articles/PMC8041379/ DOI 10.1021/acs.molpharmaceut.1c00053.
- PRIDE: https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD014131 .
- **15 biological human liver donors, 5 probe compounds, donor-matched hepatocytes (HH) and liver microsomes (HLM), global proteomics.** Five drug CLint observations per donor are repeated measurements, NOT five independent people.
- Original PMC AWS PDF and protein XLSX successfully downloaded and SHA-256 verified in Omics2PK CI PR #10/#13. Supporting PDF has 21 native-text pages; proteomics Data S1 contains separate Hepatocytes and Microsomes sheets.
- **Outstanding:** the exact donor × compound × preparation × original numeric CLint matrix has not been verified as an extractable table/paired outcome. A figure caption or protein file is NOT proof of tabulated clearance. Issue #11 blocks model fitting until that crosswalk is demonstrated.
- **Scope:** real in-vitro clearance, not clinical PopPK. Cannot evaluate SWEET directly without RNA.

## B. Transcriptomic bridge with observed enzyme metabolism (not yet clinical PK)

### B1. Peng et al. 2025 — **6 PHH preparations** + single-cell transcriptomics and CYP substrate disappearance/metabolites
- Official public RNA: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289636 (six GEO samples; supplementary CSV only cell-level Seurat metadata, NOT RNA gene counts).
- Article and Supplement 1 (Tables S1–S7 including S4/S5): https://pmc.ncbi.nlm.nih.gov/articles/PMC12532551/ ; https://doi.org/10.1186/s13619-025-00258-6 .
- Same six primary human hepatocyte batches examined by scRNA-seq and LC-MS/MS activity of six CYP substrate/metabolite pairs, including bupropion/CYP2B6, diclofenac/CYP2C9, testosterone/CYP3A4; **substrate clearance rate (SCR) / metabolite formation rate (MFR)** (laboratory metabolic function, NOT in-vivo CL or AUC). Original DOCX Tables S4/S5 verified as donor-specific six-CYP functional tables with exactly the six GEO batch IDs. Gene expression counts are NOT supplied in GEO CSV; functional units, censoring and LoQ still need manual audit.
- **Important:** author data availability states *other original data available on reasonable request*. Do not assert full individualized analytical matrices outside GEO/Supplement 1 are public. N=6 is proof of data-flow/mechanism, **not** a generalizable machine-learning model or substantive external validation.
- **ABSTAIN current:** GEO files contain Barcode, Seurat clustering, UMAP, PCA, nGene/nUMI rather than CYP expression. First recover actual gene matrices from SRA/FASTQ or authors. Then use biological donor CV, not cell-level CV.

### B2. Yang et al. 2010 — human liver cohort, RNA expression + multiple CYP activities
- Paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC2909567/ ; https://pubmed.ncbi.nlm.nih.gov/20538623/.
- Public liver microarray from related cohort: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9588 (**427** publicly profiled GEO liver samples).
- Original integrated study characterized **466** human liver samples including gene expression, SNPs and measured CYP enzymatic activity, but 466 integrated study records and 427 GEO microarray records are **not automatically the same fully overlapping participant set**. GSE9588 explicitly says clinical data were requested but **not provided** by submitter. HLC: https://www.synapse.org/Synapse:syn4499 ; linked restricted genetic accession https://www.ega-archive.org/studies/phs000253 .
- **Outcome-access status:** 398 CYP activity-measured subjects are reported within the 466-subject study, and an unrestricted per-sample quantitative CYP-activity table with exactly matching GEO identifiers is NOT yet confirmed. It may require dataset-specific Synapse/dbGaP permission or original supplemental data. Do not fit against paper-level averages or infer matching identifiers.
- **Potential value:** much larger RNA→functional drug metabolism demonstration, with batch/site confounding and genotype baseline if legally accessible.

## C. Actual patient drug concentrations plus molecular measurement, **but individual paired public matrices not confirmed**

### C1. Tacrolimus, 137 liver recipients / 693 plasma samples — untargeted metabolomics with measured trough/dose
- Article: https://pmc.ncbi.nlm.nih.gov/articles/PMC9659571/ ; DOI https://doi.org/10.3389/fphar.2022.1045843 .
- Paper documents same-study **693** matched plasma profiles from **137** liver recipients, tacrolimus whole-blood C0, dosing and clinical context, and a longitudinal mixed model. This is a serious clinical molecular–PK candidate in principle.
- **No confirmed public sample-level matrix linking the full untargeted intensities, C0, dose, repeated subject ID and sampling timestamp.** Paper's data statement directs requests to authors. The journal article, figures and identified metabolites alone are insufficient for external held-out validation.
- Request authorized **deidentified donor-disjoint** molecular feature / concentration / dose and specimen-time matrices; use repeated-measures split by 137 people, not 693 observations. C0/dose prediction is not equivalent to full concentration-time PK model.

### C2. Genomics + real tacrolimus TDM: GEO and dbGaP (genotype/PGx baseline, NOT transcriptomics)
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130068 : public genomic variation array on 115 liver donor/recipient pairs, source paper also measures dose-normalized tacrolimus C0 in weeks 1–4; **same-person granular C0/dose arrays not confirmed in public GEO sample metadata**.
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196772 : SNP genotyping for 251 kidney recipients whose study assessed tacrolimus C0; granular PK case table availability requires checking.
- https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001670.v1.p1 : DeKAF Genomics 48 selected extreme tacrolimus-metabolism recipients; **authorized-access dbGaP**, not openly downloadable subjects+outcomes. This is a *genomic+TDM* restricted candidate, not open RNA.
- These sources can inform the future **clinical+PGx reference**. They cannot establish incremental transcriptomic or metabolomic value without those modalities and subject-level join.

## Deliberately excluded as quantitative PK evidence
- NCI-60 doxorubicin GI50 RNA/metabolites, although measured, is drug susceptibility/response, NOT PK.
- Tacrolimus-treated podocyte RNA GEO GSE317483: transcriptomic drug response under a declared concentration, NOT measured metabolism/PK.
- Tacrolimus withdrawal renal recipient RNA GEO GSE145408: immune response/PD/outcome without verified PK concentration-time curves.
- BIOMARGIN alone: multiomics/biomarkers but **no PK**; ABIS matching unproven.

## Ranked engineering/evidence actions
1. **P0 source/rights gate**: obtain the legitimately usable original 21-donor Table 1 numbers for the tucatinib paper; publisher automated access returned 403 and text/data mining rights are reserved. Local-only source contract and cross-validation code are implemented with synthetic tests; validate source units and pooling; fixed leave-one-original-donor-out protein→microsomal CLint benchmark with very small predeclared model set and full negative-result reporting. This is a **targeted proteomics** experiment, NOT thesis proof for transcriptomics or clinical individualized dosing.
2. **P0 transcriptomic block**: donor-matched six-CYP SCR/MFR panels confirmed, but GEO download is cell metadata only. Retrieve original SRA reads and process to gene counts, or request usable gene expression matrices, before any RNA-CYP association; n=6 remains exploratory.
3. **P1 wider transcriptomic**: qualify per-donor CYP activity outcomes and lawful exact join for 2010 HLC/GSE9588; refuse if no open outcome file.
4. **P1 separate PK**: source-audit the original Wegler 15-donor × five-substrate HH/HLM CLint supplement, and obtain author data if not tabulated.
5. **P2 clinical**: contact 2022 tacrolimus metabolomics study authors for matrices and/or verify feasibility of lawful BIOMARGIN molecular ↔ ABIS drug dose/concentration linkage; do not fabricate or fuzzy-match subjects.
6. **P3 site**: only after validated empirical proof and source/legal checks, expose static safe research lessons and tools at rberrah.github.io with FR/EN, dark/light, Resources and Other Disciplines; no patient-level uploads or public model-driven dosing.

**Audit certainty:** literature and official dataset identity/public download statements checked as of 2026-10-16; a published paper describing a paired assay is *not* proof every raw feature/outcome join is unrestricted and accessible. A green source-format CI is not a green scientific validation.


**Independent post-audit qualification, 16 October 2026:** [Source audit and first validation gates](SOURCE_AUDIT_AND_FIRST_VALIDATION_GATES_20261016.md). A green source/test CI does not mean the measured empirical model was fitted. Current RNA pilot abstains; tucatinib code awaits legitimate original numerical source data.
