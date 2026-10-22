# Source audit and first scientific validation gates — 16 October 2026

**Updated result (22 Oct 2026):** primary molecular/functional provenance is confirmed and several genuine held-donor E2 benchmarks have now been executed. Wegler 2021 is source-qualified for targeted CYP protein → measured in-vitro CLint; EPA RNA→CYP activity and tucatinib protein/activity→CLint are also qualified elsewhere in the repository. These remain **preclinical functional metabolism**, not clinical plasma PK.

## A. Peng 2025 — original transcriptome and six measured donor CYP panels

Article https://pmc.ncbi.nlm.nih.gov/articles/PMC12532551/ DOI 10.1186/s13619-025-00258-6. GEO https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289636 ; BioProject PRJNA1223864.

- PMC AWS official, CC BY, original Supplement 1 DOCX 13619_2025_258_MOESM1_ESM.docx, SHA256 deeb835fe82d7e0e883447268354b9abb8f3ca4950639be3b68b802d3c6183bf.
- **TRUE** original donor-specific two functional panels: substrate clearance and metabolite formation for CYP1A2, CYP2B6, CYP2C9, CYP2C19, CYP2D6, CYP3A4. Each 7 × 8 table lists all six PHH batches 025, 211, 330, 409, 416, 789 with original measured functional assay outcomes. Technical repeat/LoQ, units and sampling time remain methodological gates.
- GEO SOFT exactly six donor-matched samples; SHA256 443635f9a9c527ee32b1e0d3559f91b2583304c36e78cba5c972d9681a72012a.
- GEO six-file RAW.tar SHA256 7f114365ba5e977e685f9cee108def64deb0a31c96ea13d95e19c0b42c2b01e8.
- **CRITICAL audit correction:** the file names include "row_counts", but their actual columns are Seurat **cell-level metadata**, not gene count data: Barcode, graphclust, kmeans clusters, TSNE/UMAP, PC1–PC10, nGene/nUMI, mitochondrial QC. Matrices have 33–35 variables and around 8k–25k cellular rows per donor, with negative fractional PC/UMAP coordinates. **None of the six CYP gene expression features are supplied.** These data MUST NOT be normalized as UMI counts or passed to a model.
- Source DOI/GEO refer to actual scRNA experiments, but the SRA raw sequence reads or author-provided true processed gene-count matrices are REQUIRED to reconstruct CYP expression.
- The protected modeling code correctly emits ABSTAIN_GEO_HAS_SEURAT_METADATA_NOT_EXPRESSION_COUNTS and does not fit any real patient/donor RNA→activity association.
- Source n = **six original human donors**, not thousands of independent observations. Even after obtaining matrices, this is descriptive/exploratory in-vitro CYP function, not plasma PK.

GitHub Actions run of original source inspection: https://github.com/rberrah/Omics2PK/pull/18 .

## B. Tucatinib 2025 — 21 original donor clearance outcomes and quantitative CYP proteomics

Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC13095627/ DOI 10.1016/j.dmd.2025.100061; proteomics PRIDE PXD057282 https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD057282 .

- Paper Table 1 confirms **21 individual-donor** HLM observations of tucatinib apparent metabolic clearance with original CYP2C8/CYP3A proteins/activity and one **150-donor pooled reference preparation** (do NOT use the pool as donor #22).
- PRIDE holds targeted quantitative proteomics/peak-area raw data, **not** the original indexed individual-dose clinical data nor necessarily an accessible tabular clearance data file.
- Official publisher's site explicitly reserves rights (including text/data mining) and automated access to its Table 1 returned 403. EPMC fullTextXML returned 500, PMC OA AWS lookup for first version returned 404. Never bypass automated retrieval restrictions or commit publisher's tabular values without a rights/provenance review.
- The project now contains local-only, original 21-donor Table 1 schema and actual LODO fitting machinery, with a fixed protein CYP2C8+CYP3A pair, activity comparator, mean baseline, train-only standardization, fixed alpha=1 ridge, donor-paired bootstrap. **Only synthetic tests run.** No original Table 1 rows obtained/fit yet; no biological performance estimate is reported.
- Authorized command after local curation: python -m scripts.validate_local_tucatinib /authorized/source/table1.csv . Output is aggregate only; source table stays local.
- Other imatinib HLM 2025/2026 article reuses the same donor CYP phenotypes: useful **cross-substrate** challenge but NOT external donor-independent validation.

## C. Remaining cohorts/claims

- Wegler 2021 original HH/HLM samples, n=15 donors, five drugs: **source-qualified and fitted** from native Table S3. Exact same 15 donors are retained across all five drug blocks; targeted CYP protein and raw HLM/HH CLint are paired from the original table. HH shows heterogeneous positive signal (midazolam/omeprazole/diclofenac), while HLM is largely null. Global Data S1 high-dimensional proteomics is not required for this targeted result.
- Human Liver Cohort, GSE9588: 466 molecular study liver samples, of which **398 have reported enzyme activity**, 427 genotyped GEO expression records; do not silently equate all these cohorts. GEO expressly reports **clinical data not provided by submitter**, so exact indexed CYP activity↔RNA matching still needs a data access audit.
- Tacrolimus clinical metabolomics, 137 liver recipients/693 samples: genuine reported C0/dose + metabolomics in paper, but independently downloadable paired individual matrix not confirmed. Do not claim a held-out PK benchmark.
- BIOMARGIN **has no PK**; potential authorized ABIS person+time+data matching is a hypothetical feasibility assessment, not a verified dataset.

## Scientific decision gates

| Stage | Observed result | Meaning |
|---|---|---|
| E0 publication, legal/format/file identity | PASS for Peng, partial for tucatinib | Source files verified, not just abstracts |
| E1 functional CYP donor–RNA assay linkage | PARTIAL | Six exact donor IDs in two assay tables; corresponding feature **expression** absent |
| E2 RNA-to-CYP measured functional benchmark | ABSTAIN | Source GEO CSV is Seurat metadata, not gene expression |
| E2 CYP protein-to-measured in-vitro clearance | **PASS, heterogeneous** | Wegler n15×5 source-qualified LODO plus separate tucatinib evidence; positive signal is system/substrate dependent, not universal |
| E3 genuine clinical omics+PK model | NOT VERIFIED | Still need real identical patient × assay × dose × timed concentration records |
| /omics2pk/ public scientific tool | NOT DEPLOYED | Must await eligible empirical source results and rights review |

**Do not confuse software CI PASS with verified model improvement.** Null, blocked, and ABSTAIN outcomes are valid scientific audit results and must remain visible.
