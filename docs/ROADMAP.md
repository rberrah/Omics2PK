# Omics2PK — working to-do list (16 October 2026)

This list follows the earlier MultiOmics → SWEET → Pillar A → Pillar B R2 discussions. **When inspected the Omics2PK repository was empty; there was no pre-existing to-do file to read.** The plan below has therefore been reconstructed and must be reconciled with the owner's priorities.

## P0 — foundation and executable contract (initial PR)
- [x] Define research-only, no-patient-data boundary and M0–M4 scientific comparison.
- [x] Validate native MultiOmics DOWNSTREAM_EXCHANGE_V1 locally (SHA, patient IDs, features, missingness and omic coverage).
- [x] Prepare original SWEET input format from normalized RNA only; never infer networks without running a qualified method.
- [x] Add synthetic-only unit tests and offline CI.
- [x] Run the **actual MultiOmics JavaScript export producer** using its three bundled synthetic demo omics and validate the generated exchange with Omics2PK. PR #5; 8 longitudinal subjects, 49 assays, one technical RNA repeat; not yet a browser ZIP end-to-end test.
- [ ] Run an actual user-downloaded portable browser ZIP end to end (separate from the validated native JavaScript producer); avoid uploading any user or clinical data.
- [x] Propose explicit PolyForm NC / CC BY-NC-SA licensing on the draft PR, with attribution and scope/exceptions.
- [ ] Confirm institutional/coauthor rights before redistributing any third-party code. Upstream SWEET sources are not vendored.

## P1 — reproducibility and SWEET
- [x] Pin the upstream SWEET SHA (`cf064eb5f5b996c40bba3cce507384f6c0fdecb5`) in documentation; upstream code is **not** redistributed.
- [x] Implement strict SWEET edge-output reader and tabular **M1bis feature export** with synthetic tests (not network inference).
- [x] Reproduce official SWEET sample weights (80/80), one full individual raw edge set (4,753/4,753), and global summary. **Archived z-scores differ (4,753/4,753); tracked in issue #4; never claim complete equivalence.**
- [ ] Audit SWEET handling of zero variance, zeros, missing values, and the -z flag.
- [x] Prototype a clearly labeled **non-equivalent SWEET-inspired** training-frozen reference and query-only network (PR #3).
- [x] Add subject-disjoint training leave-one-out/cross-fitting and outer train/test feature construction, with negative leakage tests (PR #3, merged).
- [x] Separate genome-wide sample weights (100 official example genes) from a targeted edge panel (25 fixed-order genes), and compare against the official transductive equation. Public diagnostic: 60 train/20 held-out, one clipped weight, ~27.1% relative Frobenius difference for one test network. No outcomes or clinical validation.
- [ ] Validate a genuine out-of-sample extension mathematically and compare against published SWEET; current prototype uses selected genes for sample weighting instead of the genome-wide input.
- [ ] Compare stability to sample count, preprocessing, batch, tissue and threshold selection. Initial 60/20 public expression-only diagnostic is in merged PR #3, but no PK endpoints, no cross-fold preprocessing calibration, and no biological claims.

## Thesis Chapter 4/10 -> genuine PK proof: cross-domain gate
- [x] Pin the **15 October 2026 frozen thesis architecture** as a public research mapping without exposing potentially patent-sensitive Pillars A/B/C/D source (docs/OMICS_TO_PK_THESIS_PROOF_AND_SITE_PLAN.md).
- [x] Implement strict baseline omic-to-observed-PK patient/date/dose/concentration metadata check: RNA-only or 2–3 omics; time zones, prior-to-dose sampling and fixed assay alignment gap; private data never uploaded.
- [x] Audit source file structure for Wegler 2021 PDF and donor-matched proteome workbooks, using exact official public SHA-256 and output only aggregate metadata. This is NOT a quantitative CLint table join.
- [ ] **P0 science blocker**: identify genuine donor-level Wegler CLint outcomes and match to donor-proteome rows with source table/assay/units provenance (issue #11). If not publicly tabulated, STOP and ask authors for numbers.
- [ ] First measured in-vitro donor-held-out CLint regression with fixed protein CYP scores, no pseudoreplication; report failure if no verifiable labels.
- [ ] Separate **transcriptomics→measured PK** source (or functional assay), not just phenotype→tumor response.
- [ ] Actual clinical BIOMARGIN/authorized linked tacrolimus or MPA TDM cohort with patient-date/dose/concentration links (ethics/license gate).
- [ ] Compare clinical+PGx reference → prespecified direct RNA scores → qualified SWEET RNA → parallel protein/metabolite → optional integrated latent scores, on identical held-out patients; external validation + Pillar C adjudication.
- [ ] **After proof only**: design public learning/tool at rberrah.github.io with FR/EN, dark/light, Resources, Other Disciplines, persistent provenance and no user clinical-data upload.

## E2 REAL original measured biological phenotype→drug-metabolism proof (2026-10-16)
- [x] Tucatinib: n=21 independent original HLM donors, published source Table 1 manually transcribed and **all four published Spearman CYP↔CL correlations reproduced**; 21 leave-donor-out folds MAE mean 6.246 vs 2 CYP proteins 2.041 vs 2 CYP enzyme activities 1.940 µL/min/mg, explicitly nonclinical and original single cohort. Corrected source narrative/table inconsistency identified; donor measurements not redistributed. Detailed [tucatinib proof](TUCATINIB_21_DONOR_MEASURED_VALIDATION_20261016.md).
- [x] **RNA/PK bridge now has genuine measured data:** Liu et al US EPA original 100 human livers, exact match across two original sample IDs, n=94–100 complete paired donors in 8 transcript→microsomal CYP activity assays, genuine donor-held-out train/test runs on the original source in GitHub CI. Compared age+sex versus age+sex+RNA, 7/8 lower MAE, 4/8 uncorrected paired bootstrap intervals exclude 0. **CYP2E1 worsened**. Source/licensing/data provenance in [EPA RNA-CYP proof](EPA_100_DONOR_TRANSCRIPT_CYP_VALIDATION_20261016.md).
- [x] Independent numerical replication: pinned scikit-learn StandardScaler/Ridge used on actual original n100 liver source in a distinct script; all **8 CYP held-donor-out MAEs match the native solver to six decimals**, GitHub Actions https://github.com/rberrah/Omics2PK/actions/runs/38930709329. This is independent *estimator* parity, NOT independent donor/site external validation.
- [ ] R/Julia third-party replication beyond sklearn, and multiple-testing-aware inference, multiple-testing-aware inference, sensitivity to CYP outliers, donor age/batch, and original paper processing consistency.
- [ ] Externally replicate with separate donor cohort and different drugs/assays. E2 outcomes are *in-vitro functional metabolism*, NOT systemic clearance, measured plasma concentration or patient MIPD.
- [ ] True E3 **clinical** Omics→PopPK remains blocked until compatible public or permissioned PK/TDM molecular data contain actual dose/time/concentration for the same biological individual. Clinical+PGx PopPK required as true comparator.
- [ ] Only after review/rights qualification, expose FR/EN beginner courses and a clearly labelled **nonclinical** research demonstration at rberrah.github.io. Never use published in-vitro phenotype to suggest patient dosing.

## Measured-source audit — real source bytes read, no fake model fit
- [x] **GSE289636** actual GEO SOFT/TAR and original CC BY source DOCX: six exact PHH donors, two six-CYP functional outcome tables. SHA-pinned.
- [x] **Critical discovery:** files labelled row_counts are Seurat CELL METADATA (Barcode, graphclust, PCA, UMAP, nUMI), not any CYP gene transcript matrix. Explicit **ABSTAIN** instead of false CPM-based regression; obtain raw SRA FASTQ-derived gene counts or author matrices (issue #16).
- [x] Original tucatinib 21 single-donor+one pooled HLM paper/PRIDE provenance confirmed. Automated publisher Table 1 access forbidden/reserved; numerical donor records NOT downloaded or modeled.
- [x] Local-only 21-donor HLM clearance–CYP protein/activity LODO evaluator with frozen 2-variable ridge, baseline and paired bootstrap. **Synthetic tests only**, awaiting rights-qualified source CSV (issue #15).
- [ ] Run first **real** source-labeled molecular predictor/functional PK endpoint dataset only after exact numeric outcome+feature pairing; distinguish source-format CI PASS from biological model validation.
- [ ] For HLC GSE9588, confirm exactly matching enzyme activity records in the 398 activity-profiled subset of the 466 molecular subjects; no automatic 427 GEO overlap.
- [ ] Site learning bridge remains withheld until at least a reproducible matched quantitative result. Full audit: [SOURCE_AUDIT_AND_FIRST_VALIDATION_GATES_20261016.md](SOURCE_AUDIT_AND_FIRST_VALIDATION_GATES_20261016.md).

## 2026-10-16 public evidence inventory
- [x] Verify accessible article/dataset identity for a **21-donor tucatinib CLint + targeted CYP proteomics table** (original public Table 1, PXD057282). This now precedes Wegler as the first direct numerical CLint validation candidate.
- [x] Identify **GSE289636 (six donor PHH scRNA with paper-supplement CYP substrate clearance)** as a separate *RNA→measured functional activity* proof, not clinical PK; source numeric assay crosswalk remains to verify.
- [x] Qualify **GSE9588/HLC** as RNA+enzyme activity study (466 original liver cohort vs 427 GEO expression records), but sample-level activity join not public-confirmed; no unverified merge.
- [x] Identify **137 transplant recipients / 693 plasma samples** with clinical tacrolimus C0 and metabolomics (paper), but original per-sample combined matrices **not confirmed downloadable**, require author approval.
- [x] Record 115 liver paired and 251 renal recipient SNP GEO series and DeKAF restricted dbGaP pharmacogenomic+TDM sources as genomic baselines, NOT transcriptomics.
- [x] Incorporate user-confirmed fact: **BIOMARGIN has no PK**. An ABIS join is only an unconfirmed authorized possibility, not a declared dataset.
- [ ] Start source-provenance, held-out exploratory tucatinib 21-donor proof; do not count pooled HLM as donor.
- [ ] Qualify original numeric substrate clearance tables for six donor scRNA PHH, before fitting any RNA activity model.
- [ ] Independently check access and per-donor enzyme assays for the much larger HLC/GSE9588 RNA cohort.

Full links, licenses and scientific restrictions: [public Omics×PK data catalogue](PUBLIC_OMICS_PK_DATASET_CATALOGUE_20261016.md).

## P2 — public measured-data pilots
- [x] Public pharmacological pilot: mixOmics NCI-60 60 cell lines, 48 ABC transporters and **measured doxorubicin GI50** (chosen before predictions); M1/M1bis/M2 plus expression/PCA/Ridge, negative controls and histology adjustment. **No MAE improvement over mean or tumor type**; report in docs/RESULTS_NCI60_DOXORUBICIN_20261016.md.
- [x] Public NCI-60 cross-study RNA + metabolite GI50 pilot: 48 syntax-only matches, **47 histology-concordant** lines (MDA-MB-435 excluded). 17,987 RNA/280 metabolites. No incremental metabolomic MAE gain across fixed PCA/Ridge or XGBoost on matched lines. Research-only merged PR #9 and docs/RESULTS_NCI60_TWOOMICS_20261016.md. Not same-sample multiomics or clinical PK.
- [ ] Sensitivity to full-cohort upstream preprocessing, tumor-stratified held-out splits, feature caps and bootstrap uncertainty for the cross-study pilot.
- [ ] Yin/Battail: independently reproduce patient-specific network + response analysis (public Zenodo metadata source format audit created). A measured response, **not measured PK**. See docs/PUBLIC_PK_SOURCE_QUALIFICATION.md.
- [ ] Wegler et al.: qualify verifiable per-donor protein and measured CLint joins from official PMC supplements (source-format inventory implemented) before any fit. 15 original donors, five probe drugs; **no SWEET without RNA**. See docs/PUBLIC_PK_SOURCE_QUALIFICATION.md.
- [ ] Search additional public paired RNA+PK and RNA+PD cohorts; inventory licensing, IDs and dates.
- [ ] Ask P&T whether existing **ABIS/ISBA** clinical PK records can be *legally* joined to BIOMARGIN omics for the same recipients, with appropriate timings. Do not claim BIOMARGIN itself supplies PK.

## P3 — feature representations and comparisons
- [ ] Clinical + PGx original PopPK reference and direct expression/PCA/GSVA/MOFA2 comparators.
- [x] Pilot M1 (raw ABC expression), M1bis (experimental network), M2 (raw+network) using real NCI-60 doxorubicin GI50, 5-fold held-out and 3×5 repeated sensitivity. **All worse MAE than a tumor-class-only baseline**; not clinical or multimodal.
- [ ] Upgrade to stratified/residualized confounder controls and nested hyperparameter checks; no post-hoc positive claims.
- [ ] XGBoost full **multiomics** (M1), SWEET-only (M1bis), omics + SWEET (M2) on genuinely paired patients with verified modalities and full R/clinical outcome; not yet possible from 48 ABC qPCR variables alone.
- [ ] Low-dimensional SWEET scores without Pillar A vs drug/tissue knowledge-guided scores with Pillar A (M3).
- [ ] Graph kernels/low-dimensional graph embeddings with and without Pillar A (M4 exploratory).
- [ ] Negative controls: nonmatched subject/network, shuffled graph edges, off-drug mechanisms, null outcomes.
- [ ] Keep DNA genotype, RNA phenotype and multiomic phenotype distinct. Penalize redundant effects with existing genotype-based covariates.

## P4 — genuine PK/PD inference
- [ ] Read-only evidence interface from private Thesis Pillar A; **do not copy private code or full knowledge vault**.
- [ ] Pre-specify biological mechanism-to-PK-parameter mappings and uncertainty.
- [ ] Export numerical candidates to R2 and fit **actual NLME/SAEM on observed concentrations and dosing**; oral data identify CL/F under routine assumptions.
- [ ] Separate no-TDM prediction from C0/MAP Bayesian personalization; avoid endpoint leakage.

## P5 — frozen validation, delivery
- [ ] Freeze training-only gene selection, transformations, SWEET reference, graph construction and hyperparameters.
- [ ] Evaluate against locked external cohort; compute uncertainty, calibration, and subgroup transportability.
- [ ] Replicate reproducible open example and full manuscript methods; never claim clinical validation from synthetic runs.
- [ ] Add opt-in MultiOmics page for *qualified* public demos, via already established rberrah.github.io URL architecture. Static site only, local/authorized compute for real analyses.

### Go / no-go criteria
G0 public-code ownership/privacy; G1 valid real MultiOmics export; G2 exact SWEET reproduction and unseen-patient handling; G3 public measured molecular phenotype task; G4 genuinely paired PK-omics cohort; G5 fit and held-out model comparison; G6 frozen independent validation. M4 may be discontinued after simpler methods if its complexity is unjustified.
