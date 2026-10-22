# Omics2PK — working to-do list (16 October 2026)

This list follows the earlier MultiOmics → SWEET → Pillar A → Pillar B R2 discussions. **When inspected the Omics2PK repository was empty; there was no pre-existing to-do file to read.** The plan below has therefore been reconstructed and must be reconciled with the owner's priorities.

## Actual supplied Phapale supplementary archives qualified (2026-10-16)

- [x] All **13 authentic Phapale 2010 Wiley legacy DOC** from owner-uploaded `paphale.zip` (1,792,400 bytes, original SHA256 b180d49fcde75d24b6d1812be001dbfa648c9565d0f5bb7100e1bc7d5c299a77), complete figures 1–8, PLS image/tables 9–12, analytical/normalization/CV methods 13. Opened and verified all original SHA hashes offline; no copyrighted files committed.
- [x] Evidence distinction: predose n29×1256 ions and postdose TAC 0–72h are truly measured in the paper; the **original numeric n29×1256 matrix, individual TAC time-concentration arrays and subject-ID crosswalk are NOT in the uploaded supplements**. Scatterplots, VIP charts, and four identified metabolites are **not** usable as participant-level covariate vectors.
- [x] Source methodology: XCMS/quantile normalization → supervised two-pass VIP screen then 7-fold/20-permutation **internal** assessment; cannot infer nested independent population validation from these figures, repeat feature screening and normalization inside train folds only.
- [x] Locally reproducible **SHA-locked audit** `scripts/audit_local_phapale_supplements.py` with synthetic-only test and private antiword parsing, original files not in CI; [full audit](PHAPALE_2010_ALL_13_ORIGINAL_SUPPLEMENTS_AUDIT_20261016.md).
- [ ] Obtain from Phapale authors original pseudonymized urine ion intensities per volunteer, measured 0–72h TAC dose/time/C(t), PK endpoint calculation and shared biological subject IDs, original batch, units and missingness flags.
- [ ] **He 2022 original 5.2MB supplemental ZIP was not included** in the supplied two archives. For true E4 renal triple-visit metabolomics and C0, original individual person×D7/M1/M3 table still required. [source request issue #31](https://github.com/rberrah/Omics2PK/issues/31).

## Newly completed original source NCA and prioritized matched pre-dose omics→PK (2026-10-16)

- [x] **REAL human mavoglurant original C(t) QC and NCA**: 120 distinct original people, 198 IV single-dose occasions, 2480 EVID0 rows (do not all count as observed DV), **53 MDV1 excluded, 2 duplicate-time occasions explicitly abstained**, **196 source-qualified NCA curves with 2403 original measured concentrations** and 76 humans with 2 qualified occasions. Median observed Cmax **556.5 ng/mL**, Tmax **0.20 h**, AUC(actual first→last) **1096.62 ng·h/mL**. There is **no real measured t0**, so AUC0-last/inf, clearance and model fitting not claimed. Original source script `scripts/validate_mavoglurant_120_real_occasion_nca.py`, [source QC report](PUBLIC_PREDOSE_TAC_NCA_AND_MAVO_LONGITUDINAL_RESULTS_20261016.md).
- [x] Audit first public **prospective Omics→PK human experimental design**: Phapale *CPT* 2010 DOI 10.1038/clpt.2009.296, 29 healthy human subjects, **1256 predose urine metabolites**, true tacrolimus plasma C(t) until 72 hours. Original **13 publisher .doc supplements** identified; **joint original subject-level PK and metabolite matrices still unverified**. Clinical PK no omics predictive proof.
- [x] Audit **He et al. 2022** DOI 10.3390/molecules27113517: n109 kidney recipients, D7/M1/M3 repeated metabolomics, tacrolimus C0 before the morning dose, n1514 retained metabolites, main outcome lymphocyte PD. Generic PMC archive URL is HTML, **not actual public ZIP**: official MDPI n5.2MB supplement and original per-subject/time C0 values still needed.
- [x] Screen Burghelea 2022 n42 renal TAC serum lipidome, selected **C0 high>8 vs low<5 ng/mL** among n135, not representative continuous exposure regression or longitudinal PK. Generic PMC 249KB ZIP URL returns HTML. Do not train unqualified clinical AUC/Cmax from case/control discrimination.
- [ ] **P0 request exact Phapale 29×1256 predose source matrix and true 0–72h TAC dose/time profiles**; compare against nested, subject-held-out clinical/dose reference, no p>>n leakage.
- [ ] **P0 request He 109 patient×D7/M1/M3 metabolite source, numeric C0/dose and exact dates** to validate true E4 changing molecular states, not just simultaneous PD association. Zhu 137/693 remains another strong E4 candidate.
- [ ] **Next real model E3**: actual 120-subject mavoglurant IV C(t) 1C/2C reference from RDA with observed infusion duration, grouped whole-subject held-out cohorts and residual/uncertainty calibration; later mechanistic PBPK and Omics2PK covariates only on true matched study participants. [issue #22](https://github.com/rberrah/Omics2PK/issues/22), [issue #31](https://github.com/rberrah/Omics2PK/issues/31).

## P0 Degraeve Wiley original supplement acquired — clinical per-patient data still absent

- [x] User supplied original **3,905,302-byte** `cpt3077-sup-0001-supinfo.zip`, exact original SHA256 verified locally. Real three-file inventory: study methods no individual table; 8 PowerPoint figures / 7 embedded Prism9 OLE source graphs; original S1/S2/S3 aggregate Word tables (15/30/5 rows).
- [x] Extract scientific qualification: original 100 participants → 3 microbiome QC exclusions + 4 eligibility exclusions = 93, stool 16S and measured TAC C0, no true patient concentration-time table or subject↔93 SRA BioSample crosswalk in official ZIP.
- [x] Different original outcomes **must not be conflated:** supplement Table S3 models daily TAC **DOSE**, R² in-sample 18.5→26.3→34.9→38.1%; article Table 2 models **estimated AUC/(dose/kg)**, R² in-sample up to **52.1%**. Neither is held-out validation.
- [x] Add original SHA-gated **local-only** methods/graph/table scanner with synthetic-only source integrity and anti-copy tests. Source archive and potential patient-level GraphPad data **never committed** to public repo.
- [ ] **The only remaining clinically material source request is the author's pseudonymized original 93-patient PK table + patient↔ENA16S↔ASV count matrix crosswalk**, with LSS measured vs imputed sampling flags. See [original actual supplement analysis](DEGRAEVE_WILEY_UPLOADED_SUPPLEMENT_AUDIT_20261016.md) / [issue #28](https://github.com/rberrah/Omics2PK/issues/28).
- [ ] Omics→observed C0/dose with age+hct+PGx reference **only once real exact patient IDs exist**; PopPK CL/F and Cmax/Tmax/AUC NCA need timed original clinical concentration series, beyond the cross-sectional n93 study. E4 dynamic covariates require a separate truly repeated molecular study.

## New verified public PK source audit: 93 TAC microbiomes and 3-analyte APAP (2026-10-16)

- [x] Original `PRJNA929904` ENA index: **93 original microbiome 16S specimens**, 93 CYP3A5 nonexpressers in Degraeve CPT 2024. Clinical C0 exists, but AUC0–24 was **limited-sampling ML-estimated with cohort-median imputation**, not dense observed NCA. Source 93 microbiome IDs do **not** create a verified subject/C0/AUC crosswalk. Both publisher original supplement URL routes returned HTTP 403; request official supplemental ZIP or original author table.
- [x] Frozen original open R `pksensi::APAP` data, SHA256 08b1d33b96d41cdd12fd48b2f09ad28c68728bf047ef4c93772a97af3dca00e0, **32 time rows / 3 literature studies / 30 complete parent APAP+glucuronide+sulfate observations**. PBPK parent/metabolite time series training/reference data at study level only, NOT 32 individual patients nor broad metabolomics/subject-wise PopPK.
- [ ] Obtain authorized original data to match `renal_patient_id ↔ PRJNA929904 sample_id ↔ dose+C0+C1+C3 ↔ measured/imputed AUC flag ↔ hct/age`. Refit patient-held-out C0/dose first, no original complete NCA/CL/F claim until richer PK source.
- [ ] PBPK baseline module: fixed three-analyte APAP parent/UGT/SULT handling, report study-level fit and uncertainty, not subject external CV; compare with actual literature and then biology only on genuinely paired individual measurements.

Full report: [public 93 TAC and APAP source audit](PUBLIC_93_KIDNEY_TAC_AND_APAP_PBPK_SOURCE_AUDIT_20261016.md).

## Original-data E3 / E4 candidates (2026-10-16)

- [x] Verify **original PLOS 2026 Wang et al TAC 36-human clinical workbook**, exact SHA, actual original C0, dose mg/kg and C/D formula, one row/patient and three aged <18. Execute a small **clinical-only reference** leave-patient-out model age+sex+BMI vs constant: MAEs **92.392803 vs 92.516919 C/D units** (virtually no uplift); not RNA/omics predictive proof, no C(t).
- [x] Verify actual **GSE15158 n88 and GSE70784 n432** transcriptomic sample arrays and longitudinal visit metadata with original GEO series matrix. No numeric source-donor-matched drug PK concentration table established; Fannin 2010 original 0–12 h APAP plasma graph is not equivalent to an original individual NCA data table.
- [x] Source rank manuscript and archives [public same-person Omics+PK triage](PUBLIC_SAME_SUBJECT_OMICS_PK_SOURCE_TRIAGE_20261016.md). Mouse Degraeve 2023 n351 source measured blood C(t)+gut ABCB1, original PK animal ID matrix missing; Zhu2022 137 recipients×693 trough+metabolomic records, exact quantified personal matrix not verified; Zenodo 24×18 biological time course no PK; HuggingFace PMxBench/PKPD-Dataset are question/text corpora, not primary clinical data.
- [ ] Verify **Wang2026 metagenome SRA PRJNA1416570** sample aliases ↔ original 36 clinical `No.` identifiers and obtain source quantified microbiome/metabolome matrices with known collection time and batch before any microbiome→C0 model.
- [ ] Examine original Fannin 2010 PDF/Supplemental PK graph and numeric dose–C(t) source; do not replace with figure-derived interpolated pseudo-measurements. Distinguish microarrays expressed relative to each patient's pre-dose baseline from comparable **absolute RNA expression**.
- [ ] Contact Degraeve2023 authors for true `mouse_id,dose,time,C(t),ABCB1,RNA_sample_id`, then test actual animal-level patient-disjoint PopPK/PBPK mechanism beyond a group mean.
- [ ] If allowed, obtain Zhu2022 original 137 patient×693 repeated omics+C0/D features/dates/PGx. Until then no longitudinal tacrolimus omics PK claims.
- [ ] Only after public original full match: fit appropriately sampled exposure NCA first, then true PopPK C(t), PBPK physiological mechanisms and time-varying E4 states.

## P0 interoperability — consume real MultiOmics outputs (E2 → E3 → E4)

- [x] Audit **current pinned actual MultiOmics browser exchange V1** from source and full research ZIP tests. Verify raw matrices / assay-subject-feature identity, source-scale declarations and research analysis provenance.
- [x] Prove **loss of T0/T12 categorical timepoints** on current producer: input 8 biological subjects × 2 specimen visits, 49 omics assays, but 49/49 empty exported collection times. A technical RNA duplicate is one repeated assay, **not a ninth subject or third visit**.
- [x] Local read-only source-verified Omics2PK readiness report and tests: numeric time/unit, repeated specimens, incomplete provenance, no real dose/C(t)/individual PK endpoints; **fails closed** without fabricating temporal alignment.
- [x] Add strict source-provenanced optional `individual_signatures` schema for later RNA/CYP/pathway/MOFA/SWEET features per actual assay. Group logFC/q values or transductive networks cannot masquerade as held-out individual covariates.
- [x] Add strictly measured `Cmax/Tmax/AUC0last` NCA endpoint input QC with actual source hash, physical units, original dose and time window, independent subject/occasion IDs; additional study-clock verification required to link to preceding omic specimen.
- [ ] **MultiOmics producer owner** to preserve original `timepoint` and `technical_replicate` (not inferred units); real producer commit / molecular provenance; optionally export fixed individual scores. Explicit cross-repo [MultiOmics issue #20](https://github.com/rberrah/MultiOmics/issues/20), not a direct source edit by Omics2PK.
- [ ] After producer fix, rerun actual browser ZIP end-to-end and confirm original RNA-only and multiomics patient/sample/time relationship survives; then evaluate on **genuinely matched longitudinal clinical PK + omics**, not synthetic multiomics combined with unrelated Theoph patients.
- [ ] NCA signal stage: separately test Cmax/Tmax/AUC0-last as dose/time-window-specific outcomes with actual PK and baseline molecular states; account for sampling time differences and multiplicity. **NCA endpoints are NOT intrinsic PopPK structural parameters.**
- [ ] E3 model stage: fit original dose-C(t) through Pillar B/approved MIPD-Engine, clinical+PGx reference then frozen biological covariates on unchanged subject-disjoint external folds; compare actual PK prediction, CL/F, V/F, absorption, PK uncertainty. PBPK requires mechanistic transfer and physiology, not transcript→CL direct equivalence.
- [ ] E4: dynamically evolving molecular state prior to each prospective prediction landmark, patient-disjoint + temporal anti-leakage, transportability assessment, Pillar C abstention. NewPK remains downstream comparator only.

Full source-level audit: [REAL_MULTIOMICS_OUTPUT_AUDIT_AND_PK_ADAPTER_GAPS_20261016.md](REAL_MULTIOMICS_OUTPUT_AUDIT_AND_PK_ADAPTER_GAPS_20261016.md).

## P0 NEW priority — actual longitudinal PK, then time-varying omics

- [x] Implement **PK_LONG_V1** for actual dose and C(t) events: dosing-history completeness, BLQ/LOQ, t0 predose, observed vs simulated, matrix/units and temporal order. **A nearest-dose-only join is forbidden for multi-dose PK.**
- [x] Implement molecular state sampled before each prediction landmark with fixed lag and future-omics rejection. Post-dose RNA can support later prediction only, not retrospective clearance prediction.
- [x] Reproduce source-pinned **R Theoph 12 true patients × 11 observed concentrations (132)** with correct oral mg/kg→mg doses, true sample timing, source SHA-256 and NCA0-last descriptive exposure (no CL/F or AUCinf).
- [x] Source-audit original **theo_sd/theo_md/warfarin/mavoglurant** RDA using pyreadr in ephemeral CI, with observed vs simulated, PK vs PD, EVID and OCC distinctions. Freeze warfarin 515 file rows vs 519 reference-text discrepancy.
- [x] Discover **mavoglurant real 120 people/198 single-dose IV occasions**, 78 subjects with two occasions; 2,480 recorded concentrations, no continuous multi-dose within occasion. Source ng/mL confirmed by nlmixr2 vignette; calendar gap between occasions missing.
- [ ] Add source-level PK event adapter for actual **mavoglurant 198 occasion profiles** (time reset per OCC; do not count 198 biological subjects), then model with qualified appropriate baseline in Pillar B.
- [ ] Convert **warfarin PK and PCA-PD** into distinct endpoint streams; confirm original 515 vs documented 519 entries.
- [ ] Research original matched longitudinal drug plasma PK + repeated omics: Fannin 2010 acetaminophen 6 treated/3 controls with 0–12h drug C(t) and delayed serial RNA/metabolome, actual per-subject matrices not yet available; RNA measured 48h cannot forecast earlier drug Cmax/AUC.
- [ ] Next: phase dose/C(t) data into B structural/SAEM model comparisons, respecting actual subject split and same biological patient over time; evaluate clinical+PGx then prespecified omic-state improvement on truly matched records.

**Source and scientific limitations:** [actual longitudinal PK audit](ORIGINAL_LONGITUDINAL_PK_SOURCE_AUDIT_20261016.md) · [longitudinal biological state roadmap](LONGITUDINAL_PK_AND_DYNAMIC_OMICS_ROADMAP_20261016.md).

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
- [x] Reproduce official SWEET sample weights (80/80), raw edges and global summary. **Archived z-score mismatch root cause resolved:** historical upstream code used `raw - mean/std`; current code uses the correct `(raw-mean)/std`. Both archived 4,753-edge patient files match the legacy precedence bug exactly and are rejected as current numeric references.
- [x] Qualify the SWEET `-z` behavior and archived-output provenance; current parenthesized formula is the only allowed reference in Omics2PK.
- [ ] Audit remaining SWEET handling of zero variance, zeros and missing values independently of the resolved archive bug.
- [x] Prototype a clearly labeled **non-equivalent SWEET-inspired** training-frozen reference and query-only network (PR #3).
- [x] Add subject-disjoint training leave-one-out/cross-fitting and outer train/test feature construction, with negative leakage tests (PR #3, merged).
- [x] Separate genome-wide sample weights (100 official example genes) from a targeted edge panel (25 fixed-order genes), and compare against the official transductive equation. Public diagnostic: 60 train/20 held-out, one clipped weight, ~27.1% relative Frobenius difference for one test network. No outcomes or clinical validation.
- [x] Validate a genuine out-of-sample **query-anchored v1** construction mathematically against published transductive SWEET while keeping the reference frozen. Public 60/20 diagnostic uses all 100 official genes for sample weighting and a fixed 25-gene edge panel: v1 is closer than v0 for 19/20 held-out subjects (mean relative Frobenius difference 0.26729 vs 0.27257), but remains explicitly non-equivalent because the official reference is transductive.
- [ ] Compare v1 stability to preprocessing, batch, tissue and threshold selection. **Reference-size/composition sensitivity is now qualified outcome-free on real GSE9588.** Broad 512-probe weighting and order-independent SHA-256 references remove most earlier shift. Nested reference-size sensitivity gave clamp fractions **4.51% (n50), 3.13% (n75), 3.67% (n100), 0% (n150)**. Four additional predeclared salted SHA-256 n150 references show **0–1.805%** clamp fractions: 2/5 total n150 references have zero clamps, the others have 1, 3 or 5 clamps. Therefore n150 is materially better calibrated but **not** a universal zero-clamp threshold. **Train-only preprocessing infrastructure is now implemented and CI-qualified** for source-normalized passthrough, train-only imputation/z-scaling and per-sample CPM/log1p raw-count handling; held-out perturbations cannot change the preprocessing fingerprint or training features. Remaining gate: source-specific preprocessing/batch/tissue sensitivity on authentic cohorts; generic batch correction remains intentionally unimplemented.

## Thesis Chapter 4/10 -> genuine PK proof: cross-domain gate
- [x] Pin the **15 October 2026 frozen thesis architecture** as a public research mapping without exposing potentially patent-sensitive Pillars A/B/C/D source (docs/OMICS_TO_PK_THESIS_PROOF_AND_SITE_PLAN.md).
- [x] Implement strict baseline omic-to-observed-PK patient/date/dose/concentration metadata check: RNA-only or 2–3 omics; time zones, prior-to-dose sampling and fixed assay alignment gap; private data never uploaded.
- [x] Audit source file structure for Wegler 2021 PDF and donor-matched proteome workbooks using exact official public SHA-256; source PDF is native text and does not require OCR.
- [x] **Wegler donor-level crosswalk qualified from original Table S3**: strict parser recovers exactly 15 donors × 5 drugs, same donor set across drug blocks, HLM/HH CYP protein concentrations, Kpuu and source clearance equations. `eq.1` HLM CLint,mic and `eq.4` HH CLint,hep are independently validated against the article's rounded midazolam aggregates; 2 bufuralol HLM ND rows preserved.
- [x] First measured in-vitro donor-held-out CLint regression with fixed corresponding CYP protein: per-drug LODO, fixed Ridge alpha=1, no tuning. **HH:** protein improves descriptive held-donor log-MAE materially for midazolam (Δ −0.272), omeprazole (−0.472) and diclofenac (−0.091), with paired bootstrap intervals below zero; bupropion ~neutral, bufuralol worse. **HLM:** pooled change negligible (Δ −0.008) and every per-drug bootstrap interval crosses zero. This is real E2 in-vitro evidence, not clinical PK. See `docs/WEGLER_SOURCE_STRUCTURE_RESULTS_20261016.md`.
- [ ] Separate **transcriptomics→measured PK** source (or functional assay), not just phenotype→tumor response.
- [ ] Actual clinical BIOMARGIN/authorized linked tacrolimus or MPA TDM cohort with patient-date/dose/concentration links (ethics/license gate).
- [ ] Compare clinical+PGx reference → prespecified direct RNA scores → qualified SWEET RNA → parallel protein/metabolite → optional integrated latent scores, on identical held-out patients; external validation + Pillar C adjudication.
- [ ] **After proof only**: design public learning/tool at rberrah.github.io with FR/EN, dark/light, Resources, Ecosystem, persistent provenance and no user clinical-data upload.

## E2 REAL original measured biological phenotype→drug-metabolism proof (2026-10-16)
- [x] Tucatinib: n=21 independent original HLM donors, published source Table 1 manually transcribed and **all four published Spearman CYP↔CL correlations reproduced**; 21 leave-donor-out folds MAE mean 6.246 vs 2 CYP proteins 2.041 vs 2 CYP enzyme activities 1.940 µL/min/mg, explicitly nonclinical and original single cohort. Corrected source narrative/table inconsistency identified; donor measurements not redistributed. Detailed [tucatinib proof](TUCATINIB_21_DONOR_MEASURED_VALIDATION_20261016.md).
- [x] **RNA/PK bridge now has genuine measured data:** Liu et al US EPA original 100 human livers, exact match across two original sample IDs, n=94–100 complete paired donors in 8 transcript→microsomal CYP activity assays, genuine donor-held-out train/test runs on the original source in GitHub CI. Compared age+sex versus age+sex+RNA, 7/8 lower MAE, 4/8 uncorrected paired bootstrap intervals exclude 0. **CYP2E1 worsened**. Source/licensing/data provenance in [EPA RNA-CYP proof](EPA_100_DONOR_TRANSCRIPT_CYP_VALIDATION_20261016.md).
- [x] Independent numerical replication: pinned scikit-learn StandardScaler/Ridge used on actual original n100 liver source in a distinct script; all **8 CYP held-donor-out MAEs match the native solver to six decimals**, GitHub Actions https://github.com/rberrah/Omics2PK/actions/runs/38930709329. This is independent *estimator* parity, NOT independent donor/site external validation.
- [ ] R/Julia third-party replication beyond sklearn, and multiple-testing-aware inference, sensitivity to CYP outliers, donor age/batch, and original paper processing consistency.
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
- [x] Wegler et al.: original Table S3 now provides a source-verified same-donor targeted CYP-protein ↔ measured CLint join for 15 donors × 5 probe drugs; donor-LODO executed with heterogeneous HH/HLM results. **No SWEET/RNA comparison is possible from Wegler because matched RNA is absent.** Global Data S1 proteome expansion remains optional and should not be post-hoc fitted at p≫n.
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
