# Omics2PK

## Real public hepatic multiomic source qualification and Omics2PK benchmark reruns (18 October 2026)

An automated, **executed** DOI/PMCID/original-supplement audit of six human liver RNA/protein/metabolic-function papers and GEO GSE9588 now documents actual source accessibility (not just published cohort sizes). No new public **same-donor RNA+protein+measured metabolic substrate clearance** numerical matrix could be source-qualified from ASPET/ACS/Wiley official endpoints; original Synapse HLC **curatedPhenotype `syn4629`** and curatedExpressionLiver `syn88644` entities do exist, but anonymous file-content API returned **HTTP 403** for both (and for the genotype/parent files); **do not claim RNA vs RNA+protein superiority**. GEO GSE9588 **original 39,302 microarray probes × 427 GSM sample columns** was downloaded and verified numerically; an actual 75-reference + 352-query frozen SWEET-inspired demonstration ran but produced **238 distinct graphs for 352 projections**, exposing weight clamping/nonindividualized-signature risks. Original genuine EPA n100 targeted CYP RNA→activity, Latham n21 hepatic CYP proteins→tucatinib/imatinib metabolism, and Rotroff n32 predose metformin→clinical drug AUC/CLr benchmarks were rerun in CI and reproduce prior findings. [Detailed source-by-source actual audits, GEO QC, real model results and scientific gates](docs/HEPATIC_PUBLIC_DATA_EXECUTED_SOURCE_AUDIT_20261018.md).

## Identified public human-liver transcriptome/protein/CYP metabolism source targets (18 October 2026)

Independent original source discovery prioritized published *same-donor* hepatic RNA/protein/observed substrate metabolism: **Zhang n100 ten CYPs**, **Tanner n360 CYP2A6**, **UGT2B17 n455 with RNA n230**, **Ohtsuki n17**, Jamwal n74, and HLC/Sage **GSE9588 n427 RNA** with separately published n398 CYP activities. **No complete original numerical three-layer source dataset has yet been retrieved and verified across same donor IDs; no RNA-vs-RNA+protein advantage may be claimed**. All original-data gates, source DOI and dataset links, model endpoints, subcohort denominators and retrieval priorities: [hepatic public source catalogue](docs/HEPATIC_PUBLIC_RNA_PROTEIN_DRUG_METABOLISM_SOURCE_CATALOGUE_20261018.md).

## Targeted CYP3A4 RNA versus protein for hepatic midazolam metabolism — source gate (18 Oct 2026)

A same-donor comparator for **74 human livers** now exists for CYP3A4 mRNA alone, CYP3A4 protein alone, clinical context, and combinations against **measured in-vitro midazolam 1′-hydroxylation**. Synthetic contract tests validate donor-heldout fitting, source checks and output privacy. **Actual 74 original RNA/protein/activity values were not recovered from the official source endpoints; therefore NO real RNA-vs-multi-omic numerical model was fitted or claimed.** See [source audit and original-only experiment protocol](docs/JAMWAL74_CYP3A4_RNA_PROTEIN_MIDAZOLAM_ORIGINAL_SOURCE_AUDIT_20261018.md). This comparison uses targeted qPCR/protein, not comprehensive transcriptomics/proteomics, and predicted paper PBPK is never observed midazolam PK.

## Hepatic metabolism on two real drugs: protein vs activity and fusion (18 October 2026)

Two truly original source-matched studies of 21 human liver microsome donors measured CYP2C8/CYP3A protein abundance, enzyme-activity assays and **tucatinib metabolic clearance**; 17 of the same donors also had quantifiable **imatinib metabolic clearance**. With identical donor-heldout validation, adding functional probe activity to CYP proteomics reduced tucatinib MAE **2.041→1.783 µL/min/mg**, and imatinib MAE **14.166→13.620**. Small n, paired intervals include no gain for clearance, no independent-donor external validation. **This is protein + functional phenotyping, NOT RNA + protein multiomics**: same-donor transcriptome is unavailable and no RNA-only comparison is asserted. Real PMC article source downloaded transiently, verified published imatinib CYP rank correlations, code and outputs in [qualified hepatic metabolism report](docs/HEPATIC_TWO_DRUG_CYP_PROTEIN_ACTIVITY_FUSION_REAL_20261018.md), original scientific [GitHub Actions result](https://github.com/rberrah/Omics2PK/actions/runs/39147377727).

## First genuine public clinical E3 biomarker-to-measured-PK benchmark (18 October 2026)

**Rotroff et al. metformin original n33 (DOI 10.3389/fphar.2016.00135) is SOURCE-QUALIFIED and actually modeled.** The official original publisher supplements contain 33 exact-person identifier matches between 33 predose A/B/C metabolic GC-TOF triplets, **384 predose measured molecular features**, and true original person-level metformin AUClast0–24 and CLr renal clearance (**32 complete individuals per PK endpoint**). Public source workbooks are SHA256 pinned; no subject data copied into GitHub. Actual true source numerical run: [E3 original measured metformin PK CI](https://github.com/rberrah/Omics2PK/actions/runs/39145479441).

32-person heldout-fold Ridge benchmark: age+sex+weight vs same clinical baseline plus three **train-only** unsupervised metabolomics PCs. **AUC MAE 2826.05 → 2825.02 (essentially zero benefit; uncertainty includes no gain), CLr MAE 115.11 → 118.14 (worse)**, errors in original spreadsheet numerical units pending physical unit verification. These are authentic clinical PK measurements, but the tested omic augmentation **is not superior**. No independent external validation or SWEET/XGB model is claimed. [Full documented source join, clinical E3 model, results and limits](docs/E3_REAL_METFORMIN33_PREDOSE_384_OMICS_AUC_CLR_HELDOUT_20261018.md).

## Real human-liver molecular-to-metabolism proof (E2, 18 October 2026)

On the ORIGINAL EPA n100 human liver cohort (source DOI 10.23719/1519197), eight donor-matched qPCR CYP transcript → measured microsomal enzyme-activity endpoints were independently validated using original biological donor heldout folds. Adding genuine same-donor RNA beyond age/sex reduced MAE in **seven of eight**; after **199 same-sex donor RNA mismatching permutations per enzyme**, **four of eight** original CYP assays (CYP1A2, CYP2B6, CYP2D6, CYP3A4) showed a predeclared one-sided empirical permutation p=0.005, below Bonferroni p=0.00625. The CYP2E1 negative result (worse MAE) is retained. **This validates E2 functional metabolism association, never clinical individual plasma PK/PopPK.** [Source-verified original-data results and reproducible negative-control methodology](docs/EPA100_REAL_199_PERMUTATION_NEGATIVE_CONTROLS_20261018.md).

## E3 public clinical omics–PK source availability (18 October 2026)

The five processed Silveira/Duarte rosuvastatin n40 Google Drive CSVs are **private**, confirmed by the user. Do not pursue new public-link download attempts or claim published AUC/Cmax or person-level omics matching from original EBI mzML alone. The original Wegler rosuvastatin **n54** published supplementary XLSX was actually obtained via official Europe PMC and verified to contain **54 donor protein abundance/allometric records, but no individual AUC/Cmax columns**. Both original clinical datasets are therefore **E3 source-blocked**, despite genuine original raw mass spectra or protein matrices. Focus next on verified independent E2 human liver omics→activity/CLint datasets; no E2 source may be promoted to clinical PopPK evidence. [Detailed source verdict](docs/E3_PUBLIC_DATA_PIVOT_SILVEIRA_CLOSED_TO_WEGLER54_20261018.md).

## E3 rosuvastatin: original author GitHub and public EBI compared (18 October 2026)

The author's complete GitHub tree has the 21MB original notebook, README and requirements, **not the five processed CSVs**. Those five dataset URLs are hosted externally on Google Drive in the notebook. EBI MTBLS3128 stores source LC-MS files, but those are not equivalent to the ProgenesisQI/MetaboAnalyst feature matrices. A complete original-cell publisher XLSX audit found **no recognized 40-person individual AUC/Cmax table**; a source-verified EBI ZIP inventory then examined **all 156/156 archive central directories** (3,174 internal files) and found **zero named CSV/XLSX/tabular-format members**, without reading full raw contents. The original numeric predose molecular signals and PK same-person crosswalk remain unavailable: **no measured clinical E3 effect is claimed**. [Actual cross-source comparison and limitations](docs/ROSUVASTATIN40_GITHUB_VS_EBI_FIVE_PROCESSED_CSV_VERDICT_20261018.md).

## Public EBI MetaboLights MTBLS3128 — local data download (18 October 2026)

The original public **ISA metadata and raw LC-MS mzML** can be downloaded via official EBI anonymous FTP using \`python scripts/download_mtbls3128_original_public.py --list-only\`, then \`python scripts/download_mtbls3128_original_public.py\` for original ISA metadata only; optional bounded \`--mzml-max-files N\` downloads raw mass spectra. Files and SHA256 manifests stay locally under Git-ignored \`private/\`; **none is bundled in GitHub or Actions**. [Instructions and limitations](docs/MTBLS3128_OFFICIAL_PUBLIC_DATA_LOCAL_DOWNLOAD_20261018.md). This still does **not** provide the 40 matched original numeric rosuvastatin AUC/Cmax targets or original processed metabolite CSVs.

## E3 rosuvastatin: private original-source held-human-out evaluator (18 October 2026)

A **new local-only exploratory n40 reference** now evaluates author-verified rosuvastatin **AUC/Cmax** against a fixed, subject-disjoint clinical/constant versus clinical+predose-metabolome ridge comparator, with train-only imputation/normalization and paired uncertainty. It requires the original **one-person-per-row source**, pinned SHA256, separately documented source authenticity and four explicit manual attestations. **The authentic n40 quantitative PK ↔ predose metabolome files are not yet obtained/qualified; no human E3 results are claimed.** [Private execution protocol, source gates and limitations](docs/ROSUVASTATIN40_PRIVATE_E3_HELDOUT_PROTOCOL_20261018.md). No patient-level rows or predictions are uploaded to CI.

## New Omics2PK research website and original SMX/n26 benchmark (17 October 2026)

A **standalone editorial research homepage** is prepared under [site/](site/) with the same PMx Explain / MultiOmics visual conventions, real dropdowns **Learn / Analyze / Other Disciplines / Ressources**, FR/EN and dark/light; static preview passes CI. [Site and first GitHub Pages deployment guide](docs/OMICS2PK_WEBSITE_BRAND_PARITY_AND_PUBLISH_20261017.md). **Do not claim the intended URL is live until Pages has actually been enabled and verified.** Public site has *no hosted model, patient uploads or dosing tool*.

A different original **human NAT1/NAT2** cohort (Lanphier et al. DOI 10.1016/j.dmd.2025.100195) has now been **privately numerically reanalyzed** from published Table 1 with **26 distinct PHH donors**, n21 measured NAT activities/proteomics, n20 NAT2 genotypes and strict 15-way common donors. On **identical 15 donors** the baseline / activity / protein / genotype LODO MAEs are respectively **0.041914 / 0.037981 / 0.043823 / 0.039271 pmol/min/million cells**; all exploratory paired-bootstrap uplift intervals include zero. Four article Spearman coefficients are independently reproduced; prior peptide ↔ r label mapping was corrected. [Measured source audit, both estimators and limitations](docs/SMX_NAT_26_ORIGINAL_MEASURED_LODO_RESULTS_20261017.md). **No original individual donor rows or publisher Table 1 reproduced.** Live machine-downloadable original JATS still inaccessible; published-Table-1 private curation is not a public CI source download.

A new **clinical E3 lead** is [Tee et al. 2022 metformin n17](docs/METFORMIN17_CLINICAL_OMICS_PK_OPEN_SOURCE_AUDIT_20261017.md) with 255 original plasma concentration observations and pre-dose/peak LC-MS metabolomics **in the same clinical study**. Published n17 PK Table S5 and actual personal m/z intensity↔PK matrix linkage **have not both been verified accessible**, so this is a *candidate*, not an Omics→clinical-PK prediction result.

## New confirmed source-level sensitivity: eight related hepatic CYPs (17 October 2026)

On the **same 100 genuine source-verified U.S. EPA human livers**, eight held-donor-out RNA→observed microsomal CYP activity models were re-evaluated using a **shared-donor paired 5,000-resample bootstrap** with **Bonferroni-adjusted 99.375% percentile intervals**. Seven of eight have favorable RNA point MAE; four have nominal unadjusted 95% intervals excluding zero, but only **three survive the eight-endpoint adjustment**: **CYP1A2, CYP2B6 and CYP3A4**. **CYP2D6 does NOT survive** multiplicity correction, and **CYP2E1 is directionally worse**. Actual original source checksum and exact original estimator MAE parity checked before analysis; no donor rows committed. [Full numerical source evidence and caveats](docs/EPA_100_LIVER_EIGHT_CYP_BONFERRONI_REAL_RESULTS_20261017.md). This establishes exploratory **RNA→in-vitro enzyme activity**, not omics→clinical PK.

## Latest source qualification: same-donor two-drug metabolism and adaptive NLME (16 October 2026)

Two independently published in-vitro drug experiments — **tucatinib and imatinib** — reuse the **same 21 individual human-liver-microsome donor lots and CYP phenotypes**. Their original Tables 1 make a direct same-person, different-substrate study possible; **only 17 imatinib apparent clearance values are quantifiable**, while all 21 provide N-desmethyl-imatinib formation peak-area ratios. The original published tables have now been manually source-verified and independently cross-checked against the eight original published Spearman correlations. **Observed donor-held-out MAE for tucatinib n21:** baseline 6.245714 vs two-CYP protein 2.041447; **imatinib original quantified n17:** 28.869894 vs 14.166086 µL/min/mg. N-desmethyl-imatinib peak-area-ratio n21 improves from 0.146683 to 0.082222 (different units). Four imatinib no-substrate-depletion lots remain unquantifiable and are never replaced by zero. No original donor rows were rehosted. [Full original measured reanalysis and uncertainty limits](docs/MEASURED_TKI_TWO_DRUG_ORIGINAL_RESULTS_20261017.md), [source identity audit](docs/TWO_DRUG_SHARED_DONOR_PUBLIC_SOURCES_20261016.md).

A separate public **sulfamethoxazole/NAT1–NAT2** primary hepatocyte study now adds 26 genuine original donors with observed NA-SMX formation, 21 NAT activity probe donors and 20 NAT2-genotyped donors. Its genotype/protein/activity versus clinical/demographic local-only benchmark is implemented and source-gated, but **no real-source SMX predictive performance is claimed yet**. [Original source scope and study plan](docs/PUBLIC_SMX_NAT_26_DONOR_QUALIFICATION_20261017.md).

On the separate original **mavoglurant** clinical concentration dataset, adaptive Gauss-Hermite integration of a narrow 2C IV IIV(CL) NLME pilot stabilized the held-out joint predictive NLL at **3.908502 per biological person for 7/11/15/21 nodes**, unlike prior unadapted GH; fixed 2C RMSLE **0.383383** and adaptive NLME **0.387256**, hence no evidence of improved unconditioned population point predictions. The source-hash-verified fivefold benchmark and synthetic tests [passed in CI](https://github.com/rberrah/Omics2PK/actions/runs/38973121688). [Numerical results and limitations](docs/MAVOGLURANT_ADAPTIVE_GH_NLME_ORIGINAL_RESULTS_20261016.md). **This is not independent clinical Omics→PK evidence or a qualified SAEM/FOCEi model.**

## Mavoglurant IV held-out structural PK comparator (preparatory)

A pre-registered **1C vs 2C finite-IV-infusion pooled fixed-effect reference** now operates on the **original SHA-verified 120-human, 196-qualified-occasion, 2,403-measurement mavoglurant source**. Five human-disjoint folds keep repeat occasions together, with observed RATE/AMT and ng/mL units protected, and do not permit subject-specific post-hoc fitting. The reproducible workflow reports only held-out aggregate log prediction errors. **This is not SAEM/NLME PopPK, not an Omics→PK effect and not external omics validation**. See [protocol and interpretation](docs/MAVOGLURANT_ORIGINAL_IV_HELDOUT_STRUCTURAL_BASELINE.md). Its fivefold source-verified numerical validation has completed successfully; see the linked original-source results and the newer NLME qualification above.

**Research-only alpha. Not clinically validated.**

Omics2PK connects the reproducible outputs of [MultiOmics](https://github.com/rberrah/MultiOmics) to research on patient-specific molecular networks (SWEET), drug-specific knowledge from Thesis Pillar A and PopPK/PKPD covariate discovery in Pillar B (route R2).

## New: transition to longitudinal PK (October 2026)

The immediate goal is now **real original patient concentration–time profiles and dosing histories**, followed by **time-varying RNA/protein/metabolite covariates** on matched subjects, rather than additional cross-sectional enzyme correlations. The source-gated PK_LONG_V1 contract verifies observed versus simulated events, predose sampling, BLQ/LOQ, full dose history and prospective molecular sampling timestamps. Original **Theoph**: 12 humans, 12 doses, 132 actual C(t) observations, Cmax/Tmax/AUC0-last descriptive NCA. Separate pinned original RDA audit covers 120 humans, 198 distinct PK occasions and 2,480 measured concentration rows for **mavoglurant**, 32 subjects with warfarin PK/PD mixed records, and explicitly simulated theophylline multidose data. **None of these original PK series contains omics**; they are safe PK architecture/clinical-series qualifications, not proof of transcriptome-driven precision dosing.

See [Longitudinal Omics→PK roadmap](docs/LONGITUDINAL_PK_AND_DYNAMIC_OMICS_ROADMAP_20261016.md) and [original measured PK dataset audit](docs/ORIGINAL_LONGITUDINAL_PK_SOURCE_AUDIT_20261016.md).

## MultiOmics source-interoperability audit: E2 molecular evidence → E3 PK/NCA

We have reproduced the **actual** `rberrah/MultiOmics` V1 browser exporter against its current pinned public source. The native 8-subject, 16-specimen, 49-assay demo **loses 100% of the original categorical T0/T12 visit metadata** on downstream export, because the producer copies `collection_time` but not the input `timepoint`. Raw RNA, unknown provenance scales and missing PK dosing histories remain correctly blocked for clinical interpretation. Omics2PK now offers `python -m scripts.audit_multiomics_pk_readiness <local_exchange_dir>` for locally verifying the original V1 output, plus a strict optional `individual_signatures` contract and measured `Cmax/Tmax/AUC0last` NCA target contract. No attempt is made to derive patient PK from group fold changes, an omics-only ZIP, or unverified chronological labels.

[Full native output audit and producer change request](docs/REAL_MULTIOMICS_OUTPUT_AUDIT_AND_PK_ADAPTER_GAPS_20261016.md). The agreed cross-repo changes are filed as [MultiOmics issue #20](https://github.com/rberrah/MultiOmics/issues/20); this repository does not alter MultiOmics source. Omics→NCA→PopPK/PBPK evidence remains future until genuinely observed, subject-linked data and time-validated covariates exist.

## Public patient-linked omics and PK original-source hunt (2026-10-16)

A first **real individual clinical tacrolimus exposure target** is now available from Wang et al. 2026 PLOS ONE [clinical S1 original supplement](https://doi.org/10.1371/journal.pone.0343817): **36 liver transplant recipients with dose, actual whole-blood C0, C/D, age, sex, height and weight**, including **three aged <18**. The source spreadsheet SHA-256 and original C/D physical dose/weight consistency are checked in CI. A **clinical-only** 36-fold held-patient-out reference (age+sex+BMI fixed ridge vs constant original-source baseline) showed nearly identical MAE; **no omic measurements were joined or fitted**. NCBI SRA metagenome accession PRJNA1416570 is public, but an exact source-defined crosswalk to this clinical table and the original metabolite matrix **remain unverified**.

Original **Fannin 2010 GEO GSE15158 (88 RNA assays)** and **Bushel/Fannin 2016 GSE70784 (432 RNA assays)** were live audited; both have authentic serial blood gene expression but neither currently supplies source-verified individual clinical PK concentration-time tables. Fannin's published early serum APAP graph is not a primary numeric curve; its post-dose RNA is a downstream response, not an admissible baseline predictor of early AUC. Other source candidates include Degraeve et al. 2023 n351 genuine mouse TAC PK blood samples + gut ABCB1 (unverified animal RNA/PK crosswalk), Zhu et al. 2022 n137/693 clinical TAC C0+metabolomics (original patient matrix not openly verified), Zenodo 24-participant 18-visit metabolomics (no original drug PK), plus source-qualified negative HuggingFace training corpora.

**Full audit and strict evidence decisions:** [Public Omics×PK source triage](docs/PUBLIC_SAME_SUBJECT_OMICS_PK_SOURCE_TRIAGE_20261016.md). The repository will not stitch unrelated omics and PK cohorts together or treat a C0/D endpoint as PopPK CL/F/AUC.

## Public PK/omics evidence gate: 93 human TAC microbiomes; parent-metabolite PBPK

In the original [Degraeve et al. 2024 tacrolimus renal cohort](https://doi.org/10.1002/cpt.3077), the genuine ENA `PRJNA929904` source contains **93 distinct original 16S samples** for 93 clinical CYP3A5-nonexpressing transplant recipients. Actual **C0 was observed**, while the paper's AUC0-24 was **machine-learning estimated** from sparse concentrations (with missing-point cohort-median substitution). Original individual C0/AUC ↔ bacterial specimen crosswalk has NOT been obtained, and automated official supplementary ZIP requests returned HTTP 403. Thus these are **not** ready for a patient-validated Omics→AUC study.

Separately, the public [R pksensi APAP dataset](https://search.r-project.org/CRAN/refmans/pksensi/html/APAP.html) contains **32 source rows from three literature studies, including 30 complete parent-acetaminophen / glucuronide / sulfate observations**, with real drug metabolite concentration-time signals. They can support **PBPK multi-analyte population-curve checks**, but there are NO original patient identifiers or independent RNA/metabolomic feature matrices in this RDA.

See [strict original-source audit and E3/PBPK limitations](docs/PUBLIC_93_KIDNEY_TAC_AND_APAP_PBPK_SOURCE_AUDIT_20261016.md). No independent E3 clinical Omics→PopPK gain is claimed.

## Degraeve 93 renal TAC original source supplement: acquired, raw PK/16S crosswalk still missing

The **exact authorized original Wiley** `cpt3077-sup-0001-supinfo.zip` was inspected locally (SHA256 `03aab5deb292a03753073e493d7ae49896884be2cde22d5c2d1c9f3037488cc2`; **not redistributed in public GitHub**). The three original files provide clinical/16S methods, 8 figures with GraphPad source objects, and 3 aggregate tables. There is **no documented original 93-row patient clinical C0/C1/C3 ↔ SRA/ASV matrix crosswalk**. This resolves the former publisher-403 archive access blocker, but NOT the patient-level biomarker→PK data blocker.

Important independent endpoints: supplementary **Table S3 is TAC daily *DOSE* (not AUC), n93 in-sample R² 38.1%** with age, hematocrit, ASV1508, ASV664; the main paper's **Table 2 is sparse ML-estimated AUC/(dose/kg), n93 in-sample R² 52.1%**. Neither establishes held-out clinical Omics→PK predictive gain. See [actual local source audit](docs/DEGRAEVE_WILEY_UPLOADED_SUPPLEMENT_AUDIT_20261016.md), [crosswalk issue #28](https://github.com/rberrah/Omics2PK/issues/28). Local-only archive inspector: `python -m scripts.audit_local_degraeve_supplement /path/to/original.zip` (CI uses synthetic documents only).

## Measured longitudinal human PK source: 196 actual NCA visits (October 2026)

The first original **measured patient-level multi-occasion exposure descriptors** have now been numerically evaluated on source-pinned `nlmixr2data::mavoglurant`: 120 independent humans, 198 distinct IV single-dose occasions, 196 eligible NCA timecourses (2 duplicate-time source occasions abstained), **2,403 MDV0 positive source concentrations** after excluding **53 MDV1 records**, with **76 people contributing two eligible visits**. The 196 published study-occasion series yield median Cmax **556.5 ng/mL**, observed Tmax **0.20 h**, and area between actual observed first and last sample **1,096.62 ng·h/mL**. There is **no original quantified t=0 concentration**, so no invented AUC0 or AUCinf or CL. Model training and external held-subject prediction are still pending. Reproducible original analysis in `scripts/validate_mavoglurant_120_real_occasion_nca.py`, CI and [scientific audit](docs/PUBLIC_PREDOSE_TAC_NCA_AND_MAVO_LONGITUDINAL_RESULTS_20261016.md).

New paired biological/PK source high priority: **Phapale et al. 2010 29 healthy volunteers**, tacrolimus up to 72h **after** 1,256 urine metabolite signals collected **before** treatment; same-person original matrix and exact numeric time profiles await author/Wiley supplement inspection. **He et al. 2022 109 renal recipients** have triple D7/M1/M3 blood metabolomics and trough TDM predose, but per-visit numeric source/ID export is not yet verified; primary PD target is T cell activity, not PK. `He 2022` original supplement URL currently responds with HTML to a generic PMC attachment request, not raw ZIP. Do not overstate these studies' numerical E3/E4 reproducibility. See [public source audit](docs/PUBLIC_PREDOSE_TAC_NCA_AND_MAVO_LONGITUDINAL_RESULTS_20261016.md) and [issue #31](https://github.com/rberrah/Omics2PK/issues/31).

## User-provided complete 2010 tacrolimus pharmacometabolomics supplements (audited October 2026)

The user supplied **all thirteen genuine Phapale et al. 2010 CPT supporting .doc files**, SHA256 of the original ZIP `b180d49fcde75d24b6d1812be001dbfa648c9565d0f5bb7100e1bc7d5c299a77`. All original scientific figure/table/method files were privately reviewed. They reproduce PLS/LC-MS methodological evidence, including quantile normalization, VIP selection, internal 7-fold validation and 20 permutations, and four named candidate metabolites. **However, the 29-person×1,256 predose urine data matrix, individual 0–72 h tacrolimus C(t) numerical profiles, and same-subject omic/PK identifiers are NOT provided**, so Omics→clinical NCA/PopPK performance remains **unverified**. A local SHA-gated source audit script and synthetic-only CI are included; **the copyrighted original ZIP and any clinical observations are never committed**.

See [all 13 original supplements and methods scientific audit](docs/PHAPALE_2010_ALL_13_ORIGINAL_SUPPLEMENTS_AUDIT_20261016.md) and [exact remaining Phapale/He data request #31](https://github.com/rberrah/Omics2PK/issues/31). The second uploaded ZIP was a byte-for-byte duplicate of previously audited Degraeve 2024 methods, not He et al. 2022's longitudinal metabolomics.

## What is executable today

- Offline inspection of the **MultiOmics DOWNSTREAM_EXCHANGE_V1** directory: checks SHA-256, study-data access declaration, assay/sample/subject and feature-ID matching, finite numeric matrices and recorded processing scale.
- Conversion of suitable **normalized transcriptomic** matrices into the input files expected by the original SWEET implementation: expression.txt, patient.txt, gene.txt. The converter rejects unaddressed missing data, repeated subjects, incompatible namespaces and mixed tissues.
- Strict quality control of **already-computed SWEET** edge lists and extraction of descriptive network features for the proposed M1bis arm (network-wide edge statistics, per-gene signed/absolute strength). No mechanistic or CYP activity interpretation.
- Experimental training-frozen **SWEET-inspired** network proposal (distinct from upstream SWEET), with subject-level LOO/cross-fitted and held-out features, optionally separating genome-wide sample weights from targeted edge genes. Public 80-sample expression-only diagnostic, but no PK endpoint or validated causal interpretation.
- CI that executes the **real native MultiOmics exporter** on its bundled three-omic synthetic demo: 8 longitudinal subjects, 49 assays, raw RNA strictly blocked from SWEET.
- Research-only **real observed pharmacological benchmark** using 60 NCI-60 cell lines and predeclared doxorubicin GI50, with raw/PCA/network and fixed XGBoost (M1/M1bis/M2), negative shuffled-network control and tumor-type-adjusted baselines. **No observed primary MAE improvement; tumor-type-only is strongest.** Not clinical PK or verified biological mechanism.
- Additional cross-study observed GI50 study: **47 tumor-lineage-concordant NCI-60 cell lines**, 17,987 RNA features and 280 metabolites, with fixed RNA-only, metabolite-only, joint and tumor-adjusted Ridge/XGBoost comparisons. **No measured benefit of adding metabolomics to RNA in primary MAE**; one reclassified cell line explicitly excluded. These assays are NOT same-time or same-specimen multiomics.
- Public Wegler 2021/ Yin-Battail 2025 source-format inventory: bounded retrieval of PMC supporting metadata, workbook sheet names, and Zenodo record information. **No fitted CLint, donor-level link or reproduced SWEET clinical result is claimed**.
- **Thesis-locked Omics→PK evidence gate:** local, strictly fail-closed subject/time/drug/dose/concentration pairing for **RNA alone, proteomics alone, metabolomics alone, or 2–3 modalities together**; explicit cross-assay window. Strict pre-dose covariate timing prevents post-dose biological measurement leakage. Only aggregate QC is returned — **no PopPK fit or proof of clinical improvement**.
- **Wegler source-schema audit:** reads open 2021 supporting PDF and compartment-specific proteomics workbook temporarily in GitHub Actions with official PMC AWS URLs and fixed SHA-256, emits only structure/keyword counts, no donor data or invented CLint rows.
- Unit tests using synthetic data; no private patient files or cloud PK inference in GitHub Actions.

**Not implemented yet:** a *scientifically validated* inductive SWEET method (the implemented experimental method is explicitly *SWEET-inspired*), drug-specific mechanistic patient scores, Pillar-A overlays, XGBoost/PopPK fits, and clinical validation. The pinned upstream regression check reproduces raw example networks but flags an official z-score archive mismatch (issue #4). Treat the outputs as input preparation only.

## Run locally

Requires Python >=3.11.

    python -m pip install -e .
    omics2pk validate-exchange /path/to/MultiOmics-archive/exchange
    omics2pk prepare-sweet /path/to/MultiOmics-archive/exchange ./sweet-input
    # After independently running official SWEET locally:
    omics2pk summarize-sweet ./sweet-edge-outputs ./sweet-input/patient.txt ./sweet-features
    python -m unittest discover -s tests -v

Never commit the generated SWEET inputs for controlled or clinical subjects: they contain molecular measurements and subject-linked assay IDs.

## Protocol and milestones

- [SWEET network outputs and scientific caveats](docs/SWEET_OUTPUTS.md)
- [Experimental frozen-reference networks and leakage controls](docs/INDUCTIVE_SWEET_RESEARCH.md)
- [Native MultiOmics producer integration validation](docs/NATIVE_MULTIOMICS_INTEGRATION.md)
- [Public measured pharmacology evidence matrix](docs/EVIDENCE_MATRIX.md)
- [Doxorubicin observed-response negative benchmark results](docs/RESULTS_NCI60_DOXORUBICIN_20261016.md)
- [NCI-60 cross-study 47-line RNA + metabolite GI50 protocol (PR #9 merged)](docs/NCI60_TWOOMICS_BENCHMARK.md)
- [Observed cross-study NCI-60 two-omics negative benchmark results (47 cell lines)](docs/RESULTS_NCI60_TWOOMICS_20261016.md)
- [Wegler measured CLint source qualification and Yin SWEET source gates](docs/PUBLIC_PK_SOURCE_QUALIFICATION.md)
- [Thesis-aligned Omics→PK evidence ladder and future site plan](docs/OMICS_TO_PK_THESIS_PROOF_AND_SITE_PLAN.md)
- [Wegler 2021 PDF/proteome source-schema results and data availability blocker](docs/WEGLER_SOURCE_STRUCTURE_RESULTS_20261016.md)
- [**First real transcript→drug-metabolizing enzyme activity validation, 100 donor-matched livers (8 CYP outcomes)**](docs/EPA_100_DONOR_TRANSCRIPT_CYP_VALIDATION_20261016.md)
- [**First real CYP proteomics→microsomal tucatinib clearance validation, 21 donor-held-out folds**](docs/TUCATINIB_21_DONOR_MEASURED_VALIDATION_20261016.md)
- [Independent original-source audit and actual first validation gates (including mandatory RNA ABSTAIN)](docs/SOURCE_AUDIT_AND_FIRST_VALIDATION_GATES_20261016.md)
- [Qualified public Omics×PK data catalogue: in-vitro CLint, RNA functional assays, clinical TDM](docs/PUBLIC_OMICS_PK_DATASET_CATALOGUE_20261016.md)
- [Research roadmap and decision gates](docs/ROADMAP.md)
- [Comparative methods protocol: ML vs PopPK, SWEET and evidence graphs](docs/PROTOCOL.md)
- [Data security, local-only computation and publication](docs/DEPLOYMENT_AND_PUBLICATION.md)

## Publishing / licensing

The GitHub repository can be public without hosting any clinical records. New original source is offered under PolyForm Noncommercial 1.0.0 with the mandatory NOTICE attribution; authored texts are CC BY-NC-SA 4.0. See LICENSE, LICENSE-CONTENT.md and CITATION.cff. Commercial use requires separate permission from rights holders. This is source-available, **not OSI open-source**. We do not bundle private Thesis code or the upstream SWEET implementation.
