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

## P2 — public measured-data pilots
- [x] Public pharmacological pilot: mixOmics NCI-60 60 cell lines, 48 ABC transporters and **measured doxorubicin GI50** (chosen before predictions); M1/M1bis/M2 plus expression/PCA/Ridge, negative controls and histology adjustment. **No MAE improvement over mean or tumor type**; report in docs/RESULTS_NCI60_DOXORUBICIN_20261016.md.
- [ ] Qualify NCI-60 cross-study 48 syntactic / 47 lineage-concordant cell-line RNA + metabolite response benchmark, including tumor-class agreement and observed negative or positive results (draft PR #9). This is cell-line identity matching, NOT contemporaneous same-sample multiomics; no clinical PK claim.
- [ ] Yin/Battail: reproduce patient-specific network + response analysis. A measured response, **not measured PK**.
- [ ] Wegler et al.: matched liver donor proteome + measured CLint for five probe drugs; **no SWEET without RNA**.
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
