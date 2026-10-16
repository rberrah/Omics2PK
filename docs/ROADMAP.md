# Omics2PK — working to-do list (16 October 2026)

This list follows the earlier MultiOmics → SWEET → Pillar A → Pillar B R2 discussions. **When inspected the Omics2PK repository was empty; there was no pre-existing to-do file to read.** The plan below has therefore been reconstructed and must be reconciled with the owner's priorities.

## P0 — foundation and executable contract (initial PR)
- [x] Define research-only, no-patient-data boundary and M0–M4 scientific comparison.
- [x] Validate native MultiOmics DOWNSTREAM_EXCHANGE_V1 locally (SHA, patient IDs, features, missingness and omic coverage).
- [x] Prepare original SWEET input format from normalized RNA only; never infer networks without running a qualified method.
- [x] Add synthetic-only unit tests and offline CI.
- [ ] Run the actual public **MultiOmics** ZIP/exchange production output through this CLI (no copied production patient data).
- [x] Propose explicit PolyForm NC / CC BY-NC-SA licensing on the draft PR, with attribution and scope/exceptions.
- [ ] Confirm institutional/coauthor rights before merging or redistributing any third-party code.

## P1 — reproducibility and SWEET
- [x] Pin the upstream SWEET SHA (`cf064eb5f5b996c40bba3cce507384f6c0fdecb5`) in documentation; upstream code is **not** redistributed.
- [x] Implement strict SWEET edge-output reader and tabular **M1bis feature export** with synthetic tests (not network inference).
- [x] Reproduce official SWEET sample weights (80/80), one full individual raw edge set (4,753/4,753), and global summary. **Archived z-scores differ (4,753/4,753); tracked in issue #4; never claim complete equivalence.**
- [ ] Audit SWEET handling of zero variance, zeros, missing values, and the -z flag.
- [x] Prototype a clearly labeled **non-equivalent SWEET-inspired** training-frozen reference and query-only network (PR #3).
- [x] Add subject-disjoint training leave-one-out/cross-fitting and outer train/test feature construction, with negative leakage tests (PR #3).
- [ ] Validate a genuine out-of-sample extension mathematically and compare against published SWEET; current prototype uses selected genes for sample weighting instead of the genome-wide input.
- [ ] Compare stability to sample count, preprocessing, batch, tissue and threshold selection. A 60/20 public expression-only smoke test is present in PR #3; no PK endpoints.

## P2 — public measured-data pilots
- [ ] Yin/Battail: reproduce patient-specific network + response analysis. A measured response, **not measured PK**.
- [ ] Wegler et al.: matched liver donor proteome + measured CLint for five probe drugs; **no SWEET without RNA**.
- [ ] Search additional public paired RNA+PK and RNA+PD cohorts; inventory licensing, IDs and dates.
- [ ] Ask P&T whether existing **ABIS/ISBA** clinical PK records can be *legally* joined to BIOMARGIN omics for the same recipients, with appropriate timings. Do not claim BIOMARGIN itself supplies PK.

## P3 — feature representations and comparisons
- [ ] Clinical + PGx original PopPK reference and direct expression/PCA/GSVA/MOFA2 comparators.
- [ ] XGBoost full omics (M1), SWEET-only network features (M1bis), omics + SWEET (M2); matched subjects.
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
