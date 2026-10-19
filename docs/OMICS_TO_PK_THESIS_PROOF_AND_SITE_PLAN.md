# Omics2PK — thesis-locked scientific roadmap and public learning/product boundary

**Architecture source**: "Plan détaillé de thèse, architecture gelée, 15 octobre 2026" supplied by the thesis author. Public repository contains only *research interface descriptions*, not protected Pillar A/B/C/D source, patient-level records or potentially patent-sensitive implementation.

## Central research question

**Does a measured RNA, protein and/or metabolite phenotype improve PK or PK/PD prediction beyond the clinical and pharmacogenetic (PGx) reference, when sample identity, molecular collection time, dosing and PK observations are truly linked?**

**A modality is not presumed helpful.** A monotonic biological cascade (genotype → RNA → protein → metabolites → PK parameter) does *not* imply monotonically improving prediction. More features may be redundant, mediators, or confounders. A molecular correlation with *drug response* alone is not proof of improved PK prediction.

## Immutable thesis mapping

| Thesis section | Omics2PK responsibility | Gate / exclusion |
|---|---|---|
| Chapter 4: sources of biological information | Use MultiOmics research exports to form QC-qualified transcriptomic/other omic phenotypes; optional independent SWEET from RNA only. | Never force protein/metabolite matrices into SWEET as if original SWEET supported them. |
| Chapter 10: omics-informed PopPK | Test empirical incremental covariates vs existing clinical+PGx model, e.g., tacrolimus CYP3A5 and clinical covariates. | Same patients and timing, identical splits; independent PK observations and doses mandatory. |
| Chapter 13: Pillar C | Report traceable biological relevance, estimability, mediation/proxies and uncertainty for candidate omic effects. | Absent, confounded, duplicated or nontransportable mechanisms must be ABSTAIN/REVIEW, not silently accepted. |
| Chapter 12/Pillar B | Receive **qualified** covariate candidates and fit actual concentration/time/dose PK models as an external engine. | A correlation classifier or CLint assay is NOT a fitted population PK model. |
| R0/R1/R2 and D | Report evidence level and whether individual PK data exist; supply gated research evidence to the existing model qualification architecture. | In-vitro donor CLint is not a clinical R2 PopPK dataset; no automatic D/MIPD use. |

The frozen architecture remains **R0 → R1 → R2 → longitudinal individualization D**, with A/B/C as evidence, model generation and qualification responsibilities. Omics2PK supplies an **optional biological measurement branch**, not a new unreviewed master pillar.

## Evidence ladder: execute in this order

### E0. Analytical/source proof (current)
- Verify published source bytes, SHA-256, permissions and assay units.
- Explicitly distinguish measured clinical concentrations/PK parameters from in-vitro intrinsic clearance, GI50 and purely omic endpoints.
- Existing NCI-60 doxorubicin GI50 negative result: an exploratory PD proxy and **not** successful PK validation.

### E1. Biological identity and timing gate (software now in place)
- One original biological donor/patient ID; provenance of sample and species, assay compartment, RNA/protein/metabolite measurement time.
- Clinical task: **dose administration time + observed concentration time + drug identity + exact patient IDs**, not only enrollment date or "same cohort".
- Biological state covariate must be measured **before** the PK prediction landmark or explicitly modeled dynamically. The present strict baseline contract requires **pre-dose** molecular samples within a *protocol-fixed* lookback.
- Enforce time zone, technical replicates, consistent units, known concentration-to-dose event and explicit mismatches. No automated alias guessing.
- Real-world repeated dosing, infusion, metabolites and sparse TDM require a dedicated design-specific dose history and NLME contract before clinical fitting.

### E2. First genuinely measured clearance pilot (Wegler 2021, issue #11)
- The 15 donor-matched human hepatocyte / liver microsome systems, five probe drugs, protein quantities, **in-vitro CLint**.
- Complete donor × drug × preparation × method × units outcome/protein crosswalk from *original* supplementary tables. A readable proteome workbook is not sufficient.
- Donor-disjoint leave-one-donor-out, drug-specific reference and prespecified single-CYP protein predictor. Report fold errors, uncertainty, denominators.
- If original per-donor CLint is only present as figures, do not silently digitize and declare measurement truth; establish a separate low-fidelity exploratory study or obtain authors' tabular outcomes.
- The result is a **mechanistic/IVIVE link**, not proof of clinical PK prediction.

### E3. Transcriptomics-alone pilot (separate from protein CLint)
- Need a cohort with **measured RNA and matched real PK** or a well-supported RNA→protein/enzymatic CLint assay linked by subject. If unavailable, document the gap.
- Compare (1) clinical/PGx baseline, (2) expression-derived prespecified pathway scores, (3) original SWEET or independently qualified strictly inductive network scores.
- SWEET must use only TRAIN-derived reference networks and preprocessing during held-out prediction; the original transductive algorithm is not a valid test-time pipeline. Upstream archive z-score discrepancy remains issue #4.

### E4. Parallel or integrated multiomics, ONLY if data justify
- RNA + measured protein and/or metabolite phenotypes; strict sample identity and timeline alignment, no assuming simultaneous collection across public repositories.
- Fixed hypothesis-driven low-dimensional features (e.g. CYP3A activity, ABCB1 transporter-related proxies, inflammation state) before high-dimensional graphs/ML.
- Prespecified comparators and incremental Δerror / held-out calibration. Evaluate whether proteomics or metabolites add independent information rather than automatically enriching a larger predictor.
- Penalize redundant PGx→RNA→protein→metabolite cascade effects; test negative controls and remove proxies if indistinguishable.

### E5. Clinical R2 demonstration (BIOMARGIN only if access and linkage qualify)
- Demonstrate actual alignment of transcriptomic/omics specimen, patient, collection time, formulation, complete doses, measured plasma concentrations and/or observed AUC for tacrolimus or MPA, with ethics/legal access.
- Baseline: Woillard/alternative clinical + PGx PopPK model; estimate and predict on identical outer patient/site-held-out splits. Then add prespecified biological phenotype candidate scores; compare error, uncertainty, bias, transportability and subgroup stability.
- Do not claim BIOMARGIN supplies PK merely because a BIOMARGIN omic patient record exists. Confirm matching with an authorized PK/EHR data source; keep all protected individual data off GitHub and the static site.
- Evidence from E2 or E3 does **not** substitute for this clinical validation.

### E6. Thesis A/B/C qualification and external lock
- Freeze source provenance, feature selection, patient split, scaling, confounder controls, reference model and population/time validity.
- B may test PK effect sizes as candidate parameters; C qualifies direction, magnitude, biological consistency, redundant mediators, uncertainty and transportability. Explicit abstention is a valid answer.
- Independent cohort, temporal/site transport, source ontology and model version controls, no patient-level artifacts committed.
- Only qualified analysis outputs can be considered for later individualization in D.

## Clinical benchmark labels — avoid semantic drift

| Thesis Chapter 10 arm | Input | Purpose |
|---|---|---|
| Reference | Clinical PK + clinical covariates + available PGx | Incremental-effect baseline |
| M2 | Reference + independently qualified RNA/SWEET state scores | Transcriptomic incremental effect |
| M3 | M2 + *parallel* protein / metabolite state measurements | Additional modality signal |
| M4 (exploratory) | Qualified integrated pathway or network-aware multiomics | Complexity justified only if sample size supports it |

For method ablations, additionally compare reference+direct RNA scores **without SWEET**. Prefer a few predeclared functional features to thousands of network edges.

## Site integration: **after empirical qualification**, not a substitute for it

Source of truth: Omics2PK repository is a research engine/contract; MultiOmics performs exploratory/confirmatory omic QC and pathway analysis; rberrah.github.io is a static entry point, learning interface and approved public demo. Do not insert a clinical predictive service into a static GitHub Pages site or collect clinical omics/PK data on third-party endpoints.

### Future section: Ecosystem → Omics and pharmacometrics
Proposed concise path: **/omics2pk/** (final route to be reconciled with active site navigation, not introduced until gates E2–E4 are demonstrably complete).
- FR/EN, light/dark, consistent header/logo, margins, Resources and Ecosystem, return link to rberrah.github.io, accessible navigation and research-level warnings.
- Courses: (1) DNA vs RNA vs protein vs metabolite; (2) linking phenotype to CYP/transporter/organ activity; (3) CLint vs plasma CL/CL-F and IVIVE; (4) dose–concentration and PopPK; (5) time/provenance and patient identity; (6) baseline PGx vs omics incremental value; (7) SWEET network inductive vs transductive; (8) confidence, confounding, missingness and C-level abstention.
- Beginner “discover in one click” demos should use reproducible **public or explicitly synthetic** data and report what was and was not measured, uncertainty, provenance and version.
- Professional tool: validated input schema → QC report → profile/pathway scores → optional qualified SWEET → exported **PK covariate research package** (not automated dosage advice). Actual PopPK estimation stays in an authorized scientific backend only.
- No deployment claim, publication claim or clinical predictive benefit until the corresponding evidence gate is passed. For sensitive Thesis A/B/C intellectual property, publish only reviewed public interfaces after institutional IP checkpoint.

## Success criteria and negative-result discipline
A valuable result may be: PGx sufficient; measured transcriptome not better; SWEET no incremental gain; protein CLint only in vitro; no usable aligned BIOMARGIN PK measurements. All remain scientifically publishable if predeclared, traceable, externally tested as appropriate and reported without overstating clinical relevance.
