# Hybrid Fusing — shared scientific philosophy, not a shared algorithm

**Date:** 2026-11-10  
**Status:** CONCEPT / RESEARCH DESIGN ONLY — not implemented, not prespecified, not validated.  
**Home:** Omics2PK. **Possible later consumers:** NetPK-style model adaptation and MetaPK in Thesis.  
**Immediate priority:** maximise rigorous H1–H5 evidence from existing **public data**, not funding preparation.

## 1. Purpose

Investigate whether mechanistic knowledge, observed predictive information and individual context can jointly determine **what to admit, what to select, what to adapt, and whether to fuse**. This is the common idea behind:

1. **NetPK (tacrolimus) and rezafungin adaptation:** assess pharmacological plausibility and differences between source and target populations, then test which covariates, effects or parameter distributions merit adaptation.
2. **Omics2PK:** identify molecular measurements and representations that add information *conditional on* clinical and pharmacogenetic predictors; use complementary modalities only when matched in the same people.
3. **MetaPK:** reject inapplicable published PK models, then decide whether to select a model, recalibrate it, or combine predictive **distributions**, while retaining model disagreement and uncertainty.

These are **different mathematical operations**. A common decision contract is plausible; a universal numerical weighting formula is not assumed.

**Central falsifiable question:** does a mechanistically informed, context-aware admission and selection rule followed by conditional integration outperform a transparent rule-based gate plus a simple well-calibrated statistical benchmark?

The word **Fusing** is a working label, not an originality or performance claim. Stacking, Bayesian model averaging, mixtures of experts, sparse regression, model selection and mechanistic priors are mandatory conceptual comparators.

## 2. One decision sequence, three integration operators

**(A) QUALIFY / FILTER — hard eligibility before prediction.**

Check source integrity, molecule, analyte, timing, biological unit, tissue, population domain, formulation, assay/measurement reliability, covariate availability, units, endpoint and missingness. Reject impossible or unidentifiable joins, incompatible models and circular targets. An evidence gap is reported as *unknown* or *not evaluable*, never automatically treated as biological absence.

**(B) RANK / SELECT — incremental information after the reference.**

Compare candidate evidence on three distinct axes:

- Mechanistic support: strength, direction, pathway, target parameter, compartment, provenance, contradictions and applicability to the endpoint.
- Statistical contribution: **out-of-sample incremental improvement conditional on a baseline**, not marginal association with the raw outcome.
- Reliability: sample size, stability, replicability, assay quality, calibration, redundancy, cost and sensitivity to nuisance variation.

Mechanistic plausibility is a **soft prior for admissible candidates**, not permission to keep a candidate with no demonstrated benefit. Conversely, a statistically useful candidate without an established mechanism is retained as an explicitly labelled exploratory hypothesis. Neither score is silently converted into a probability or a PK parameter prior.

**(C) INTEGRATE / FUSE / ABSTAIN — the simplest valid operator.**

Possible outcomes are: baseline only; one selected candidate; sparse additive covariate effect; parameter/distribution adaptation; fusion of molecular representations; or a weighted mixture of **predictive** distributions. Integration is not obligatory. A negative, unstable or poorly calibrated fusion must revert to a simpler comparator or abstain.

### Operator boundaries

| Domain | Object selected | Permitted integration | Common category error |
|---|---|---|---|
| **NetPK tacrolimus / rezafungin** | Covariate, causal/physiological hypothesis, parameter effect, source-to-target shift | Covariate on CL/F/V/absorption, mechanistic module, distributional recalibration with identified assumptions | Treating a physiologically plausible covariate as a proven causal effect; mixing effects from incompatible models |
| **Omics2PK** | Gene, protein, metabolite, activity proxy, pathway or multi-omic representation | Sparse biomarker panel, externally fixed functional score, conditional molecular predictor or fold-safe multimodal integration | Combining unrelated donors, choosing markers on the test endpoint, confusing in-vitro function with patient PK |
| **MetaPK** | Published PK model and its domain of applicability | One best applicable model, calibrated model selection, or mixture of predictive distributions with context-dependent weights | Treating similarity of a model's development population as evidence of its predictive competence; averaging incompatible parameter vectors |

Examples of the shared **reasoning**, not reusable code: tacrolimus CYP3A5/haematocrit/food mechanisms; albumin and illness-related differences for rezafungin; an omic signal that improves on a CYP activity-score baseline; selecting a PK model appropriate to age, care setting, renal function and available covariates.

## 3. Minimal formalisation

For context x, candidates j and endpoint y:

1. **Eligibility:** G_j(x) is a documented gate (eligible / rejected / unknown with a reason), applied before endpoint-based model selection.
2. **Reference:** B(x) is the prespecified clinical/PGx or PK-model baseline. It must be identical across fair contrasts.
3. **Incremental evidence:** Delta_j = Loss_OOF(B) - Loss_OOF(B + j), with all residualisation, preprocessing, selection, tuning and recalibration **inside training folds**. Positive Delta means improvement. This is measured evidence, not a mechanistic score.
4. **Mechanistic evidence:** M_j records signed, source-grounded biological applicability and uncertainty; it may rank candidates *before* fitting but cannot override held-out failure.
5. **Complementarity:** evaluate j conditional on already selected information, including within-family technical redundancy and cross-compartment mismatch. Performance from different cohorts does not establish complementarity.
6. **Decision:** choose a sparse admissible set, a calibrated integration operator, the baseline alone, or an explicit **NOT_EVALUABLE** result.

For a valid predictive-distribution mixture in MetaPK:

**p_fused(y | x) = sum over eligible m of w_m(x) × p_m(y | x)**, where weights are non-negative and sum to one. Its variance has **within-model** and **between-model** components; both must be reported. Fusing densities is **not** fusing structural ODEs, covariate equations or PK parameter distributions. Model recalibration and uncertainty repair are separate operations; mixing undercovered distributions does not guarantee calibration.

No universal formula such as w = mechanistic_score × association_score is authorised by this document. Any proposed weighting or gating algorithm must be benchmarked against simpler alternatives.

## 4. Empirical warnings already learned

- **MetaPK:** density-based development-population coverage weights (pi_m) failed against a simple predeclared applicability-rule filter in the distinct-domain MIMIC vancomycin library. **Population density is not predictive model competence.** Do not reinvent pi_m under the Fusing name. This result is recorded in Thesis, including its limitations.
- **Omics2PK H2:** STRING/RWR/CollecTRI mechanism-guided candidates do not automatically outperform direct transcripts or statistical selection; a mechanistic network can be wrong, incomplete or uninformative.
- **Omics2PK H1:** ARIANA metoprolol has an exploratory gain from targeted solanidine metabolites; venetoclax/4beta-hydroxycholesterol (n=12) has no detected gain. Marker utility is drug- and endpoint-dependent.
- **Omics2PK H4:** wide representations may fail under severe p/n; a positive permutation null and training-only differential shrinkage can mimic gain. Test conditional selection against **baseline residuals**, rather than rewarding predictors correlated with what the baseline already explains.
- **Omics2PK H5:** tissue/preparation and sampling context matter. Wegler's encouraging hepatocyte findings were not confirmed in a subsequent 15-donor independent replication; uncertainty that refits models must be distinguished from bootstrap of fixed out-of-fold errors.
- **EATRIS-Plus:** genuine multi-omic matching supports methodological tests, **not** pharmacological or clinical-PK claims in the absence of a drug endpoint.
- **No clinical-PK confirmation:** the public evidence currently contains exploratory clinical successes and failures, not independent confirmation of a multi-omic clinical benefit.

These observations motivate **testing** a hybrid selector. They do not establish its superiority.

## 5. Mandatory falsification programme — within H1–H5, no new H6

**H1 — incremental information.** Compare the same matched individuals, endpoint and folds with baseline clinical + genuine PGx, statistical molecular selection, mechanistic-only selection and the hybrid. New clinical data sources or biologically independent cohorts take precedence over further tuning on previously inspected cohorts.

**H2 — mechanism-guided selection.** Test a genuine hybrid against: (i) statistical-only selection, (ii) network-only selection, (iii) simple externally fixed biology score, (iv) univariate/elastic-net reference, and (v) randomly rewired, direction-shuffled or deliberately wrong mechanism priors. An advantage must exceed those nulls and survive independent evaluation.

**H3 — parsimony.** Compare gain-versus-panel-size and k80/k90 **only where an actual baseline gain exists**. Include the cost of obtaining the molecular measurement and stability of selected features. A feature that merely reproduces the baseline is not valuable incremental information.

**H4 — reliability region.** Vary p/n, signal-to-noise, missingness, redundancy, bias, correlation and number of irrelevant features. Compare raw-target selection against fold-safe residual/conditional selection, including the selection of nuisance correlates. Refit the **entire** pipeline under permutation. Correct leakage in imputation and retain prior frozen results unchanged.

**H5 — biology and transportability.** Check molecule, compartment, tissue, dosing/sampling clocks, quantitative versus semi-quantitative assay, functional proxy and temporal stability. Do not transfer an enzyme-activity claim to clinical AUC/CL/F without same-subject evidence.

For **PK-model distribution fusion** (future MetaPK evaluation), additional comparators are required: best applicable model, explicit rule-based filter, uniform averaging after filtering, stacking, BMA where its assumptions are met, conditional mixture-of-experts and calibration-only repair. Compare log-score, point error, calibration/coverage **and interval width** on held-out people/sites. Ensure models and weights never use the evaluated concentration unless the task expressly allows past-patient TDM and the prediction landmark is later.

### Rules for every experiment

- Freeze questions, candidate universe, preprocessing, folds, endpoints, losses and success/failure criteria before evaluating outcomes.
- Split by **biological person**; stratify or hold out laboratory/site and time when required. Never join different cohorts into pseudo-subjects.
- Fit baseline, residuals, feature selection, weights, PCA, imputation, recalibration and model tuning inside training folds. Cross-fit residuals if used as inputs for a subsequent selection stage.
- Report uncertainty that incorporates the model-refitting procedure where feasible, multiplicity, calibration and negative controls. Non-detection is not equivalence.
- A mechanism prior must not be learned from the held-out endpoint; a clinical gate cannot be chosen from observed errors.
- Separate **selection**, **conditional covariate/parameter adaptation**, **representation fusion** and **predictive-distribution fusion** in every result.
- Never claim originality merely from the Fusing label or an ensemble that reproduces existing stacking/mixture-of-experts behaviour.

## 6. Minimal interfaces for a *future* common framework

This is a **conceptual data contract**, not a request to create a shared library.

**CandidateRecord**
- Stable ID, kind (covariate / omic marker / functional score / PK model), original source and checksum, drug/endpoint, units, specimen and biological unit.
- Mechanism, signed evidence and conflicting evidence; applicability conditions, required inputs and missing-data behaviour.
- Evidence tier (published / reproduced / exploratory / independently tested), assay quality and limitations.
- For PK models only: predictive distribution, original model assumptions and information on residual, interindividual and parameter uncertainty.

**DecisionRecord**
- Prediction landmark, eligible set with exclusion reasons, reference model and evaluation context.
- Selected variables or models and the reason for each decision; performance and incremental contribution estimated without leakage.
- Integration operator and parameters/weights, uncertainty decomposition, calibration, fallback and NOT_EVALUABLE reasons.
- Source fingerprints, protocol version, test cohort fingerprint and full provenance.

Possible shared layer: **evidence + applicability + decision audit**. Keep distinct domain-specific operators for covariate effects, omic integration and predictive mixtures.

## 7. Sequencing and repository boundaries

**Now — Omics2PK:** continue its public-data H1–H5 programme and audit existing MultiOmics/Omics2PK selectors and ablations. Propose the **smallest falsifiable hybrid experiment**; do not add an architecture because the concept note exists.

**Later — cross-project review:** compare tacrolimus NetPK, rezafungin adaptation and the MetaPK rule-filter/weighting failures to Omics2PK's actual experiments. A separate transversal review agent may then decide whether a truly shared mathematical decision engine is supported.

**Do not:** import Thesis/MetaPK/NewPK code into Omics2PK; alter frozen studies; claim a new H6; present simulations as biological validations; infer efficacy from a graph prior; or pursue funding design as the present scientific objective.

## 8. Decision gate for implementation

Implement only after a written protocol identifies (a) a real eligible public dataset or a strictly bounded synthetic falsification test, (b) at least one simple baseline that could win, (c) a distinctive mechanism beyond existing methods, and (d) a measurable failure condition.

**Default outcome until such a test exists: concept retained, no new production code.**
