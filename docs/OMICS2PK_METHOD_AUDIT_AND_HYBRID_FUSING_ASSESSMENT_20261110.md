# Method audit, a minimal definition of Hybrid Fusing, and an assessment of its novelty

**Version 1.0.0 — 2026-11-10.** Evidence type: **audit and design, nothing executed**.
**No code was imported from Thesis, NetPK, NewPK or MetaPK.** Those repositories
were read for conceptual comparison, as the commissioning note asks, and
nothing from them is copied, adapted or depended upon.
**No frozen protocol or result is modified.**

---

## 0. The assessment, stated first

**Three of the four components of Hybrid Fusing already exist, two of them
inside this ecosystem, and the one that is genuinely new has a close precedent
that failed.**

| Component | Status | Where |
|---|---|---|
| mechanistic selection, outcome-blind | **exists and is tested** | `mechanistic_graph_selection.py`; H2a |
| **hybrid selection: knowledge proposes, data disposes** | **exists, is tested, and does not beat a univariate filter** | `hybrid_mechanistic_selection.py`; H2d |
| representation fusion across layers | **exists**, for non-PK endpoints | MultiOmics: MOFA2, DIABLO |
| prediction fusion with constant weights | **exists** as `stacking` | MetaPK |
| prediction fusion with per-individual weights | **exists** as `bma` | MetaPK |
| prediction fusion with weights from **covariates** | does not exist here; **is mixture-of-experts** | literature, 1991 |
| **weights informed by a mechanistic prior** | **does not exist anywhere in the ecosystem** | — |
| a declared fallback to selection | exists in spirit (abstention) | `hybrid_mechanistic_selection.py` |

**So the defensible novelty claim is narrow:** *a prior on the fusion weights
derived from mechanism, with a committed fallback to selection.* Everything
else is re-description.

**And that narrow claim has a precedent that failed, measured.** MetaPK's π_m
weighted published PK models by the density of their design population at the
patient — a knowledge-derived weight needing no target data. Against a plain
rule filter it lost by **Δ log-score −0.855, 95 % CI [−0.887, −0.824], 100 % of
bootstrap draws below zero**, with calibration 14 points worse. The diagnosis
transfers directly: *a density varies as 1/σ, so it rewards a **concentrated**
design population, not an **applicable** model.* The omic analogue is immediate
— **network proximity measures how well studied a gene is, not how predictive
it is** — and this project already knows it, which is why
`mechanistic_graph_selection` carries degree strata in its permutation null.

**Recommendation, before any implementation:** run the headroom gate in
`OMICS2PK_FUSING_REFUTATION_PROTOCOL_20261110.md`. If it does not clear, the
honest conclusion is to adopt **stacking with a rule filter** as the
comparator *and* the solution, which is what MetaPK's own evidence already
points to.

---

## 1. Audit of Omics2PK — 66 modules, and a category that was empty until 2026-11-11

### 1.1 Selection, richly developed

| Module | What it does |
|---|---|
| `mechanistic_graph_selection.py` | the executable chain: STRING v12 at combined score ≥ 700, column-stochastic transition, **random walk with restart** (restart 0.30, tol 1e-11), seeds from the endpoint label, from Reactome pathways or from CTD substrate relations, then a **permutation FDR with degree strata** (10 strata, 500 resamples) so a hub is not selected for being a hub |
| `hybrid_mechanistic_selection.py` | knowledge shortlists `3k`; training donors retain by **absolute partial rank correlation residualised on the baseline**, calibrated against a permutation null at 5 % FDR; **abstention allowed**; agreement with the statistical arm reported as containment, Jaccard **and** exact equality, because one number conflates them |
| `functional_omics.py` | signed CollecTRI transcription-factor activities as features (F-series) |
| `inductive_sweet.py`, `frozen_sweet_pipeline.py`, `sweet_features.py` | a SWEET-inspired inductive graph, explicitly **not** the original method |
| comparators | univariate Spearman filter, LASSO entry order, XGBoost gain (`hlc_real_xgboost_validation.py`), elastic net, repeated random selection, full panel |

### 1.2 Evaluation machinery, unusually complete

Fold-safe preprocessing (`transcriptome_preprocessing.py`), donor-disjoint
repeated nested cross-validation, leave-one-donor-out, paired bootstrap over
the resampling unit, **permutation nulls that refit the model**, permutation
FDR without a resolution floor, Bonferroni within declared families, a
**refitting out-of-bag interval** (`plasma_biomarker_validation.py`), Gram-matrix
PCA for fold-safe directions at p ≫ n, and **fourteen modules whose names declare them gates, contracts or
readiness checks**, all fail-closed: they refuse rather than approximate.

### 1.3 Fusion: the category was empty, and now holds two members — both of which lose

**As audited on 2026-11-10:** searched across `src/`, **`stacking` 0,
`ensemble` 0, `model_averag` 0, `mixture` 0, `meta_learner` 0.** The
twenty-five occurrences of "fusion" were `hepatic_protein_activity_fusion.py`,
where the word means *joining two measurements on the same donors* — a
data-availability question — and whose own docstring says "'proteins+activity'
is **NOT** multi-omics".

**Omics2PK compared arms and had never combined them.** That was the gap a
Fusing proposal would fill, and it was a real gap.

**Updated 2026-11-11.** `src/omics2pk/fusion_headroom_diagnostics.py` now
contains **one** member of the category, and the distinction matters:

| | Status |
|---|---|
| `nested_constant_stacking` — **constant weights fitted on one partition and scored on another** | **this is stacking**, implemented. It exists to measure what stacking *recovers*, and it is reported as a diagnostic, but it is not merely a ceiling and should not be described as one. **Amended 2026-11-11:** it is nested only in the prediction layer. Because that layer was computed once over the whole cohort, the weights were fitted on individuals whose base-model training sets contained the rows the weights were then scored on. `nested_fusion_evaluation` refits the arms inside every outer fold and supersedes it; measured, the difference is a factor of about **2.9** on the one contrast where the old version found a gain |
| `nested_fusion_evaluation` — `S`, `U`, `K` and `X` on shared outer folds, arms refit inside each | the category now holds a **contextual gate** as well, and it does **not** beat stacking on HLC: 0 of 30 cells reach the declared margin. A mechanistic prior on fusion weights still does not exist, and under the protocol's own condition it is not authorised |
| `constant_weight_oracle` | a hindsight ceiling, not a method |
| `individual_convex_oracle`, `individual_selection_oracle` | hindsight ceilings, not methods |
| a **gate**, w(z_i) | ~~**does not exist**~~ **built and tested 2026-11-11** in `nested_fusion_evaluation.py`, a softmax over ridge-predicted per-arm losses. It reaches the declared margin on **0 of 30** cells on HLC, so the category is no longer empty and no longer promising |
| a **mechanistic prior on weights** | **does not exist** |

So the honest statement is now: *comparator `K` of the protocol's arm list is
built and measured; arms `X` and `F` are not.* The executed figures are in
`OMICS2PK_FUSION_HEADROOM_RESULTS_20261111.md`, where constant-weight stacking
recovers a detectable gain on 3 of 13 cells and **degrades by 46 % on one**.

### 1.4 What the audit implies about H2

The hybrid rule the commissioning note proposes to build is, in its selection
half, **already built** — including the detail that yesterday's sizing
simulation independently re-derived as a requirement: selection on the
**residual after the baseline** rather than on the raw outcome
(`partial_rank_strength`, with the rank transform applied *before* the
residualisation precisely so that "a variable that merely repeats a baseline
covariate" cannot keep an unearned partial correlation).

The corroboration is worth stating: two independent routes — a destroyed-signal
experiment on HLC in 2026-11-08, and a power simulation on 2026-11-10 — reached
the same design rule. **The Axis A and cohort-sizing pipelines do not follow it**
and select on the raw outcome; that is a known inconsistency, recorded here, and
it is not fixed by editing frozen runs.

---

## 2. Audit of MultiOmics — representation fusion already exists

`/Site PK/Base/MultiOmics/multiomics-engine/`, R, public repository.

| Method | Implementation |
|---|---|
| **MOFA2** | unsupervised multi-block factor analysis. `mofa2_input_preflight` **refuses clearly underpowered sample counts** and refuses a **silent intersection of unmatched subjects**, "otherwise the effective population can depend on which views were supplied" |
| **DIABLO** (mixOmics multiblock sPLS-DA) | supervised multi-block, `keepX` grid tuned **inside repeated CV before the final fit**, with `run_diablo_outer_holdout`, `diablo_selection_stability` and `prepare_diablo_holdout_split` |
| differential | limma, DESeq2, edgeR, voom, with `compare_rnaseq_methods` |
| longitudinal | lme4, with a random slope attempted **only when the design supports it** |
| other | fgsea, Cox |

Its stated policy is directly relevant: **"sensitivity methods are compared,
never merged by averaging p-values/q-values."**

**Consequence for the proposal.** "Fusion of several complementary
representations" is **not a new idea in this ecosystem** — MOFA2 and DIABLO are
exactly that, with selection stability and outer hold-out already wired in.
What they are not applied to is a **pharmacokinetic or enzyme-activity
endpoint**, and they carry **no mechanistic prior**. Any Fusing experiment on
omic representations must therefore include **DIABLO as a comparator**, or
explain why a supervised multi-block method with fold-safe tuning is not the
relevant baseline.

---

## 3. Conceptual comparison with MetaPK — read, not imported

`rberrah/Thesis`: `src/pkdt/exposure/metapk.py`,
`docs/these/chapitre15_metapk_resultat_negatif_pi_m.md`,
`docs/METAPK_R1_NOTE_COUVERTURE.md`, `scripts/metapk/*`.

### 3.1 What it already implements

Four schemes over a library of published population-PK models, **one of which
is not a method**:

| Scheme | Definition |
|---|---|
| `uniform` | equal weights — **the control, known to lose** |
| `bma` | weights ∝ each model's fit to **this patient's earlier levels**, per prediction, on the Laplace evidence |
| `stacking` | one weight vector maximising the predictive density of the held-out level; concave in the weights, solved by **EM**; accepts **per-observation importance weights**, which is how a copula correction for "the monitored are not the treated" enters |
| `oracle` | the best model per prediction, **using the answer** — the ceiling, not a method |

And the quantity passed downstream is a **width, not a point**, by the law of
total variance:

> Var[AUC] = Σ w_m σ_m² (within) + Σ w_m (μ_m − μ̄)² (**between**)

The between-model term "is the one a single model cannot produce and the one
that matters": four vancomycin models span **3.25-fold** on AUC for one
reference patient.

### 3.2 The three measured lessons that constrain any Fusing design

**(i) There is a ceiling, and gains must be reported against it.** On 29,349
held-out vancomycin predictions the best model *per prediction* was roberts
34.32 %, pkjust 34.23 %, goti 31.45 % — no model dominates. The oracle reaches
MdAPE **11.14 %** against **15.68 %** for the best single model: **4.54 points
is the entire recoverable headroom**, and "any gain must be reported against it
rather than against zero."

**(ii) Averaging without weighting degrades.** The geometric mean of the same
three models gives MdAPE **16.17 %** — *worse* than the best single model. This
is why `uniform` exists as a control rather than as a candidate.

**(iii) A knowledge-derived weight lost to a rule filter, badly.** π_m against
F: **Δ log-score −0.855 [−0.887, −0.824]**, 100 % of draws below zero,
90 %-coverage 0.715 → 0.573. What **held** was the protection — the preterm
model got weight exactly zero. What **failed** was the weighting: π_m gave 0.51
of the weight to a neurosurgery model that was running *without two of its
covariates*, because its published margins (mean ± SD) are narrower than a
competitor's printed ranges. **A density rewards a concentrated design
population, not an applicable model.**

Its decisions were: π_m stops and is **not** repaired on the evaluation data; the
**rule filter becomes the default**; a *quantile*-based applicability measure is
registered **without a run**, to be tested only prespecified and on a different
target.

Also recorded: on that library, adding two further models "n'a rien apporté sous
aucune pondération mesurée", and the log-score **penalises an absurd component
weakly** — `uniform` had the best L5 log-score while producing a 546 mg/L
interval, its 0.99 coverage "bought by width". **Any Fusing evaluation must
therefore score calibration and interval width, not a log-score alone.**

### 3.3 What MetaPK does **not** have

1. a **mechanistic** axis in the weights — its axes are predictive fit and
   population coverage;
2. **mechanistic fusion of the PK equations** rather than of their predictive
   distributions;
3. any application to **omic** predictors.

Those three are the only places a new method has room, and the commissioning
note already distinguishes (2) from distribution fusion, correctly.

---

## 4. A minimal mathematical definition of Hybrid Fusing

The note's key requirement is to distinguish **selection**, **fusion of
representations** and **fusion of predictions**. They are different objects and
conflating them is how an architecture becomes unfalsifiable.

### 4.1 Objects

For individuals *i = 1…n* with outcome *y_i*, a baseline design *b_i* (clinical
covariates and genotype), candidate features *x_{ij}*, *j = 1…J*, and a
**context** vector *z_i* (available covariates, tissue, sampling time,
assay provenance — anything observable **without** the outcome):

- a **mechanistic prior** π_j ∈ [0,1], computed **without the outcome**;
- a **conditional contribution** γ_j, measured on training folds **after
  residualising both sides on the baseline**;
- base predictors *f_m*, *m = 1…M*, each a (representation, model) pair.

### 4.2 The three operations, defined separately

**A — Filter.** 𝒜 = { j or m : admissible }, where admissibility depends on
drug, tissue, compartment, analytical quality and required covariates, and is
**independent of y**. This is MetaPK's rule filter F, and it is the only
component with a measured success behind it.

**B — Select.** A map 𝒮 : (π, γ) ↦ S ⊆ 𝒜. The existing hybrid is the special
case
> S = { j ∈ shortlist(π, 3k) : γ_j exceeds a permutation-FDR cutoff }, |S| ≤ k,
> with S = ∅ allowed.

**C₁ — Fuse representations.** A map *R* : (x_{i,S}) ↦ u_i ∈ ℝ^d, d ≪ |S|,
fitted on training folds. MOFA2 and DIABLO are instances; Gram-matrix PCA as
used in Axis A is an instance with d = 2 or 5.

**C₂ — Fuse predictions.** ŷ_i = Σ_m w_m(z_i) f_m(x_i), with Σ_m w_m = 1,
w_m ≥ 0, and *w* fitted on training folds only.

### 4.3 Where each existing method sits

| Method | w |
|---|---|
| selection | w = e_{m\*}, m\* chosen on training folds |
| uniform averaging | w_m = 1/M |
| stacking | w constant, fitted to maximise out-of-fold predictive density |
| BMA | w_m(i) ∝ p(data_i \| m), i.e. depends on **individual i's own observations** |
| mixture of experts | w(z_i) from a gating function of **covariates** |
| **Hybrid Fusing, as proposed** | w(z_i; π), a gating function of covariates **with a mechanistic prior**, plus a declared fallback to w = e_{m\*} |

### 4.4 The one formal statement that makes it testable

Write L_{im} for individual *i*'s loss under arm *m*. Then:

- **best single arm** L\* = min_m (1/n) Σ_i L_{im};
- **selection ceiling** (per-individual oracle) L_sel = (1/n) Σ_i min_m L_{im};
- **combination ceiling** L_cvx = min_{w ∈ Δ} (1/n) Σ_i ℓ( Σ_m w_m f_m(x_i), y_i ),
  the best achievable by *any* convex combination with weights chosen **knowing
  the answers**.

L_sel bounds every method that picks one arm per individual. It does **not**
bound averaging: if arm A errs +1 and arm B errs −1, the oracle's absolute
error is 1 and their average's is 0. L_cvx is the bound that covers stacking,
BMA, mixture of experts **and** Hybrid Fusing.

> **If L\* − L_cvx is not detectably greater than zero, no member of the
> fusion family can improve on choosing one arm, and Hybrid Fusing is refuted
> on that cohort and endpoint without being implemented.**

That is the refutation test. Its protocol, its comparators and its verdicts are
in `OMICS2PK_FUSING_REFUTATION_PROTOCOL_20261110.md`.

---

## 5. Critical analysis of novelty

### 5.1 What is not new

| Claimed element | Prior art |
|---|---|
| knowledge proposes, data disposes | **this project's own** `hybrid_mechanistic_selection.py`, tested on 9 endpoints × 8 panel sizes |
| conditional contribution over a baseline | same module, `partial_rank_strength` |
| abstention when nothing survives | same module |
| fusion of complementary representations | **MOFA2 and DIABLO**, in the sibling MultiOmics repository, with fold-safe tuning and selection stability |
| one weight vector fitted out-of-sample | **stacking** — in MetaPK, by EM |
| weights depending on the individual | **BMA per prediction** — in MetaPK |
| weights depending on covariates | **mixture of experts** (Jacobs, Jordan, Nowlan & Hinton, 1991) |
| fall back to a single model when fusion does not help | the standard practice stacking already permits, since w = e_m is in the simplex |
| discarding inapplicable models by declared rules | **MetaPK's rule filter F**, which beat the knowledge-derived weighting |
| passing a distribution rather than a point | MetaPK's within/between decomposition |

### 5.2 What is new, stated as narrowly as the evidence allows

1. **A mechanistic prior on fusion weights.** No implementation in either
   repository weights a *predictor or model* by mechanistic plausibility. In the
   literature the nearest neighbours are informative-prior BMA and
   knowledge-guided ensembles; within this ecosystem it is unoccupied.
2. **A single object spanning omic predictors and PK models.** MetaPK fuses
   models; MultiOmics fuses representations; nothing fuses across both. Whether
   that unification is *scientifically* useful or merely architecturally tidy is
   exactly what must be tested, and tidiness is not a result.
3. **Mechanistic fusion of PK equations** rather than of their outputs. This is
   genuinely different in kind from everything above — and it is also the item
   with no data to test it, since it would require the structural models
   themselves to be composable and a cohort on which to compare compositions.

### 5.3 Risks, in order of how likely each is to sink the project

**(1) Re-labelling.** The largest risk is that Fusing becomes stacking with a
gating function and new vocabulary. **Mitigation:** the protocol requires
stacking, BMA and mixture-of-experts as *fitted comparators on the same folds*,
not as citations, and a prespecified statement of which parameter makes Fusing
differ from each.

**(2) The π_m failure mode, transposed.** A mechanistic prior can measure
**study intensity** rather than relevance, exactly as a coverage density
measured concentration rather than applicability. Network degree is a proxy for
how much a gene has been published about. **Mitigation:** the prior enters
**only** through a degree-stratified, outcome-blind score — the one
`mechanistic_graph_selection` already computes — and a **degree-permuted prior**
is a mandatory negative control.

There is a second, independent instance of the same hazard already measured in
this project, and it argues the mitigation is not optional. Selecting candidate
features by their correlation with the **outcome** makes the informative
fraction of the selected set rise to **1.0** at redundancy 0.85 — its highest
value exactly where the layer is least useful — while selecting on the
**residual after the baseline** gives **0.234** there, which is the honest
figure. A criterion can be most confident where it is least informative. A
mechanistic prior is a criterion of that kind.

**(3) Unfalsifiability through the fallback.** An architecture permitted to fall
back to "select one arm" can never lose, and therefore can never be shown to
win. **Mitigation:** the fallback decision must be taken **inside the training
folds** and the committed choice evaluated out of sample. A method that decides
after seeing the test loss is not a method.

**(4) Degrees of freedom.** Filter, prior, shortlist width, FDR level, gating
family, fallback rule: six places to tune. At n = 47 this is enough to fit the
test set by iteration. **Mitigation:** every one is frozen in the
prespecification, and the project's existing discipline — a permutation null for
*every* reported contrast, a refitting interval, declared family size — applies
unchanged.

The refitting interval matters more for a fusion method than for anything built
so far, and the size of the correction is measured: on six contrasts a
refitting out-of-bag interval was **1.78** to **2.09** times wider than the
stored-error one for a **one-predictor** arm, and **3.76** to **6.32** times
wider for a **two-predictor** arm. A fusion method has *M* weights to estimate.
Whatever the right factor is at *M* arms, it is larger than at two, and a
stored-error interval would make a fusion gain look far more certain than it is.

**(5) Sample size.** A gating function w(z) needs enough individuals to
estimate per-arm competence as a function of z. The available clinical cohorts
are n = 47 and n = 12. **Mitigation:** the only cohorts where a gating function
is estimable at all are HLC (≈397 donors, 9 endpoints) and EATRIS-Plus (79–114
individuals, 12 layers). The clinical cohorts can test the **selection**
fallback and nothing more.

**(6) The endpoint is scarce, and reuse costs status.** Every cohort that could
host a Fusing experiment has already been analysed by this project, so every
result is capped at exploratory. **Mitigation:** none available. It is a limit
to state, not to engineer around.

### 5.4 The honest prior on the outcome

Stated before any run, so it cannot be adjusted afterwards:

- **Most likely:** the headroom gate shows a small but non-zero combination
  ceiling on HLC, constant-weight stacking recovers most of it, and the
  mechanistic prior and the gating function add nothing detectable. Fusing then
  reduces to stacking, and the correct report is that stacking is the method.
- **Second most likely:** no detectable headroom anywhere, and the whole family
  is refuted on the available data.
- **Least likely, and the only outcome that would justify the architecture:** a
  gating function on context beats constant-weight stacking **and** a
  degree-permuted prior performs worse than the real prior.

This ordering follows from H2: a univariate filter that reads the endpoint
already beats everything this project has built, and fusion is a harder
estimation problem than selection.

---

## 6. Change log

| Version | Date | Change |
|---|---|---|
| 1.0.0 | 2026-11-10 | First issue. Audit of Omics2PK (63 modules at the time of writing) and MultiOmics; conceptual comparison with MetaPK by reading only; minimal definition; novelty narrowed to one element; six risks with mitigations. |
