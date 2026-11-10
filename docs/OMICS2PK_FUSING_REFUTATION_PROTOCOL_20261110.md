# A test that can refute Hybrid Fusing, the comparators it needs, the data that exists, and the order of work

**Version 1.0.0 — 2026-11-10.** Evidence type: **protocol, nothing executed**.
Companion to `OMICS2PK_METHOD_AUDIT_AND_HYBRID_FUSING_ASSESSMENT_20261110.md`.
**No frozen protocol or result is modified. No implementation has begun.**

The commissioning note's condition is that implementation starts only after a
test exists that can refute the interest of Fusing. This document is that test.

### Relationship to `docs/HYBRID_FUSING_CONCEPT.md`

That document states the **decision philosophy**; this one states the
**experiment that can refute it**. They were written in parallel and reached
several of the same conclusions independently, which is worth recording because
agreement arrived at twice is stronger evidence than agreement arrived at once:

| Both documents independently require | Where |
|---|---|
| "Fusing" is a working label, **not** an originality or performance claim | concept §1; assessment §5.2 |
| stacking, BMA, mixtures of experts and mechanistic priors are **mandatory comparators** | concept §1; this protocol §2 |
| three distinct operators — qualify, select, integrate-or-abstain | concept §2; assessment §4.2 |
| integration is **not obligatory**; abstention and reversion to a simpler comparator are results | concept §2C; this protocol §3 |
| **"do not reinvent π_m under the Fusing name"** | concept §4; assessment §0 and §5.3(2) |
| conditional selection must be tested **against baseline residuals**, not against the raw outcome | concept §4 (H4); assessment §1.4 and §5.3(2) |
| fusing densities is **not** fusing structural equations | concept §3; this protocol §5.2 (E7) |
| within- **and** between-model variance both reported | concept §3; this protocol §4 |
| `NOT_EVALUABLE` is an explicit permitted result | concept §3.6; this protocol §5.2 |
| no formula of the shape *w = mechanistic × statistical* is authorised | concept §3; assessment §5.3(1) |

**What this protocol adds to the concept document**, and the reason it exists
separately:

1. an **audit** establishing that the hybrid selection rule already exists in
   this repository, has been tested on nine endpoints, and is detectably worse
   than a univariate filter in 21 of 72 cells;
2. a **computable ceiling**, L_cvx, that bounds every member of the fusion
   family and can therefore refute the architecture **before** it is built;
3. **data verdicts**: of seven candidate experiments, one is fully evaluable
   and three are `NOT_EVALUABLE`;
4. a **roadmap whose first implementation step is the fifth**.

Where the two documents differ in emphasis, the concept document governs the
philosophy and this one governs what is run.

---

## 0. Why a refutation test is possible before any code

Any method that combines *M* arms into a prediction — stacking, Bayesian model
averaging, mixture of experts, Hybrid Fusing — produces a convex combination
Σ_m w_m f_m(x_i). The best such combination, with weights chosen **knowing the
answers**, is a computable ceiling. If that ceiling is not detectably better
than the best single arm, **no method in the family can help**, and the
architecture is refuted on that cohort and endpoint without being built.

The ceiling is computable from arms **this project has already fitted**, on
folds it has already frozen. That is the whole point: the gate costs one run
over existing arms and can stop the programme.

MetaPK established the practice and the reason: its oracle reaches MdAPE
11.14 % against 15.68 % for the best single model, so **4.54 points is the
entire recoverable headroom, and any gain must be reported against it rather
than against zero.**

---

## 1. Gate 1 — headroom. The test that can refute the family

### 1.1 Quantities

On a cohort with frozen folds and *M* ≥ 3 arms already fitted, collect
per-individual out-of-fold predictions ŷ_{im} and losses L_{im}:

| | Definition | Bounds |
|---|---|---|
| **L\*** | min_m mean_i L_{im} — the best single arm, chosen **on training folds** | the thing to beat |
| **L_sel** | mean_i min_m L_{im} — per-individual oracle | any method choosing one arm per individual |
| **L_cvx** | min over the simplex of mean_i ℓ(Σ_m w_m ŷ_{im}, y_i) | **stacking, BMA, mixture of experts and Hybrid Fusing** |

L_sel does **not** bound L_cvx: if one arm errs +1 and another −1, the
per-individual oracle's absolute error is 1 and their average's is 0. Both are
reported; **L_cvx is the governing ceiling.**

### 1.2 The verdict rule, fixed here

Paired bootstrap over the resampling unit (donor or patient), 2,000 draws,
Bonferroni over the declared family of cohorts tested.

| Condition | Verdict |
|---|---|
| the interval for **L\* − L_cvx** excludes zero on the favourable side | `HEADROOM_PRESENT` — proceed to Gate 2 |
| it contains zero | **`REFUTED_NO_HEADROOM`** on that cohort and endpoint |
| fewer than 3 admissible arms, or no per-individual predictions recoverable | **`NOT_EVALUABLE`** |

**Two cautions that are part of the rule, not commentary.**

*L_cvx is fitted on the evaluation data*, so it is an optimistic ceiling by
construction — it is an upper bound on what any honest method could reach, which
is exactly what a refutation test needs. It must never be quoted as a
performance.

*A ceiling is not a promise.* `HEADROOM_PRESENT` says a combination could in
principle help; it says nothing about whether weights estimable from training
folds would recover any of it. Gate 2 asks that.

### 1.3 Diagnostics reported with Gate 1

- **Error decorrelation.** The matrix of per-individual error correlations
  between arms, and its median off-diagonal value. Combination gain is driven by
  decorrelation; a median above 0.95 means the arms are one arm.
- **Best-arm shares.** The proportion of individuals for which each arm is
  best. MetaPK's 34.3 / 34.2 / 31.5 % was the signal that no model dominated;
  a 95 / 3 / 2 % split means selection is the answer.
- **Uniform-average loss**, as a control. MetaPK measured that the geometric
  mean of three models was **worse** than the best single one. Averaging without
  weighting is a known way to lose, and it is reported so that any weighted
  result is read against it.

---

## 2. Gate 2 — can weights be learned at all?

Only if Gate 1 clears. Nested, with weights fitted on training folds only and
evaluated out of sample.

| Arm | Weights |
|---|---|
| `S` selection | w = e_{m\*}, m\* chosen on training folds — **the arm to beat** |
| `U` uniform | w_m = 1/M, the control known to lose |
| `K` stacking | w constant, fitted out-of-fold |
| `B` BMA | w_m(i) ∝ each arm's fit to individual *i*'s **own other observations** — available only where individuals have repeated observations |
| `X` mixture of experts | w(z_i) from a gating function of **context covariates**, no mechanistic prior |
| `F` Hybrid Fusing | w(z_i; π), the gating function **with** the mechanistic prior, plus a fallback decided in training folds |

**Primary contrast of the whole programme: `F` vs `K`.** Not `F` vs `U`, and not
`F` vs the baseline. Beating uniform averaging is already known to be easy, and
beating the baseline is not the question. If `F` cannot beat **constant-weight
stacking**, the gating function and the mechanistic prior have contributed
nothing and the honest report is that **stacking is the method**.

**Secondary contrasts, each answering one specific question:**

| Contrast | Question |
|---|---|
| `K` vs `S` | does any fusion beat selection? |
| `X` vs `K` | does context-dependence help, without mechanism? |
| `F` vs `X` | does the **mechanistic prior** add to context-dependence? |
| best of all vs **L_cvx** | how much of the ceiling was recovered? |

The last is the one MetaPK's practice requires, and reporting a gain without it
is reporting a number without its scale.

---

## 3. Gate 3 — the negative controls, declared before any run

Every one is mandatory. Each is designed to make Fusing look good when it is not.

| Control | What it destroys | What must happen |
|---|---|---|
| **degree-permuted prior** | permutes π within degree strata, keeping the degree distribution | `F` must lose its advantage over `X`. If it does not, the prior is measuring **node degree**, i.e. how well studied a gene is — the π_m failure mode transposed |
| **scrambled prior** | permutes π freely | `F` must not beat the real prior |
| **wrong mechanism** | seeds the network from a deliberately incorrect enzyme | `F` must detect it and fall back |
| **permuted context** | permutes z across individuals | `X` and `F` must lose their gating advantage |
| **decoy arms** | adds arms fitted on noise | weights on decoys must go to zero; the loss must not degrade |
| **width-matched comparison** | gives `S` and `K` the same number of free parameters as `F` | a wider arm shrinks harder under an inner-CV penalty and scores better with no information. Measured null medians in this project run to **+0.012 ΔR²** |
| **refitting permutation null** | permutes the added information and **refits everything** | computed for **every** reported contrast, not only the favourable ones. Three of this project's contrasts were withdrawn on exactly this ground |
| **refitting interval** | resamples the unit and refits inside each draw | a non-refitting interval is measured here to be **1.78–2.09×** too narrow for one predictor and **3.76–6.32×** for two |

**And the anti-unfalsifiability rule.** The fallback to selection must be
decided **inside the training folds**, and the committed choice evaluated out of
sample. A method allowed to fall back after seeing the test loss can never lose
and therefore can never be shown to win.

---

## 4. Calibration and width, not a log-score alone

Where the output is a distribution rather than a point — any PK application —
the following are reported together, because MetaPK measured that a log-score
**penalises an absurd component weakly**: uniform weighting had the best
log-score on its five-model library while producing a **546 mg/L** interval, its
0.99 coverage "bought by width".

Required: the primary loss; coverage at 90 %; median interval width; and the
**within / between** variance decomposition, so that disagreement between arms
is reported as disagreement and not as confidence.

---

## 5. The data that actually exists for this

Verdicts under the project's standing rule: individual-level, same-individual
measurement, no cross-cohort fusion, no authorisation required.

### 5.1 Evaluable

| # | Experiment | Cohort | Verdict | Why, and what it limits |
|---|---|---|---|---|
| **E1** | **Fusing-omic, Gates 1–3** | **HLC ≈ 397 donors**, transcriptome, 9 enzyme-activity endpoints, **8 arms already fitted on frozen folds** | **`EVALUABLE`** | the **only** cohort where a gating function is estimable. Endpoint class `functional_in_vitro`: cannot promote a clinical claim. Already analysed, so **exploratory** ceiling |
| **E2** | Fusing of **representations** across layers | **EATRIS-Plus 79 / 55 / 114**, 12 layers, same individuals | **`EVALUABLE`** for a **declared non-pharmacological target only** | no pharmacological endpoint exists. The 12 layers make "fuse or select" a real question; **DIABLO must be a comparator** |
| **E3** | Gate 1 on a **clinical** endpoint | **ARIANA n = 47**: statistical `M1`/`M3`, mechanistic `M4c`/`M4`/`M5c`/`M5` and a hybrid arm, all on identical folds | **`EVALUABLE` for Gate 1 only** | n = 47 cannot estimate a gating function. Tests whether headroom exists at all on a clinical PK endpoint |
| **E4** | Gate 1 on a three-layer cohort | **Lee / Jackson n = 15**, 4 arms | **`EVALUABLE` for Gate 1 only** | at n = 15 the refitting interval already spans [−0.33, +0.74]; expect `REFUTED_NO_HEADROOM` or `NOT_EVALUABLE` |

### 5.2 Not evaluable, with the reason

| # | Experiment | Verdict | Reason |
|---|---|---|---|
| **E5** | **Fusing-PK: combining published PK models on real observations** | **`NOT_EVALUABLE`** | no admissible public source carries **dose, sampling time and concentration** per individual at a size allowing ≥2 PK models to be compared. ARIANA has **one** observation per patient, so no PK model is identifiable; venetoclax has 12 patients with single concentrations. MetaPK's own evidence used MIMIC-IV, which requires **credentialed access** and is therefore excluded here |
| **E6** | **Independent replication of ARIANA** | **`NOT_EVALUABLE`** | searched again on 2026-11-10. **PMC13106024** (venlafaxine + solanidine + CYP2D6, n = 20, Koreans) publishes PK **by phenotype and genotype group**, and its Table S1 is grouped too; the Europe PMC supplement endpoint returns only images. **PMC13609000** (untargeted plasma metabolomics, n = 686, Kyoto J-MICC) publishes **no individual rows**, its endpoint is **genotype-derived** metaboliser status rather than a measured activity, and its six candidate compounds are **MSI level 3** — tentative identifications, not standards-confirmed |
| **E7** | **Mechanistic fusion of PK equations** | **`NOT_EVALUABLE`** | distinct in kind from distribution fusion, as the commissioning note says. Would need composable structural models **and** a cohort on which to compare compositions. Neither exists in admissible public data |

### 5.3 What the data verdicts mean for the programme

**The programme is testable at exactly one place: E1.** Everything else is
either a headroom check on a cohort too small to learn weights, a
non-pharmacological feasibility study, or not evaluable at all.

That is a narrow base, and it should be said plainly rather than spread across
four experiments to look broader. It also sets the realistic ambition: **E1 can
establish whether mechanism-informed fusion beats stacking for predicting an
in-vitro enzyme activity from a transcriptome.** It cannot establish anything
about pharmacokinetics.

### 5.4 The search continues, with these targets

Recorded so the next search is not a repeat of this one:

1. a public cohort with **dose, time and concentration** per patient — the
   binding gap for E5. Candidate channels not yet exhausted: Dryad and OSF
   clinical-pharmacology depositions, journal data-availability statements in
   *CPT:PSP* and *BJCP*, and PK datasets attached to open population-PK
   publications;
2. a second cohort with an **endogenous CYP2D6 marker and a substrate's
   exposure per individual** — for E6;
3. a cohort with **two molecular layers and a drug endpoint** per individual —
   nine such cohorts measured it, three published rows, none has two layers
   *and* a drug;
4. **MetaboLights and PRIDE studies with a drug-administration design**, queried
   by study factor rather than by keyword.

---

## 6. Roadmap, in the order the evidence justifies

No step begins before the previous one reports. Each may end the programme.

| Step | Work | Ends the programme if |
|---|---|---|
| **0** | *(done)* H1–H5 state, method audit, definition, novelty, this protocol | — |
| **1** | **Gate 1 on E1, E3, E4.** One runner over existing frozen arms; no new model, no new architecture | all three return `REFUTED_NO_HEADROOM`. Then report that **selection is the answer** and stop |
| **2** | **Gate 1 diagnostics**: error decorrelation, best-arm shares, uniform control | the median error correlation exceeds 0.95 — the arms are one arm |
| **3** | **Gate 2 on E1**, with `S`, `U`, `K` only. Stacking against selection | `K` does not beat `S`. Then **stacking is refuted too**, and the project's existing selection work stands as the conclusion |
| **4** | **Gate 2 on E1**, adding `X`. Context-dependence without mechanism | `X` does not beat `K`. Then the gating idea is refuted and only the prior remains to test |
| **5** | **Gate 2 on E1**, adding `F`, with **all** Gate 3 controls | `F` does not beat `X`, **or** the degree-permuted prior performs as well as the real one |
| **6** | **E2**: representation fusion on EATRIS-Plus against DIABLO, declared non-pharmacological target | no fusion beats the best single family |
| **7** | Only if steps 3–5 survive: write the method up, and *then* consider whether the PK form is worth pursuing — with E5 still `NOT_EVALUABLE` until new data appears |

**Nothing in steps 1–2 is an implementation of Fusing.** They compute ceilings
and diagnostics from arms that already exist. The first line of Fusing code
belongs in step 5, and only if steps 3 and 4 have both survived.

### 6.1 What will be reported if the programme stops at step 1 or 3

A negative methodological result, which the commissioning note asks to be
treated as valid: *on the largest cohort this project can use, combining the
predictions of mechanistic, hybrid and statistical selectors offers no
recoverable headroom over choosing the best one.* Together with H2's finding
that no selector beats a univariate filter, that would be a coherent and
publishable conclusion — and it would redirect effort from method development
to the data gap, which is where H1 says the binding constraint is.

### 6.2 Duplication explicitly avoided

| Would duplicate | Instead |
|---|---|
| `hybrid_mechanistic_selection.py` | **reuse it unchanged** as the source of π and of the shortlist rule |
| `mechanistic_graph_selection.py` | reuse its degree-stratified permutation FDR and its ordering |
| the H2a/H2d folds, donors, losses and bootstrap | import them, as H2d already imports from H2a |
| MOFA2 / DIABLO | call them a **comparator** for E2, not a thing to rebuild |
| MetaPK's `stacking` / `bma` | **re-derive independently** in this repository. No code is imported, and the EM for a concave weight problem is standard |
| the frozen Axis A, Axis B, ARIANA and H2 runs | untouched; this programme adds arms, never edits a frozen result |

---

## 7. The one-sentence refutation criterion

> **Hybrid Fusing is refuted if, on the only cohort where its weights are
> estimable, either the convex ceiling is not detectably above the best single
> arm, or constant-weight stacking recovers the ceiling as well as a
> mechanism-informed gating function does — and a degree-permuted prior
> performs as well as the real one.**

---

## 8. Change log

| Version | Date | Change |
|---|---|---|
| 1.0.0 | 2026-11-10 | First issue. Three gates, six comparators, eight negative controls, seven data verdicts of which three are `NOT_EVALUABLE`, and a seven-step roadmap whose first implementation step is the fifth. |
