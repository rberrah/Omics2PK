# Axis B2 — independent replication and a plasma-metabolite increment, on the Lee / Jackson CYP3A panel

**Version 1.0.0 — 2026-11-10.** Frozen before execution.
Module: `src/omics2pk/plasma_biomarker_validation.py`.
Runner: `scripts/run_plasma_biomarker_validation.py`.
Source qualification: `docs/OMICS2PK_LEE_JACKSON_CYP3A_PANEL_QUALIFICATION_20261110.md`.

---

## 0. What this is, and what it is not

A **new public cohort** was found on 2026-11-10 that carries, for the same
fifteen hepatocyte donors, a per-donor hepatic protein measurement, a per-donor
plasma metabolite measurement, and a per-donor measured enzyme activity. That
combination does not exist in any other source this project has qualified.

This protocol does two things with it:

1. **Puts the project's one positive functional claim at risk** in donors and a
   laboratory independent of Wegler.
2. **Asks a question the project has never been able to ask on real data:**
   does a *plasma* metabolite add predictive information to a *hepatic* protein
   measured in the same individual?

It is **not** a pharmacokinetic study, it does **not** modify any frozen
analysis, and it introduces **no new predictive architecture** — every fitting
primitive is imported from the frozen `protein_clint_pilot`.

### 0.1 This is not a blind prespecification, and saying so is part of it

The source was read before this document was written. At the time of freezing
the following were already known, from the article's own text and from the
qualification audit:

| Already known at freeze | Value |
|---|---|
| CYP3A4 abundance vs midazolam 1′-hydroxylation, rank, n = 15 | +0.7286 |
| CYP3A4 abundance vs ibrutinib M37, rank, n = 15 | +0.9679 |
| Plasma ratio vs midazolam 1′-hydroxylation, rank, n = 15 | +0.5286 |
| Plasma ratio vs CYP3A4 abundance, rank, n = 15 | +0.4071 |
| Published: M37 and M25 vs plasma ratio, n = 14 | 0.62 and 0.67 |

So the **sign** of every relationship was known before the arms were fixed, and
this protocol cannot claim to have predicted it. What was **not** known, and is
what the execution decides, is:

- whether any of it survives **out-of-sample** leave-one-donor-out prediction,
  which no published analysis of this source has performed;
- whether a gain survives a **permutation null that refits the model**, which
  matters because a wider arm shrinks harder under an inner-CV penalty and so
  scores better with no information at all;
- whether a gain survives an interval that **refits inside every draw**;
- whether the plasma metabolite adds **beyond** the protein, which depends on
  the joint structure and not on either marginal correlation.

A rank correlation of +0.73 and an out-of-sample improvement are different
claims. The first is in the source already; only the second is evidence that a
measurement would be worth making in a new cohort.

---

## 1. The cohort

Fifteen cryopreserved primary human hepatocyte donors (BioIVT), from
Lee J, Fallon JK, Smith PC, Jackson KD, *Clin Transl Sci* 2023;16(2):279–293,
PMCID PMC9926076, CC BY-NC-ND.

| Layer / outcome | Where | Unit |
|---|---|---|
| CYP3A4, CYP3A5 abundance | supplement Table S2 | pmol/mg protein |
| plasma 4β-HC, cholesterol, ratio | **main-text Table 1** | ng/ml, mg/ml, ×10⁴ |
| midazolam 1′-hydroxylation | supplement Table S3 | pmol/10⁶ cells/min |
| ibrutinib M37, M34, M25 | supplement Table S3 | pmol/10⁶ cells/min |
| age, sex, race, BMI | supplement Table S2 | — |

**One cohort, not two.** Lee et al., *Clin Transl Sci* 2024 (PMC11645475)
reports venetoclax metabolism in the **same** fifteen donors and the same
twenty microsome lots. Two drugs measured in the same donors are one cohort
with two substrates. The venetoclax paper is therefore a source of comparator
numbers and never a second cohort.

**Twenty microsome lots exist and cannot be used.** Supplement Table S1 gives
abundance for twenty single-donor microsome lots with **no activity** for them.
The microsomal arm that Wegler could test is not testable here, which matters
because Wegler's microsomal result was the negative one.

---

## 2. The circularity check, done before anything was run

The Wegler pilot found two of three candidate endpoints circular: equation 8
divides the clearance by the probed CYP and multiplies back by the homogenate
CYP concentration, and equation 9 adds K<sub>puu</sub>. Both would have made a
CYP-abundance predictor predict itself.

**This source has no such quantity.** The outcome is a formation rate
normalised **per million cells**; the predictor is an abundance **per mg of
protein**. The outcome is divided by no CYP quantity, so no arm here is
circular.

That difference has a cost, and it is recorded rather than hidden: the two are
expressed per **different denominators**, so each donor's protein-per-cell
ratio enters the relationship as an uncontrolled term. It weakens the
comparison. It does not make it circular.

`FORBIDDEN_ENDPOINTS` names all four predictor columns, and
`check_imported_machinery` asserts that the guard still refuses each of them,
so a future hand-edit cannot quietly make an endpoint out of a predictor.

---

## 3. The four arms

| Arm | Predictors |
|---|---|
| `V0_mean` | none — the training-fold mean |
| `V1_hepatic_protein` | CYP3A4 abundance |
| `V2_plasma_metabolite` | plasma 4β-HC / cholesterol ratio |
| `V3_protein_plus_plasma` | both |

Two predictors at most, at fifteen donors. No component extraction, no
selection among many candidates, no interaction term.

**Declared sensitivity, not a second primary:** `V1s_cyp3a_total` and
`V3s_cyp3a_total_plus_plasma` replace CYP3A4 with CYP3A4 + CYP3A5, because
CYP3A5 also 1′-hydroxylates midazolam. It carries no multiplicity adjustment of
its own and may corroborate the primary reading, never replace it.

---

## 4. The three contrasts

| Name | Contrast | Question |
|---|---|---|
| `replication_V1_vs_V0` | protein vs mean | does the enzyme's own abundance predict its probe reaction out of sample, independently of Wegler? |
| `incremental_value_V3_vs_V1` | **primary** — both vs protein | does a plasma metabolite add to the hepatic protein in the same donors? |
| `substitution_V2_vs_V0` | plasma vs mean | does the plasma metabolite alone beat no information? |

`incremental_value_V3_vs_V1` is the **primary** contrast because it is the one
no existing public analysis answers, and because it is the structural question
the funding case rests on: whether a second layer adds beyond the first.

The three are **not** interchangeable. A positive `substitution` with a null
`incremental_value` would mean the plasma metabolite carries information that
the protein already carries — redundancy, not complementarity. That reading is
declared here so it cannot be presented later as a complementarity result.

### 4.1 Endpoints

Primary endpoint **midazolam 1′-hydroxylation**, because it is the same probe
reaction as the Wegler contrast being replicated. Secondary endpoint
**ibrutinib M37**. Both on log₁₀, errors restored as a geometric fold error and
in pmol/10⁶ cells/min.

The endpoint is an **activity rate at a single substrate concentration**. It is
not V<sub>max</sub>/K<sub>m</sub>, not a fitted intrinsic clearance, and not on
the same scale as Wegler's equation 7. The replication is of a **relationship**
and a **direction**, never of a number.

### 4.2 Multiplicity

Three contrasts × two endpoints = **six**. Bonferroni at that family size, by
`declared_family_size = 6`, which the runner asserts equals the design it
actually executes. The sensitivity arms sit outside the family and are labelled
as such.

---

## 5. Validation and preprocessing

- **Leave-one-donor-out.** Fifteen folds, one donor each.
- The ridge penalty is chosen by an **inner leave-one-out among the training
  donors only**, over the frozen grid `(0.001, 0.01, 0.1, 1.0, 10.0, 100.0)`.
- Standardisation constants come from the **training donors only**, in every
  fold and inside every bootstrap draw and permutation replicate.
- The held-out donor's endpoint is never seen by any step that produces its
  prediction.
- No imputation: a missing or non-positive endpoint raises rather than being
  dropped. The qualification audit confirmed all fifteen donors are complete on
  all four fields, so no rule needs to fire.

---

## 6. Controls

### 6.1 The permutation null, for every contrast

For each contrast, the **added** block's donor labels are permuted, the
comparator's predictors and the endpoint stay paired, and **everything is
refitted** — standardisation, penalty, coefficients — inside each of 200
replicates.

The null is computed for **every** reported contrast, not only the ones that
look significant. The Axis A run established why: a wider arm shrinks harder
under an inner-CV penalty and therefore scores better with no information, so
the null is **not** centred on zero and a bootstrap-significant gain lying
inside a wholly positive null is an artefact. Three Axis A contrasts were
withdrawn on exactly that ground.

Replicate values are reported as **deciles**, never as a 200-number list.

### 6.2 Two intervals, and which one governs

| Interval | Refits inside each draw | What it prices |
|---|---|---|
| `stored_error_bootstrap` | **no** | the choice of donors used for scoring |
| `refitting_out_of_bag_bootstrap` | **yes** | the choice of donors learned from **and** scored |

The first is the Wegler pilot's estimator, kept so the two cohorts are read on
the same instrument. The second answers the objection raised against the Wegler
pilot: that resampling stored out-of-fold errors without refitting does not
capture the uncertainty of learning on fifteen individuals.

The refitting estimator draws fifteen donors with replacement, refits both arms
on the draw — penalty selection included — and scores **only the donors the
draw left out**, so no donor is ever both learned from and scored. About a
third of the cohort is scored per draw, which makes the interval wide. That
width is the honest answer at fifteen donors and not a defect of the estimator.

**The refitting interval governs.** Where the two disagree, the stored-error
interval is the optimistic reading and must not be quoted alone. The artefact
records `the_two_intervals_disagree` for every contrast.

### 6.3 The verdict rule, fixed before execution

`SUPPORTED` requires **all three**:

1. the **refitting** interval excludes zero;
2. its point estimate is an improvement (negative, i.e. less error);
3. the permutation rank p-value is ≤ 0.05.

Anything else is `NOT_SUPPORTED`. The stored-error interval alone is never
sufficient.

---

## 7. What a result here can and cannot support

| Reading | Permitted? |
|---|---|
| the enzyme→activity relationship reproduces in a second independent cohort | **yes**, if `replication_V1_vs_V0` is SUPPORTED |
| a plasma metabolite adds to a hepatic protein for predicting hepatic activity | **yes**, if `incremental_value_V3_vs_V1` is SUPPORTED |
| a plasma metabolite could replace a liver measurement | **no** — would need an equivalence design, which fifteen donors cannot carry |
| anything about predicting a concentration, an AUC or a clearance in a patient | **no** — there is no dose and no drug concentration in this source |
| the plasma biomarker works in living patients | **no** — this is procurement plasma from deceased donors |
| an absence of gain shows the plasma metabolite is uninformative | **no** — at n = 15 an interval containing zero means no difference was **detected** |

**Absence of a difference is not equivalence.** At fifteen donors the
deliverable is the estimate and its width, not a significance verdict.

---

## 8. Relation to the frozen Wegler result

Nothing in the Wegler pilot is modified. Its result stands as published and
this cohort is reported beside it, never pooled with it:

- different donors, different laboratory, different year, different
  quantification;
- a **fitted intrinsic clearance** there, a **single-concentration formation
  rate** here;
- both are hepatocyte systems and midazolam is the probe in both.

Two cohorts agreeing in direction is a replication of the direction. It is not
one estimate measured twice, and the two numbers are never averaged.

---

## 9. Freeze table

| Item | Value |
|---|---|
| seed | 20261110 |
| validation | leave-one-donor-out, 15 folds |
| inner penalty selection | leave-one-out among training donors |
| ridge grid | 0.001, 0.01, 0.1, 1.0, 10.0, 100.0 |
| primary endpoint | `midazolam_1_hydroxylation` |
| secondary endpoint | `ibrutinib_m37` |
| primary arm | `V3_protein_plus_plasma` |
| primary comparator | `V1_hepatic_protein` |
| replication contrast | `V1_hepatic_protein` vs `V0_mean` |
| substitution contrast | `V2_plasma_metabolite` vs `V0_mean` |
| permutation replicates | 200, model refitted under the null |
| stored-error bootstrap draws | 2000 |
| refitting out-of-bag draws | 2000 |
| declared family size | 6 |
| target scale | log₁₀ |
| forbidden endpoints | all four predictor columns |
| status ceiling | `exploratory_independent_replication` |

---

## 10. Paths

- module `src/omics2pk/plasma_biomarker_validation.py`
- runner `scripts/run_plasma_biomarker_validation.py`
- artefact `reports/data/plasma_biomarker_validation.json`
- per-donor table: **outside this repository**, since the source is
  CC BY-NC-ND and this project commits no individual-level rows.

---

## 11. Change log

| Version | Date | Change |
|---|---|---|
| 1.0.0 | 2026-11-10 | Frozen before execution. §0.1 records that the source's own correlations were known at freeze and that this is therefore not a blind prespecification. |
