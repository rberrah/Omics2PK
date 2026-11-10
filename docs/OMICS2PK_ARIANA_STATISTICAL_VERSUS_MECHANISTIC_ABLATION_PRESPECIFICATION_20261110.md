# ARIANA prespecification — the three-way statistical ablation, and statistical against mechanistic

**Version 1.0.0 — frozen 2026-11-10.** Objective: H1 family, secondary.
Status ceiling: `exploratory`. Cohort: `ariana_metoprolol_47`.

This protocol fixes a **contrast family and its multiplicity adjustment** over
arms that already exist. It fits **no new model**, introduces **no new arm**,
changes **no frozen fold**, and touches **neither** the frozen primary contrast
(`M3_clinical_pgx_omic` versus `M1_clinical_pgx`) **nor** the frozen key
secondary (`M5_mechanistic_pgx_omic` versus `M4_mechanistic_pgx`).

---

## 1. Honesty about what is and is not prespecified

The commissioning note asks for a statistical ablation of genotype alone,
metabolomics alone, and both, on the same 47 patients, the same folds and the
same endpoint. **Those arms are already implemented and already computed.**
`scripts/run_h1_ariana_incremental_omic_value.py` carries
`M1a_clinical_genotype`, `M1b_clinical_comedication`, `M2_clinical_omic` and
`M3a_clinical_genotype_omic`, added as declared secondaries on 2026-11-07, and
their contrasts are in the committed `reports/data/h1_ariana_trough.json`.

Two consequences, stated here rather than discovered later:

- **Family A below is not a blind prespecification.** Its contrasts were
  visible in a committed report before this document was written. What this
  document adds is the declared family, the multiplicity adjustment the H1
  runner's own note requires before any claim, and a written decision rule.
  Family A may therefore be reported as *an adjusted reading of an existing
  computation*, never as a newly prespecified test.
- **Family B is new**, in that no document has reported a paired contrast
  between a statistical arm and a mechanistic arm. But the arm-level mean
  absolute errors of both families were already published in
  `reports/data/h1_ariana_trough.json` and
  `reports/data/h1_ariana_mechanistic_hybrid.json`, so the *direction* of every
  Family B contrast was foreseeable when this document was written. Family B is
  therefore prespecified as to its **family membership, estimator, adjustment
  and decision rule**, and not as to its direction.

Running a second, differently folded benchmark on these 47 patients was
considered and **rejected**: the commissioning note asks for the existing
folds, and a fresh ablation would have been a fifth pass over a cohort that has
already been analysed four times, with a different partition, answering the
same question less comparably. Not multiplying models is itself an instruction
in the note.

---

## 2. Why all three arms keep the clinical baseline

The note asks for a *fair* comparison between genotype alone, metabolomics
alone, and both. Fair means the same baseline, and the baseline must contain
the **daily dose**: the endpoint is an observed pre-dose metoprolol
concentration, and a concentration model without dose is not a weaker model but
a wrong one. All three arms therefore sit on the frozen clinical block
(`Age`, `male sex`, `weight`, `BMI`, `Daily Dose Metoprolol`), and "alone"
means "this block and nothing else of its kind".

A literally predictor-only arm — genotype with no dose — is **not** computed,
because it would not be a comparator for anything in this protocol.

This also matches the project's standing requirement: the informative contrast
is baseline-plus-feature against baseline, never feature against baseline.

---

## 3. The endpoint and what it is not

`Trough concentration metoprolol`, an **observed pre-dose concentration** in
ng/mL. It is **not** an AUC, **not** a clearance and **not** a steady-state
exposure metric, and it is never called one. The prohibited predictor
`lnMR OH-metoprolol/metoprolol` appears in no arm. There is one observation per
patient, so no MAP-Bayesian estimate conditions on the quantity it predicts.

---

## 4. The declared contrast family

Eleven contrasts. The family is fixed here, before the runner executes, and the
adjustment uses the family size **11**.

### Family A — the three-way statistical ablation (6)

| Key | New arm | Baseline | Question |
|---|---|---|---|
| `A1_genotype_vs_clinical` | `M1a_clinical_genotype` | `M0_clinical` | does genotype add to clinical? |
| `A2_comedication_vs_clinical` | `M1b_clinical_comedication` | `M0_clinical` | does the inhibitor flag add? |
| `A3_metabolomics_vs_clinical` | `M2_clinical_omic` | `M0_clinical` | does metabolomics add? |
| `A4_both_vs_clinical` | `M3a_clinical_genotype_omic` | `M0_clinical` | does the pair add? |
| **`A5_pgx_on_top_of_metabolomics`** | `M3a_clinical_genotype_omic` | `M2_clinical_omic` | **does PGx add information the metabolomics does not already carry?** |
| `A6_metabolomics_on_top_of_genotype` | `M3a_clinical_genotype_omic` | `M1a_clinical_genotype` | the reverse direction |

`A5` is the contrast the commissioning note names. It is reported with its
adjusted interval and with the direction of its point estimate stated
explicitly, in whichever direction that falls.

### Family B — statistical against mechanistic, same folds (5)

| Key | New arm | Baseline | Panel |
|---|---|---|---|
| `B1_statistical_vs_mechanistic_omic_matched` | `M2_clinical_omic` | `M5c_mechanistic_omic_only` | `oh_and_ssda` (matched: two ratios each) |
| `B2_statistical_vs_mechanistic_with_genotype_matched` | `M3a_clinical_genotype_omic` | `M5_mechanistic_pgx_omic` | `oh_and_ssda` |
| `B3_statistical_vs_mechanistic_genotype_only` | `M1a_clinical_genotype` | `M4_mechanistic_pgx` | none |
| `B4_the_two_baselines` | `M0_clinical` | `M4c_backbone_only` | none |
| `B5_the_unmatched_panel` | `M2_clinical_omic` | `M5c_mechanistic_omic_only` | `oh_only` (one ratio against two) |

`oh_and_ssda` is **already a declared panel** in
`src/omics2pk/ariana_mechanistic_arms.py` and already a reported sensitivity
cell of the mechanistic run, so using it here introduces nothing new. `B5`
exists so that the cost of the panel choice is visible rather than buried: the
mechanistic primary panel is `oh_only`, the statistical arms carry both ratios,
and comparing them directly would otherwise compare representations and feature
counts at the same time.

### 4.1 What Family B does **not** measure

A mechanistic arm here fits **one** free parameter on these 47 patients and
borrows every other constant from the literature. A statistical arm fits a
ridge penalty and six to eight coefficients on the same 47 patients. Family B
therefore compares *"transferred structure plus one free parameter"* against
*"estimation from 47 patients"*. **It does not isolate the representation.**
Any sentence that reads a Family B result as "the mechanistic representation is
better" is forbidden by §7; the permitted reading names both ingredients.

---

## 5. Estimator, folds, losses

Nothing here is new. Everything is imported from the frozen runners rather
than reimplemented, and the runner verifies the imported constants before it
starts.

| Item | Value |
|---|---|
| Patients | 47, one row each |
| Folds | `H1.make_folds`, 5 outer × 2 repeats, 3 inner, seed `20261107` |
| Fold identity | the same fold list for every arm in both families |
| Model | the frozen `ridge` for statistical arms; the frozen mechanistic arms as published |
| Primary loss | mean absolute error on the observed concentration |
| Co-criteria | log MAE and 10 %-trimmed MAE, with the frozen reversal rule |
| Bootstrap | 2,000 paired draws over **patients**; a fold repeat is never a second patient |
| Multiplicity correction | Bonferroni across the family of **11**, applied to the bootstrap percentiles |
| Exploratory threshold | ≥ 10 % relative MAE improvement **and** the adjusted interval excluding zero |

**Bonferroni, not Holm.** The H1 runner's note asks for Holm. The repository's
adjusted-interval machinery (`bonferroni_interval`) is Bonferroni, and
Bonferroni is uniformly more conservative: an interval that excludes zero after
Bonferroni would also exclude it after Holm. Using the stricter available
adjustment is therefore sound, and it reuses tested code instead of adding an
untested path. The choice is recorded here rather than explained afterwards.

---

## 6. Decision rules

1. A contrast carries a claim only if its **adjusted** interval excludes zero
   **and** the relative improvement reaches 10 %.
2. An unadjusted interval that excludes zero while the adjusted one does not is
   reported as **inconclusive within the family**, never as a positive finding.
3. An interval that includes zero means no effect was **detected**, never that
   there is none.
4. The whole protocol is capped at `exploratory`. ARIANA has been analysed
   repeatedly by this project; nothing here is an independent validation, a new
   cohort, or evidence transferable to another drug.
5. **No result from this protocol may be presented as a clinical-PK proof.**
   It is a comparison of representations on one observed trough concentration
   in 47 patients from one study.

---

## 7. What may not be done after seeing the results

- No new arm, panel, baseline, fold, loss, threshold or family member.
- No re-reading of Family B as a statement about representation alone (§4.1).
- No promotion of an unadjusted interval to a claim.
- No restatement of the frozen primary contrast or key secondary.
- No sentence calling the endpoint an AUC, a clearance or an exposure metric.
- A defect found after this freeze is recorded in §9, not repaired silently.

---

## 8. Freeze

| Item | Value |
|---|---|
| Family size | 11 |
| Fold seed | `20261107` |
| Outer folds × repeats × inner | 5 × 2 × 3 |
| Bootstrap draws | 2,000 |
| Adjustment | Bonferroni |
| Relative MAE threshold | 0.10 |
| Matched mechanistic panel | `oh_and_ssda` |
| Unmatched mechanistic panel | `oh_only` |
| Note's named contrast | `A5_pgx_on_top_of_metabolomics` |
| Status ceiling | `exploratory` |

`scripts/run_ariana_statistical_versus_mechanistic_ablation.py` reads this
table back and refuses to run if a value differs.
`tests/test_ariana_statistical_ablation.py` fails if this document and the
runner disagree.

---

## 9. Change log

| Version | Date | Change |
|---|---|---|
| 1.0.0 | 2026-11-10 | First freeze. Records in advance that Family A is an adjusted reading of an already-visible computation and not a blind test (§1); that a second differently folded benchmark was considered and rejected (§1); why all three ablation arms keep the clinical and dose baseline (§2); the matched-panel requirement and the `oh_only`-versus-`oh_and_ssda` asymmetry (§4); that Family B cannot isolate the representation because the arms differ in how many parameters they fit on this cohort (§4.1); and that Bonferroni replaces the note's Holm because it is the stricter adjustment and the tested one (§5). |
