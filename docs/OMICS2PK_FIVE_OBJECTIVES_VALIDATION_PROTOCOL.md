# Omics2PK — five-objective validation protocol

**Version:** 1.2.0
**Frozen:** 2026-11-07
**Status:** frozen before execution of the Phase 2–5 runs. Any later change requires a new version number, a date, and a statement of what was changed and why. Thresholds must not be revised after seeing a result.

No threshold, criterion, comparator or primary contrast has been changed since the initial freeze. Version 1.1.0 changes only which sources are admissible, and does so before any run on the affected source — see §12.

## 0. Question and standard of proof

> Can mechanism-guided omics selection improve pharmacokinetic prediction while limiting the risks of high-dimensional molecular data?

The aim is **not** to show that NetPK wins. It is to locate where mechanism-guided selection adds value and where it fails. Negative results, ties and abstentions are results and are retained.

A hypothesis is reported with exactly one of:

| Status | Meaning |
|---|---|
| `SUPPORTED` | A prespecified criterion in this document was actually met. |
| `PARTIALLY_SUPPORTED` | Met on the primary criterion but failed a mandatory robustness criterion, or met in one cohort and not in another. |
| `NOT_SUPPORTED` | The criterion was tested under adequate conditions and failed. |
| `INCONCLUSIVE` | Tested, but the design cannot discriminate (e.g. interval width exceeds the margin). |
| `ABSTAIN` | Not testable: the required endpoint, matching or access does not exist. |

`SUPPORTED` is never awarded on a qualitative impression.

## 1. Admissible data

Governed by `registry/PUBLIC_DATA_EVIDENCE_REGISTRY.json`, validated by `src/omics2pk/public_data_evidence_registry.py`.

A source is admissible for confirmatory analysis when its download gate is one the requester satisfies **unilaterally**. What decides this is the gate, not the tier label:

| Gate | Example | Admissible |
|---|---|---|
| none | file served from a public URL | yes |
| self-signed terms of use | Synapse `SelfSignAccessRequirement` | yes, since v1.1.0 |
| platform certification quiz | Synapse certified-user requirement | yes, in principle; none in use |
| committee decision on an application | Synapse `ManagedACTAccessRequirement`, Estonian Biobank | **no** |
| contacting an author | "data available on reasonable request" | **no** |

The admitted gates involve no third party reviewing, deciding or being contacted, which is what the public-only rule forbids. Accepting data-use terms is a unilateral act that grants access by itself; it is not an authorization request.

**Admissible sources are still reported as gated when they are.** The registry summary field `gated_but_admissible` lists them, and a reader must be able to see which conclusions rest on a source retrievable only after registration. Redistribution restrictions are unaffected: HLC files stay out of this repository and only aggregates are committed.

**Corrected in v1.2.0.** The freeze-time text claimed ARIANA was the only public dataset with a clinical PK endpoint matched to omics *and* pharmacogenetics. That was wrong: Jarmusch releases an individual pharmacogenomic table for 14 participants (CYP2C19, CYP3A4, CYP3A5, CYP2D6) alongside midazolam CL/F, so PK and PGx are matched there too. What is missing in Jarmusch is the *omic* join — the per-donor feature matrix `GNPS_quant.csv` is absent from the authors' repository, and the subject crosswalk is unresolved.

The accurate statement:

> ARIANA n=47 is the only public cohort in which clinical exposure, pharmacogenetics and **named, quantified** biomarkers are matched at the individual level today.

Jarmusch carries clinical PK and PGx matched per individual, with an untargeted omic layer not yet joined per donor. Rotroff (metformin AUC/CLr) carries clinical PK. HLC, EPA, Wegler and the TKI donor panel carry **functional in vitro** endpoints. NCI-60 carries a **cell-line** endpoint. Müller 2026 (n=20, potato meal) carries **biomarker-only** observations with no drug administered.

**These three endpoint classes are never pooled into a common effect estimate.** A functional result is reported as the *functional analogue* of the corresponding hypothesis, never as the hypothesis itself.

### 1.1 Endpoint definitions, frozen

| Cohort | Primary endpoint | Class | Observations per subject |
|---|---|---|---|
| ARIANA n=47 | observed pre-dose metoprolol trough concentration, ng/mL | `pk_clinical` | 1 |
| HLC | CYP substrate turnover activity in liver microsomes, nmol/min/mg | `functional_in_vitro` | 1 per enzyme |
| EPA | targeted CYP activity | `functional_in_vitro` | 1 per enzyme |
| Wegler | in vitro hepatic intrinsic clearance | `functional_in_vitro` | 1 per drug |
| Jarmusch | midazolam apparent clearance CL/F | `pk_clinical` | 1 |
| Rotroff | metformin AUC and renal clearance | `pk_clinical` | 1 |
| NCI-60 | compound GI50/IC50 | `pd_cell_line` | 1 per line |

The ARIANA endpoint is a **trough concentration**. It must not be called AUC, clearance or CL/F.

### 1.2 Prohibited predictors, frozen

- In ARIANA, the **OH-metoprolol/metoprolol metabolic ratio** shares its denominator with the target concentration. It is admissible **only** as a validation target for an activity estimate, never as a predictor in any arm. Enforced by a preflight check that it appears in no design matrix.
- No post-dose molecular measurement may be used as a pre-dose predictor without an explicit temporal qualification.
- A GNPS or library name is not a standard-confirmed structure and is recorded with its identification level.

## 2. Validation machinery, common to every arm

1. **Partitioning.** Repeated subject-disjoint nested cross-validation. Outer 5 folds × 2 repeats; inner 3 folds for every tuning or selection decision. Fold assignment is produced by one shared generator and is **identical across all arms of a contrast**. All records of a subject or donor stay in one outer fold.
2. **Training-only discipline.** Imputation, standardisation, selection, reference construction (including any omic reference level), hyperparameters and locking of mechanistic baselines use training folds only. A test subject contributes to nothing before its own prediction.
3. **Seeds.** Fixed and recorded. Any experiment whose conclusion depends on a noise realisation is replicated over seeds and reported as a distribution, never as a single run.
4. **MAP-Bayesian estimation.** With a single concentration per subject, conditioning the individual random effect on that concentration and then predicting the same concentration is circular; the measured optimism factor is 3.6. **MAP-BE is therefore excluded from every ARIANA arm by construction.** It may appear only in a simulation bench with a known truth, or on a cohort supplying an independent calibration concentration. A simulated MAP-BE result is never placed in the same table as a real observation without an explicit label.
5. **NLME.** A population likelihood is not identifiable from one observation per subject. No ARIANA arm claims an NLME fit.
6. **Uncertainty.** Patient-level paired bootstrap of the loss difference. With n = 47 these intervals are **descriptive**, not confirmatory, and are labelled as such everywhere. Repeated folds do not increase the number of independent subjects; the independent unit count is reported alongside every interval.

### 2.1 Loss functions, frozen

Primary loss: **MAE** in ng/mL (per the objective definition).

Because ARIANA concentrations span roughly 7–572 ng/mL, an arithmetic mean error is dominated by the largest values. Two **mandatory robustness co-criteria** are therefore frozen alongside it:

- 10 % trimmed MAE;
- MAE on the natural-log scale.

Secondary, reported always: RMSE, median absolute error, median log error, geometric bias, P20, P30, calibration slope.

**Reversal rule.** A positive H1 or H2 claim requires the primary criterion to be met **and** no sign reversal on *both* robustness co-criteria. If the primary improves while both robustness criteria reverse, the status is `PARTIALLY_SUPPORTED` at best, and the influence analysis of §2.2 is mandatory in the report.

### 2.2 Influence, frozen

For every contrast that meets its primary criterion, report: the number of subjects improved and worsened on each scale; the share of total net gain contributed by the five most-improved subjects; and the gain recomputed within quartiles of the observed endpoint. A gain concentrated in a few subjects is disclosed as such and does not by itself invalidate the result, provided the direction holds across quartiles.

## 3. Model families, and what is fitted

Arms are compared **at equivalent available information**. Complexity is matched to sample size.

| Arm | Information | Estimator |
|---|---|---|
| M0 | clinical | Ridge; XGBoost |
| M1 | clinical + PGx | Ridge; XGBoost |
| M2 | clinical + omic | Ridge; XGBoost |
| M3 | clinical + PGx + omic | Ridge; XGBoost |
| M4 | mechanistic structure + PGx | fixed-coefficient mechanistic model; reduced PBPK |
| M5 | mechanistic structure + PGx + omic | fixed-coefficient mechanistic model; reduced PBPK |
| P1 | PopPK with **published** covariate coefficients | 1 free parameter |
| P2 | PopPK with covariate coefficients **estimated** on training folds | 3 free parameters |
| P3 | P2 plus the omic | 4 free parameters |
| HY | PopPK structural prediction + ML on the residual | hybrid, new in this phase |

Every arm declares its **number of free parameters fitted on the evaluated cohort** and, for every fixed constant, its value and literature source. Models relying on published parameters are tabulated separately from models estimating parameters on the cohort.

A reduced PBPK transferring Open Systems Pharmacology parameters is labelled a **transferred, simplified model**, never a reproduction of the whole-body model. Evaluating it on the profiles used to build the source model is a **parameter-transfer check, not independent external validation**.

## 4. H1 — incremental value of omics

**Hypothesis.** A mechanistically relevant molecular biomarker adds predictive information beyond clinical covariates and genotype.

**Primary contrast, frozen (one):** ARIANA, **Ridge M3 versus Ridge M1**, primary loss MAE, identical folds.

*Why Ridge and not the mechanistic model.* Ridge is the neutral comparator: it estimates its coefficients from the cohort and so cannot borrow strength from literature-fixed constants, and a regularised linear model is the appropriate complexity at n = 47. Choosing the mechanistic contrast as primary would let a favourable prior do the work. The mechanistic contrast is the designated **key secondary**.

**Key secondary contrast:** mechanistic **M5 versus M4**, same folds, same losses.

**Other secondary:** XGBoost M3 vs M1; PBPK M5 vs M4; P3 vs P1; HY vs P1.

### Criteria

| Outcome | Condition |
|---|---|
**Exploratory success** | ≥ 10 % relative MAE reduction versus the clinical + PGx arm of the **same estimator**, with a descriptive 95 % paired interval for the loss difference excluding zero, and no reversal on both robustness co-criteria.
**Replication success** | Same direction on a biologically independent public dataset with a comparable endpoint, with no new post-hoc selection.
**Failure** | No gain, or degradation, under correctly controlled evaluation.
**Abstention** | No genuine individual endpoint, or no omic↔PK matching.

**Recorded at freeze time, before any run, and re-verified in v1.2.0:** the replication criterion for the CYP2D6/solanidine axis **cannot be met with public data**. The Müller 2024 cocktail study, Wollmann n=839 and Kiiski n=314 are all `ABSTAIN`, and the Estonian Biobank is excluded. The Müller 2026 potato-meal study is open and is now registered, but **it administers no drug**, so it carries no exposure to predict and cannot replicate a PK association. An exploratory success on ARIANA therefore cannot be upgraded to confirmatory in this phase, and must not be presented as such. Replication will be sought only as the **functional analogue** of H1 on HLC/EPA/Wegler, which is explicitly weaker evidence.

**ARIANA conclusions are exploratory by construction, not merely by sample size.** This project has already analysed this cohort repeatedly, with model choices made while results were visible. Freezing the protocol stops further result-guided decisions; it cannot convert a re-analysis into an independent validation. The registry records this as `prior_analysis_exposure` with `max_claimable_status: exploratory`, and no ARIANA result may be reported as confirmatory however favourable it is or however strictly the frozen rules were followed.

## 5. H2 — reliability of mechanism-guided selection

**Hypothesis.** Selection grounded in documented biology improves the relevance and/or robustness of selection relative to purely statistical methods.

Comparators, all on identical folds and at **identical panel dimension**:

1. NetPK prespecified mechanistic selection;
2. SHAP ranking with XGBoost;
3. Elastic Net / LASSO;
4. nested univariate filtering;
5. random panel of the same dimension, replicated;
6. **deliberately wrong or biologically incompatible mechanistic panel**;
7. no selection: all admissible variables.

The mechanistic selection must be derived from knowledge **independent of the PK responses of the evaluated dataset**.

Per variable, the record states: name, molecular identifier, tissue, hypothesised mechanism, level of evidence, decision (retained / excluded / uncertain), reason, and observed predictive contribution after evaluation.

### Criteria

- Final MAE at equal panel dimension; selection frequency over repeated folds; Jaccard stability between selections; false-positive rate on null variables; causal-variable recovery in synthetic known-truth experiments; ability to reject the deliberately wrong panel; abstention rate when evidence is insufficient.
- **Methodological success:** NetPK is non-inferior to the best statistical comparator within a prespecified **5 % relative MAE margin**, *and* shows better stability or better selection specificity.
- **Superiority** requires a favourable paired difference whose uncertainty is compatible with a real improvement. If the interval is wider than the margin, the status is `INCONCLUSIVE`, not equivalence. A non-significant difference is never reported as equivalence.
- **Failure:** NetPK repeatedly retains incorrect pathways, or degrades beyond the margin.

**Mandatory honesty check, frozen.** The selection must be *produced by the stated procedure*, not written into the code after inspecting results. The Random-Walk-with-Restart network ranking is currently a **prioritisation** step whose final retain/exclude decisions are hard-coded. Until an executable rule maps the ranking to the decision, the mechanism is reported as **knowledge-guided selection, not automated discovery**, and H2 cannot be `SUPPORTED` on the autonomy claim. A separate test must verify that the documented rule, applied mechanically, reproduces the retained set.

## 6. H3 — parsimony and minimal information

**Hypothesis.** A small, mechanistically relevant panel retains most of the attainable predictive gain.

Produce, per dataset, performance against number of retained variables for budgets `k = 0, 1, 2, 3, 5, 10, 20, 50, 100, …`, never exceeding the number of qualified variables available. Compare NetPK, SHAP, Elastic Net, repeated random selection and the full panel.

`Gain(k) = MAE(baseline) − MAE(panel k)`.

The reference gain for the percentage denominator is defined **exclusively inside the training folds**. The best test result must never be used to set the target.

Report `k80` and `k90` (smallest panel retaining 80 % / 90 % of the reference gain), the number of variables each method needs, the analytical cost of the panel, and variability across folds and seeds.

- **Success:** NetPK reaches k80/k90 with fewer variables than the comparators, without substantial predictive loss.
- **Abstention:** if no panel improves on the baseline, `k80`/`k90` are undefined. A parsimony score is never manufactured in that case.

## 7. H4 — when omics become harmful

**Hypothesis.** Degradation comes not from the extra information itself but from estimating and exploiting it with finite samples.

Conditions tested: number of candidates; training-set size; variables independent of the target; variables correlated but non-causal; redundancy across transcript / protein / metabolite; signal-to-noise ratio; chemical mis-annotation; wrong biological mechanism; tissue mismatch; sampling-time mismatch.

Real datasets are used wherever they permit the test. **Synthetic-noise experiments use at least 100 independent repetitions**, with original data and targets unmodified.

Report: MAE change versus baseline; probability of selecting at least one decoy; sensitivity and specificity under controlled truth; stability against p/n; the full distribution over repetitions; probability of degradation exceeding 5 % and 10 %; robustness to removing one modality.

- **Success:** reproducible degradation conditions are identified, and whether NetPK limits them better than the comparators is quantified.
- **Frozen prohibition:** no universal `p/n ≈ 2` boundary may be asserted. The deliverable is an **empirical reliability region**, possibly dependent on model, noise and correlation structure. A previous single-realisation sweep is superseded and must be re-run with ≥ 100 repetitions before any boundary statement.

## 8. H5 — biological, analytical and temporal robustness

**Hypothesis.** Biomarkers close to the biological function are more useful when origin, measurement and timing are compatible with the PK mechanism.

Layers studied separately: genotype; total RNA versus functional isoforms; protein; enzyme activity; functional metabolites; inter-omic combinations.

- **ARIANA:** reproducibility of the OH-solanidine and SSDA ratios; differences by CYP2D6 activity score; quantified versus semi-quantitative identification. ARIANA itself documents neither absolute sampling time nor repeated sampling, so those items are not testable *in ARIANA*.
- **Müller 2026 (n=20, standardized potato meal) — the three previously abstained items are now testable.** Registered in v1.2.0; 20 participants, all CYP2D6 genotyped, 9 plasma timepoints over 0–29 h, 180/180 cells present for each of solanidine, 4-OH-solanidine and SSDA. This lifts the `ABSTAIN` on:
  - **sampling-time effect** — the ratio is computed at each of the 9 timepoints and its dependence on sampling instant is quantified;
  - **food effect** — the source records "fasted on potato products for 3 days" for 15 of 20 participants against 5 who were not;
  - **intra-individual temporal stability** — within-subject dispersion of the ratio across the 9 timepoints, reported per participant.

  Four constraints on how this cohort may be used, all recorded before any run:
  1. **It replicates nothing about PK.** No drug is administered. It tests whether the marker behaves reproducibly, never whether it predicts exposure.
  2. **It is not an independent laboratory.** Sarömba, Stingl and Just are authors of both this study and ARIANA. The participants are different people, so donors are not double-counted, but shared method bias would produce agreement just as a true effect would. Any concordance is reported with this stated.
  3. **Limit-of-quantification imputation is already present in the source** and is flagged *only by bold font*: 27 of 180 solanidine cells (15.0%), 15 of 180 for 4-OH-solanidine (8.3%), 9 of 180 for SSDA (5.0%). A value-only read would silently treat imputed floors as measurements, so the bold flags are extracted and imputed cells are reported separately. **Participant 1, the cohort's only poor metaboliser, is below the limit of quantification at all 9 timepoints for both metabolites** — the subject the marker most needs to identify has a fully imputed ratio.
  4. **Activity scores are not comparable with ARIANA's without alignment.** This cohort scores `*41` as 0.25 (post-2021 CPIC); the Rüdesheim 2020 mapping used by the mechanistic module is built on the grid (0, 0.5, 1.25, 1.5, 2, 3) and predates that change. The mapping is not applied across the two cohorts until the convention is explicitly aligned, and the alignment is reported as a decision, not performed silently.

  The source additionally quantifies the **precursor dose** in the meal (α-chaconine, α-solanine, solanidine across 4 aliquots), permitting a precursor-normalised analysis, and includes **stool incubations** showing solanidine is formed only with an active microbial community. The latter is mechanistically double-edged and both directions are tested: microbiota-dependent precursor supply is a reason a measured ratio can beat a genotype estimate, and a reason its reproducibility may be poor.
- **HLC / EPA / Wegler:** RNA → protein → activity → CLint comparison; enzyme-specific variables versus broad panels; redundancy and complementarity across layers. HLC has no matched protein, so the chain cannot be completed within one cohort and must not be assembled across cohorts.
- **Jarmusch:** plasma versus urine; day 1 versus day 9; retention time, annotation level and GNPS networks; blank detection; exclusion of chemically insufficient or contaminated signals.

Criteria: correlation with an independent functional measure; 95 % intervals; intra-individual coefficients of variation where data permit; ranking stability; robustness to matrix and time; identification quality (confirmed / probable / putative / unknown); incremental value after clinical adjustment; consistency with documented enzymology.

- **Success:** a reproducible functional association in an independent dataset with a sufficiently qualified measurement.
- **Frozen limitation:** an independent biological correlation does not by itself demonstrate a PK predictive gain, and is never reported as one.

## 9. Multiplicity

Five objectives × several datasets × several models × several endpoints generate many comparisons.

- **One primary contrast for the whole programme:** ARIANA Ridge M3 vs M1 on MAE (§4). Reported unadjusted.
- The designated key secondary (mechanistic M5 vs M4) is reported unadjusted and labelled secondary.
- All other H1 contrasts across arms and cohorts are **secondary with Holm adjustment** within the H1 family.
- H2–H5 comparisons are **exploratory** and labelled as such; where a family of tests shares an endpoint (for example the ten HLC enzymes), Bonferroni adjustment within that family is applied and the nominal level stated.
- Fold repetition is never counted as additional independent subjects. Every interval states the number of independent biological units.

## 10. Missing data, refusal and sensitivity

- Missing covariates are imputed inside training folds only; a covariate missing for most subjects (for example body weight in HLC) is **not** imputed as though observed, and the arm requiring it is dropped with a note.
- **Refusal criteria:** no deterministic subject-level molecular↔endpoint link; endpoint available only as a published aggregate; access requiring authorisation; omics measured after the prediction landmark without temporal qualification; identification level insufficient for the mechanistic claim.
- **Prespecified sensitivities:** ARIANA formulation (succinate extended-release primary, immediate-release tartrate as sensitivity, with dose equivalence 47.5 mg succinate ≡ 50 mg tartrate and 95 ≡ 100); transferred mechanistic constants (fraction metabolised, omic exponent, genotype/omic weighting, activity-score correction); loss scale; removal of influential subjects; alternative seeds.
- **The formulation is never chosen on the observed MAE of the 47 patients.** The arm ordering must be reported for every formulation scenario.

## 11. Reporting

Final table, per dataset and objective:

| Cohort | Objective | Endpoint class | Baseline | NetPK | Best statistical comparator | ΔMAE | Descriptive 95 % interval | Independent units | Selection | Status |
|---|---|---|---|---|---|---|---|---|---|---|

Favourable and unfavourable results are reported at the same level of detail. Figures must visually distinguish real results, simulations and unverified hypotheses.

A result obtained on a functional or cell-line endpoint is never presented in the same effect column as a clinical PK result.

## 12. Change log

| Version | Date | Change |
|---|---|---|
| 1.0.0 | 2026-11-07 | Initial freeze. Primary contrast, loss hierarchy, reversal rule, multiplicity plan, MAP-BE exclusion, public-only admissibility and the recorded impossibility of CYP2D6-axis replication on public data. |
| 1.1.0 | 2026-11-07 | §1 rewritten: admissibility is decided by the download gate rather than the tier label, and self-signed terms of use are admitted. Requested by the repository owner. Grounds: the HLC gate was inspected against the Synapse public API and is a `SelfSignAccessRequirement` on `syn4499` (id 9604913, `isCertifiedUserRequired=false`, `isValidatedProfileRequired=false`), i.e. no committee and no author contact. Committee-reviewed gates remain excluded. |
| 1.2.0 | 2026-11-07 | **Erratum, four corrections after external review.** (E1) A second Müller study was missed: the audit checked the 2024 cocktail paper and wrongly concluded no public solanidine time-course existed. Müller 2026, doi `10.3389/fphar.2026.1886427`, releases individual data openly; registered and the H5 `ABSTAIN` on sampling time, food effect and intra-individual stability is lifted (§8). (E2) The claim that ARIANA was the only public PK+omics+PGx cohort is withdrawn; Jarmusch has individual PGx for 14 participants (§1). (E3) ARIANA conclusions are stated as an exploratory re-analysis, capped at `exploratory` regardless of result (§4). (E4) "An NLME likelihood is not identifiable from one observation per subject" is corrected to the precise statement: between-subject and residual variability are not *separately* identifiable. No threshold, margin, comparator or primary contrast changed. |

### Validity of the 1.1.0 amendment

A change to admissibility is only legitimate if it cannot have been chosen for the result it produces. Three facts are recorded, and the first is enforced by a test:

1. **No HLC model had been fitted under this protocol when the ruling was made.** No HLC performance figure existed under these rules, so the ruling could not have been selected on one. `tests/test_public_data_evidence_registry.py::test_amendment_records_that_the_ruling_preceded_the_runs` fails if that claim is removed from the registry.
2. **No criterion moved.** Every threshold, margin, comparator and the primary contrast are byte-identical to 1.0.0. The amendment widens the input set; it does not lower a bar.
3. **No hypothesis changed class.** HLC's endpoint is still `functional_in_vitro`. It gains H3/H4/H5 and the functional analogue of H1; it still cannot satisfy H1, and ARIANA remains the only H1-eligible dataset. Becoming retrievable does not turn an in vitro measurement into clinical pharmacokinetics.

The direction of this amendment is worth stating plainly: it **enlarges** the evidence base in which the hypotheses can fail. HLC is the cohort where the existing finding is *unfavourable* to broad panels — a targeted CYP panel beats the broad multi-CYP/POR/transporter panel in 9 of 10 endpoints — so admitting it adds opportunities for refutation, not for confirmation.
