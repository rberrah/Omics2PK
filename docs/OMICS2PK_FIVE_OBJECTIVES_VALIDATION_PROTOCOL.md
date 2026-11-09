# Omics2PK — five-objective validation protocol

**Version:** 1.4.6
**Frozen:** 2026-11-08
**Status:** frozen before execution of the Phase 2–5 runs. Any later change requires a new version number, a date, and a statement of what was changed and why. Thresholds must not be revised after seeing a result.

No threshold, criterion, comparator, margin or primary contrast has been changed since the initial freeze. Version 1.1.0 changes only which sources are admissible. Version 1.4.0 rewrites §5 so that the mechanistic selection is computed rather than written down, and records a negative finding about what a network ranking can decide; both were done before any enzyme-activity value was read — see §12.

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

**Further corrected in v1.3.0.** Jarmusch's omic layer is *also* joined per donor: the GNPS processed output `qiime2_table.qza` (6,349 MS1 features × 1,198 sample IDs) filtered to plasma / day 1 / time 0 min gives 13 matched subjects, 12 with a day-9 occasion. So **H1 has two public clinical-PK tests, not one**, and the second has already produced a **negative** result: untargeted metabolomics was worse than measured CYP3A4/CYP3A5 variants on both days (day 1 ΔMAE +2.545, descriptive 95 % [+0.4421, +5.9959], excluding zero on the harmful side; day 9 +0.148, [−0.4273, +0.6899]). It stays `ACQUIRED` rather than `QUALIFIED` because its feature universe was built across all samples before fold splitting, which §2 forbids — a deviation that favours the omic arm, so the negative result is conservative. Jarmusch's ions are mostly unidentified, which is why the uniqueness claim above, resting on *named and quantified* biomarkers, is unaffected. Rotroff (metformin AUC/CLr) carries clinical PK. HLC, EPA, Wegler and the TKI donor panel carry **functional in vitro** endpoints. NCI-60 carries a **cell-line** endpoint. Müller 2026 (n=20, potato meal) carries **biomarker-only** observations with no drug administered.

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

**Status of this section.** Rewritten in v1.4.0, **before any enzyme-activity value was read**, to replace a hard-coded retained set with an executable chain and to record a negative methodological finding established from the network alone. Every comparator, criterion and margin frozen in v1.0.0 is carried over unchanged; what changes is that the NetPK arm is now computed rather than written down, and that the *size* of its panel is no longer claimed to come from the network. The chronology matters and is verifiable: the panels, the sample sizes and the finding below were obtained from the network, the pathway annotation and the array's probe annotation, none of which contains an endpoint.

**ARIANA cannot carry H2.** It has exactly two candidate omic variables. Choosing two of two is not selection, and no selector can distinguish itself on a panel that size, so running this comparator set there would produce a near-tie by construction. H2 runs on the wide panels.

### 5.1 The executable chain, frozen

| Step | Frozen specification |
|---|---|
| Network | STRING v12.0 human, `combined_score ≥ 700`. Identifiers collapsed onto preferred gene symbols keeping the strongest edge. Yields **16,201 nodes and 236,930 undirected edges**, recorded at qualification. |
| Propagation | Random walk with restart, restart probability **α = 0.30**, column-stochastic transition, L1 tolerance 1e-11, at most 400 power iterations. |
| Seeds, variant **S1** | The endpoint's own column name, resolved against the network: drop the dataset's documented suffix, then drop trailing single uppercase characters until the remainder is a network node. No gene list is consulted and no family is special-cased. |
| Seeds, variant **S2** | Every member of the Reactome pathways matching a frozen name pattern, equally weighted. **Identical for all ten endpoints**, so no endpoint-specific information enters at all. |
| Seeds, variant **S3** | Genes carrying a curated `metabolic processing` interaction with the **documented probe substrate**, weighted by the number of distinct supporting references. The enzyme is never named. |
| Output | An ordering of every non-seed gene. **The network supplies the order, not the size.** |
| Gene → probe | A gene's candidates are all platform probes carrying its symbol; where several exist, the one with the highest mean expression **in the training fold**. |

**S3 is available for two endpoints only, and this is a provenance limit, not a choice.** The probe substrates of eight of the ten assays are not retrievable: the source archive's readme names only the two used for the duplicated enzyme, the 2008 publication does not describe the activity assays at all, and the 2010 publication's substrate list sits in a supplement absent from the open repository, with the publisher's supplement path returning 404. Naming the other eight from memory would be an unverifiable input of exactly the kind this protocol forbids. S3 is therefore reported as an **illustration of the full chain**, never as a test.

### 5.2 Why the panel size is not delegated to the network

**Finding, established before any endpoint was read, and reported as a result rather than worked around.** A network propagation score cannot be given a calibrated significance threshold.

- An empirical p-value from `B` resamples cannot fall below `1/(B+1)`, while a Benjamini–Hochberg threshold over sixteen thousand candidate genes is of order 1e-6. The first version of the rule could therefore never retain anything, whatever the data.
- Replacing the p-value with a permutation false-discovery estimator on **node counts** removes that resolution floor, but not the deeper obstacle: every null that preserves the degree sequence also preserves, approximately, *how many* nodes reach a given proximity to a seed. Relocating a seed to a comparably connected node gives its new neighbours about the scores the real seed's neighbours have, so the count ratio sits near one wherever the walk is informative.
- Measured on the real network: the threshold rule retains **one or two genes for all ten endpoints**, and retains **the same number when deliberately given the wrong enzyme**. With many seeds it does retain — 66 to 205 genes for the pathway variants — so the diagnostic has power for *module coherence* and essentially none for *single-seed neighbour identification*.

The conclusion is not that the parameters need tuning. It is that network proximity carries information about **which** genes are partners and almost none about whether being a partner is unusual. The panel size is therefore chosen inside **training folds only**, exactly as a data-driven competitor chooses its own. Every arm then uses the endpoint in the same way, and the only difference between arms is where the ordering came from — which is the comparison H2 actually requires. The threshold rule is kept and reported as the diagnostic of this fact.

**Withdrawn in advance.** Any claim that the network autonomously decides how many variables to retain. H2 can be evaluated on the autonomy of the **ordering**; the autonomy of the **size** is not achievable by this mechanism and is not claimed.

### 5.3 Four experiments, kept separate

**H2a — does the mechanistic ordering select better?** Same donors, same folds, same predictor, same number of variables for every arm. Comparators, as frozen in v1.0.0:

1. NetPK: the executable mechanistic ranking of §5.1;
2. SHAP ranking with XGBoost;
3. Elastic Net / LASSO by coefficient magnitude;
4. nested univariate filtering;
5. random panel of the same dimension, **≥ 100 replications**;
6. **deliberately wrong mechanistic panel**: the identical chain seeded with an enzyme that does not perform the measured reaction, by a rotation fixed before running;
7. no selection: all admissible probes;
8. the target transcript alone, which is the information the endpoint's label gives away for free and is the arm NetPK must beat to claim it finds *complementary* mechanisms.

Panel size `k` is swept over the frozen grid `k = 1, 2, 3, 5, 10, 20, 50, 100` and, separately, chosen in the inner folds; both are reported. All selection, preprocessing and tuning happen in training folds; no arm sees a test fold.

**H2b — can the chain identify the mechanism from the drug, rather than being told the enzyme?** Introduced in v1.4.3 at the request of an external review, which drew the distinction the earlier text had blurred.

H2a gives the chain the enzyme that performs the measured reaction. That is **contextualised mechanistic selection, not autonomous discovery of the mechanism**, and must be called so throughout: the endpoint's own label supplies the enzyme, and the only thing the ordering can add is which *further* genes to retain. H2b is the harder and more interesting question — given a drug, can the chain find the enzymes, transporters and pathways that matter, and then select the omic variables attached to them?

The instrument already exists: the substrate-seeded variant of §5.1, which withholds the enzyme name and seeds from the documented probe substrate. Its reach on this cohort is limited to **two endpoints** by the provenance limit recorded in §5.1, so H2b is **not** evaluated as prediction there. It is evaluated instead as **mechanism recovery**, which needs no pharmacokinetic measurement at all and is therefore not limited by cohort availability:

- for a set of drugs with an independently documented primary metabolising enzyme, does the chain rank that enzyme first from the drug name alone;
- how often does it rank an enzyme that no guideline associates with the drug (false identification);
- does it abstain when the drug is absent from the interaction source, rather than returning its nearest guess.

**Ground truth must not come from the chain's own input.** The chain reads a curated chemical-gene interaction resource, so taking the reference drug-enzyme pairs from that same resource would be circular. The reference is instead an independently curated clinical consensus list. Even so, both are distilled from overlapping literature, and the benchmark therefore measures whether a public interaction database plus a public network **reproduces a clinical-guideline consensus** — retrieval reliability, false-identification rate and abstention — and **not** whether anything new is discovered. That limit is stated here so no later text can read the result as discovery.

**H2d — does hybrid selection keep the parsimony without the vulnerability?** Added in v1.4.4. Knowledge proposes a shortlist, training data decides what to keep, demote or abstain from. The rule, its criteria and its failure modes are frozen in a separate document, `OMICS2PK_H2D_HYBRID_SELECTION_PRESPECIFICATION_20261108.md`, and **not** in this one, for a reason that must travel with the result: the hypothesis was generated by the H2a destroyed-signal diagnostic on one endpoint, so it is hypothesis-generating there and confirmatory only on the other nine. That document names the generating endpoint, caps its result at `exploratory`, and states in advance the failure mode most likely to occur — that the hybrid collapses onto the univariate selector, which would mean the knowledge contributed nothing beyond shortlisting.

An external review proposed this strategy under the label "H2c". That label had already been taken, two revisions earlier, by the integration question below, so the hybrid strategy is **H2d**. The collision is recorded rather than resolved silently.

**H2c — does mechanistic integration then help?** Renamed from H2b in v1.4.3 to free the label for the question above; the content is unchanged. Identical NetPK variables, statistical integration against parameterisation as an enzyme-activity score feeding a structural model. Separating selection from integration isolates the value of the ordering from the value of pharmacological coefficients borrowed from the literature. H2c is not evaluated until H2a is complete.

**What the decisive result is, and is not.** Not that NetPK obtains the lowest MAE. The questions that bear on the method's worth are whether it can work with few variables, whether it resists a signal that is not there, and whether it recognises the cases where its own knowledge is insufficient. A favourable MAE with a high false-identification rate or no abstention capability is a weaker result than a modest MAE with both.

### 5.4 Validation design, frozen

- **Primary:** donor-disjoint repeated nested cross-validation, 5 outer folds × 2 repeats, 3 inner folds, identical folds across every arm and model, seed `20261108`. **Analysis unit, made explicit in v1.4.3:** out-of-fold results are stored per repeat, and a donor's absolute error is averaged across repeats before any loss or bootstrap. The donor stays the independent unit and a repeat is never a second donor. The first implementation wrote every fold into one vector per arm, so the second repeat overwrote the first and the reported losses came from one pass; no result had been published from it.
- **Secondary, transportability:** leave-one-collection-centre-out, 3 folds. The source readme states the cohort came from three independent tissue collection centres; the donor identifier's leading digit separates the 466 expression donors into exactly three groups of **245, 184 and 37**. This is an **inference from the identifier, not a documented field**, and is labelled as such wherever it is used. The 37-donor group is used as a held-out set only, never as a training set alone.
- **Endpoint class:** `functional_in_vitro` — measured microsomal activity in nmol/min/mg protein. Never pooled with `pk_clinical` and always reported in a separate column. The sample size is **387 to 398 donors per endpoint** after the donor join, which is verified rather than assumed.
- **Baseline covariates: age and sex only.** Corrected in v1.4.1, before any enzyme-activity value was read, for a reason independent of any endpoint: weight is recorded for **144 of 466** donors and body mass index for **124 of 466**, while age is complete and sex is missing for one donor. Keeping weight and body mass index would discard about seventy per cent of the cohort or impute the majority of their values. Both are therefore excluded and the exclusion is recorded here rather than discovered later. The one missing sex is imputed with the training-fold majority.
- **Missing values, frozen.** Expression is supplied as log ratios already clipped to [−2, 2], with about 0.3 % of cells absent. A probe missing in more than 20 % of a training fold is dropped in that fold; otherwise a missing cell takes the training-fold mean of its probe. No imputation ever reads a test fold or an endpoint.
- **Feature universe, frozen.** One feature per gene symbol: the probe with the highest mean expression **in the training fold**. Probes with no symbol are excluded. This makes the panel dimension mean the same thing in every arm, since the mechanistic ranking orders genes while the statistical comparators would otherwise order probes and could spend their budget on several probes of one gene.
- **Loss hierarchy:** as §2.1, unchanged. MAE primary, with 10 % one-sided trimmed MAE and log MAE as mandatory robustness co-criteria and the reversal rule capping status at `PARTIALLY_SUPPORTED`.

### 5.5 Null variables, frozen

To measure the false-positive rate, **200 decoy features** are added: each is a real probe whose donor order has been permuted, carrying the symbol of a gene drawn uniformly at random from the network, under a fixed seed. Retaining a decoy is an error.

**Corrected in v1.4.2, before any reported run, because the reasoning first written here was wrong.** The original text said NetPK's decoy rate must equal its rate on the corresponding real genes, since an outcome-blind ordering cannot tell a destroyed signal from an intact one. The first execution showed NetPK retaining **zero** decoys, and the explanation is not discrimination: a decoy carries the symbol of a *randomly drawn* network gene, which is almost always far from the seed and was therefore never a candidate for the top of the ordering. The decoy block measures susceptibility to a **spurious correlation**, which only an arm that reads the endpoint can have. It is retained for that purpose and the mechanistic arms' zero must not be reported as a merit.

**Destroyed-signal diagnostic, added in v1.4.2 — this is the test that does expose the weakness.** The target gene's own column has its donor order permuted. The gene keeps its place in the mechanistic ordering while carrying no information at all, so a mechanistic arm must keep choosing it and must lose its advantage, while an arm that reads the endpoint should choose something else. Reported per endpoint: the mechanistic arm's MAE at `k = 1` on the wrecked matrix, the univariate arm's MAE on the same matrix, and how many folds the univariate arm still picked the dead gene. This measures the price of being outcome-blind and is the direct link to H4.

**Corrected in v1.4.6, after running it: the diagnostic is interpretable only
on an annotation-qualified universe.** The permutation destroys the signal in
the one column that stands for the target gene. It cannot destroy it in a
*differently named* feature that measures the same transcript, and this array
contains such features. Measured on `CYP2D6_activity`:

| universe | the univariate arm's replacement | |ρ| | MAE on the wrecked matrix |
|---|---|---|---:|
| unrestricted | `CYP2D17` in **10 of 10 folds** — verdict `conflicting`, identity establishable from no public source, refused by every qualified pool | 0.45–0.51 | 96.86 |
| pool A | **nine different genes in ten folds**, all admissible | 0.27–0.32 | 108.88 |
| *(reference: the clinical baseline is 109.64, the blind mechanistic arm 110.15)* | | | |

`CYP2D17` is the name of a macaque CYP2D gene. Whatever its probe hybridises,
it tracks the same measurement as the target column closely enough to replace
it, so the unrestricted diagnostic let an arm "resist" a destroyed signal by
finding **the same molecule under another name**. Under the qualified pool that
feature is refused, the arm finds a different weak association in almost every
fold, and it performs at the clinical baseline — like the blind arm.

Three consequences, kept separate:

* **What survives unconditionally** is the *detection* claim: the univariate arm
  never selected the dead gene, in 100 folds across ten endpoints. An arm that
  reads the endpoint notices that a column has stopped carrying information.
* **What does not survive** is reading its lower error as robustness of the
  *method*. On a qualified universe the recovery disappears, so the earlier
  reading of this diagnostic as the mechanistic arm's main weakness was
  conditional on an inadmissible duplicate and must not be repeated.
* **No code is changed, and none should be.** The restriction repairs the
  diagnostic by itself, because the duplicate lives in the refused stratum: on a
  qualified universe every feature resolves to a gene and the single column
  standing for the target gene is the only copy of it. Permuting "every feature
  of the gene" would not have helped in the unrestricted universe, since the
  duplicate is not annotated as that gene — which is the annotation problem
  itself, not a parameter of the diagnostic.

**And a consequence for H2d.** Where no admissible feature carries the target's
signal, the destroyed-signal experiment cannot discriminate between arms: every
arm is at the baseline, so a repair criterion that requires an interval to
exclude zero fails for want of anything to find rather than for want of a method
that could find it. That reading is mandatory wherever the criterion is reported
unmet.

### 5.6 Analytical admissibility of the candidates, and the shared universe

Added in v1.4.5, after the primary H2a run and **before** any pool-restricted
run or any H2d run. It changes no threshold, margin, comparator or primary
contrast, and recomputes no primary number.

**The finding that forced it.** A gene symbol in this archive's annotation does
not identify the molecule measured. Checking the archive's own annotation
against the independent platform record for all 40,638 features, with
accessions resolved to genes through NCBI's public accession-to-gene map before
anything is compared, gives:

| verdict | features | what is established |
|---|---:|---|
| `concordant_human` | 22,279 | two independent annotations resolve to the same human gene |
| `single_source_human` | 4,018 | one source names a human gene; nothing agrees or disagrees with it |
| `symbol_names_another_gene` | 572 | the accessions resolve to one gene while the analysis's symbol is the approved symbol or a listed synonym of another |
| `platform_self_contradictory` | 610 | the platform's accession and its own gene assignment name different genes |
| `conflicting` | 1,473 | the sources name different genes; identity unestablished |
| `undetermined` | 10,679 | no source resolves the feature to a human gene |
| `technical_control` | 1,007 | the array's own control feature |

**Two of these refusals exist because the first graded version was still too
permissive**, and both were corrected before any pool-restricted run:

- a verdict must not attribute a measurement to the wrong gene. Where the
  analysis's symbol is the approved symbol or a listed synonym of a *different*
  human gene, the claims are mutually exclusive and the feature is refused.
  Where the symbol is simply unknown to the current gene reference — a clone or
  locus placeholder, or a symbol a spreadsheet turned into a date — it carries
  no competing claim, the accessions decide, and the measurement is kept with
  the gene's current symbol recorded beside it. The two groups are of very
  different size, 632 against 5,268, so refusing both would have discarded
  about five thousand real measurements for a naming-vintage reason;
- the platform record can contradict itself, its accession naming one gene and
  its gene assignment another. Internal coherence is checked **before** any
  agreement between sources is computed, because a record that contradicts
  itself cannot confirm the archive. 610 features, 607 of which the previous
  rule admitted.

**A declared limit that follows, and is not repaired here.** 4,929 admissible
features carry a name the gene no longer has. The analysis universe is keyed on
the array's own symbols, so a mechanistic ordering keyed on current gene symbols
cannot reach those features under those names. Re-keying the universe would make
roughly five thousand further measurements reachable by the mechanistic arms,
which is a methodological change rather than a correction; it is not done, and
the limit is recorded so it is not mistaken for an absence.

No verdict checks hybridisation specificity, because neither source distributes
the oligonucleotide sequence. "Identified" therefore means **identified by
annotation**, and the corollary is stated with equal force: **an undetermined
feature is not a wrong measurement.** It is a measurement whose biological
meaning the available annotation does not pin down, and a statistically
informative variable can remain biologically undetermined. Analytical quality
and predictive usefulness are reported separately and never merged.

**The pools, frozen.** Built at the **feature** level and only then grouped by
the symbol an analysis uses, so a gene with one admissible and one inadmissible
feature keeps the admissible one; dropping the gene whole would discard a real
measurement because a *different* feature of the same gene is badly annotated.

| pool | admits | symbols | features |
|---|---|---:|---:|
| **A** qualified human | concordant + single source | 21,634 | 26,066 |
| **A strict** | concordant only | 19,352 | 22,125 |
| **B** = A ∩ network | concordant + single source | 11,030 | 14,230 |

The admissible feature identifiers, the rule version and the SHA-256 of every
source are committed in `registry/QUALIFIED_PROBE_POOLS.json`, so each
individual decision is reconstructible rather than summarised.

**How a pool is used, frozen.**

- The pool applies to **every arm** of a comparison, never to one of them. A
  restriction applied to the mechanistic arm alone would change the comparison
  rather than qualify it.
- The feature standing for a gene is chosen per training fold (§5.4), so a
  selection recorded as a symbol is not traceable to a measurement. The feature
  actually used is recorded fold by fold and checked against the pool; a
  selected feature outside the pool stops the run instead of being reported.
- A technical control carries no transcript, so its survival into a pool is a
  construction error. Its absence is asserted by the runner, not assumed.
- **H2a and H2d are comparable only on one universe.** Any text comparing them
  must cite runs carrying the same pool file and pool name. An unrestricted
  result and a pool-restricted result are not comparable.
- The `undetermined` and `conflicting` strata are kept and reported, never
  deleted. Each analysis states which strata it admits.
- **The 200 decoy features are not pool members and never become one.** A decoy
  carries the symbol of a network gene the array does not measure, so it cannot
  intersect a pool; the identifier recorded for a decoy selection is a reserved
  marker rather than a feature. A decoy selection is a selection error and is
  counted separately from any annotation problem.

**Status of the pool-restricted runs.** The primary H2a run stands as executed
on the unrestricted universe, under protocol v1.4.3, and is not recomputed. The
pool-restricted runs are **prespecified sensitivity analyses**: pool A answers
whether the comparison holds among analytically admissible variables, pool B
answers it among variables both admissible and reachable by the mechanistic
ordering. Neither is the honest answer on its own, because each restriction
moves the comparison in both directions at once — it removes unannotated
features a statistical arm could exploit, and it removes candidates from every
arm.

**One endpoint loses its target.** `CYP1A2_activity` has **no admissible
feature**: all four of its features are `conflicting`. Its H2a result is
withdrawn, it is excluded from the pool-restricted primary tables, and the
activity measurement remains usable only as a **secondary** endpoint for
mechanisms other than its own transcript. The withdrawal holds under the
stricter and fairer rule, not only under the first one.

### Criteria

- Final MAE at equal panel dimension; selection frequency over repeated folds; Jaccard stability between selections; false-positive rate on null variables; causal-variable recovery in synthetic known-truth experiments; ability to reject the deliberately wrong panel; abstention rate when evidence is insufficient.
- **Methodological success:** NetPK is non-inferior within a prespecified **5 % relative MAE margin** against **every** statistical comparator separately, *and* shows better stability or better selection specificity. **Corrected in v1.4.3.** The earlier wording, "non-inferior to the best statistical comparator", was implemented by picking the statistical method with the lowest out-of-fold error and then contrasting against it. That uses the test folds twice and credits the winner of three draws, so the comparator looks better than any single method truly is. Each method is now contrasted separately and the margin must hold against all three, which is strictly harder than beating a test-chosen winner and needs no extra computation.
- **Superiority** requires a favourable paired difference whose uncertainty is compatible with a real improvement. If the interval is wider than the margin, the status is `INCONCLUSIVE`, not equivalence. A non-significant difference is never reported as equivalence.
- **Failure:** NetPK repeatedly retains incorrect pathways, or degrades beyond the margin.
- **Stability is not scored naively.** At a fixed `k`, NetPK's Jaccard across folds is **1.0 by construction**, because the ordering never reads the outcome. That is a structural property and not a merit, and reporting it as one would be misleading. What is compared is the data-driven arms' stability, and NetPK's stability when its `k` is chosen in the inner folds — where it can vary.
- **Status cap.** A `functional_in_vitro` result cannot promote the clinical H2 claim above `PARTIALLY_SUPPORTED`. It is reported as the **functional analogue** of H2. No public cohort currently permits the clinical-PK form of this experiment at a panel width where selection is meaningful.

**Mandatory honesty check, now enforced rather than asserted.** The selection must be *produced by the stated procedure*, not written into the code after inspecting results. Two tests enforce this: one parses the selector's source and fails if any gene symbol appears in it, so a retained set cannot be smuggled in as a literal; the other re-derives the decision from its written definition — the walk by direct linear solve, the null by the same stratified relocation, the cut by the same estimator — and asserts the module agrees. The former hard-coded retain/exclude step is removed.

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

| 1.3.0 | 2026-11-07 | **Recorded late: this row was missing while §1 already carried 1.3.0 text.** Jarmusch's omic layer is joined per donor through the GNPS processed output, so H1 has two public clinical-PK tests rather than one, and the second is already negative. No threshold, margin, comparator or primary contrast changed. |
| 1.4.0 | 2026-11-08 | **§5 rewritten: the mechanistic selection is now executable.** The hard-coded retain/exclude step is removed; the ranking is computed from a downloaded public network with seeds obtained automatically from the dataset's own column name, a public pathway annotation or a public chemical-gene resource, and two tests enforce that no gene is named in the selector and that the decision is reproducible from its written definition. Records the finding that a propagation score cannot be given a calibrated significance threshold, measured as one or two genes retained for all ten endpoints and the same number with the wrong enzyme; the panel size is therefore chosen in training folds and the claim of autonomous sizing is withdrawn in advance. Adds the decoy-variable specification, the leave-one-centre-out secondary design with its inferred-not-documented label, the eighth comparator (target transcript alone), and the statement that ARIANA cannot carry H2. Done before any enzyme-activity value was read. No threshold, margin, comparator or primary contrast changed. |

| 1.4.1 | 2026-11-08 | **§5.4 corrected before any run, for a reason independent of any endpoint.** Measured covariate completeness in the curated phenotype: weight 144 of 466 donors, body mass index 124 of 466, age complete, sex missing for one. The baseline is therefore age and sex only; keeping weight and body mass index would discard about seventy per cent of the cohort. Adds the frozen missing-value rule and the one-feature-per-gene universe, so that panel dimension means the same thing in a gene-ordering arm and a probe-ordering arm. No threshold, margin, comparator or primary contrast changed. |

| 1.4.2 | 2026-11-08 | **§5.5 corrected: the stated reason for the decoy block was wrong.** The protocol predicted that a mechanistic arm would retain decoys at the rate it retains real genes. The first execution retained zero, and the cause is not discrimination: a decoy carries the symbol of a randomly drawn network gene, almost always far from the seed and so never a candidate for the top of the ordering. The block in fact measures susceptibility to a spurious correlation, which only an arm reading the endpoint can have, and the mechanistic zero must not be reported as a merit. Adds the destroyed-signal diagnostic, which permutes the target gene's donor order so the gene keeps its rank while carrying no information: the test that does expose the cost of being outcome-blind. Corrected before any reported run. No threshold, margin, comparator or primary contrast changed. |

| 1.4.3 | 2026-11-08 | **Two implementation defects corrected after external review, before any H2a result was reported.** (1) The repeated cross-validation overwrote itself: out-of-fold errors were stored in one vector per arm, so with five folds by two repeats each donor was written twice and the reported losses came from the last repeat alone. Results are now stored per repeat and averaged per donor before any loss or bootstrap, which also keeps the donor as the independent unit. (2) The non-inferiority criterion named "the best statistical comparator", and the implementation chose it by out-of-fold error -- using the test folds twice and crediting the winner of three draws. The three statistical methods are now contrasted separately and the 5 % margin must hold against every one of them. The run in progress was stopped and discarded rather than reported. No threshold, margin, comparator or primary contrast was weakened. §5.3 also splits the second question in two at the review's request: H2a is relabelled as contextualised mechanistic selection rather than autonomous discovery, since the endpoint's label supplies the enzyme; the new H2b asks whether the chain can identify the mechanism from the drug instead, evaluated as mechanism recovery against an independently curated consensus so the ground truth does not come from the chain's own input; the integration question keeps its content and becomes H2c. |

| 1.4.4 | 2026-11-08 | **Adds H2d, hybrid selection, and resolves a second label collision.** Knowledge proposes a shortlist of three times the budget from the outcome-blind ordering; the training donors then keep, demote or abstain, by partial rank correlation against the clinical baseline calibrated on a permutation null. Prompted by the H2a destroyed-signal result on one endpoint, so the rule is frozen in its own dated document which records that provenance, caps the generating endpoint at `exploratory`, makes the other nine the confirmatory set, and declares in advance the failure mode of collapsing onto the univariate selector. The review that proposed it called it "H2c", a label already used for the integration question, so it is H2d. No H2a number is recomputed or reinterpreted, and no threshold, margin, comparator or primary contrast changed. |
| 1.4.5 | 2026-11-09 | **Adds §5.6: analytical admissibility of the candidate features, and the rule that H2a and H2d share one universe.** Prompted by two external reviews and by a result of the project's own: a gene symbol in this archive does not identify the molecule measured, so the candidate features are graded by how strongly two independent annotations agree on a human gene, with accessions resolved to genes before comparison, with a symbol that names a *different* human gene refused and an outdated name kept with the gene recorded beside it, and with the platform record checked for internal coherence before any agreement between sources is computed. Three pools are frozen, built feature by feature and only then grouped by symbol so a mixed gene keeps its admissible features, with the complete admissible lists, the refused strata and every source checksum committed. The pool applies to every arm or to none; the feature actually used is traced fold by fold and a selection outside the pool stops the run; the undetermined and conflicting strata are kept rather than deleted, and an undetermined feature is explicitly not a wrong measurement. `CYP1A2_activity` has no admissible feature under the stricter rule either, and its result stays withdrawn. The 200 decoys stay outside every pool. The limit that 4,929 admissible features carry a name the gene no longer has is declared rather than repaired, because re-keying the universe would be a methodological change. Specified after the primary H2a run and before any pool-restricted or H2d run; no primary number is recomputed and no threshold, margin, comparator or primary contrast changed. |

| 1.4.6 | 2026-11-09 | **§5.5 corrected after running the destroyed-signal diagnostic: it is interpretable only on an annotation-qualified universe.** The permutation destroys the target gene's column, not a differently named feature measuring the same transcript, and this array has one. On `CYP2D6_activity` the univariate arm replaced the destroyed column with the feature annotated `CYP2D17` in 10 of 10 folds at |ρ| 0.45–0.51 — a feature whose identity no public source establishes and which every qualified pool refuses — reaching 96.86 against a clinical baseline of 109.64. Restricted to pool A the same arm picks nine different genes in ten folds at |ρ| 0.27–0.32 and lands at 108.88, beside the blind mechanistic arm's 110.15. The detection claim stands: the univariate arm never selected the dead gene in 100 folds. The reading of its lower error as methodological robustness does not, and is withdrawn. No code changes: the qualified pool repairs the diagnostic because the duplicate lives in the refused stratum. Where no admissible feature carries the target's signal, the experiment cannot discriminate between arms, and a repair criterion reported unmet must be read as "nothing to find", not "no method could find it". No threshold, margin, comparator or primary contrast changed. |


### Validity of the 1.1.0 amendment

A change to admissibility is only legitimate if it cannot have been chosen for the result it produces. Three facts are recorded, and the first is enforced by a test:

1. **No HLC model had been fitted under this protocol when the ruling was made.** No HLC performance figure existed under these rules, so the ruling could not have been selected on one. `tests/test_public_data_evidence_registry.py::test_amendment_records_that_the_ruling_preceded_the_runs` fails if that claim is removed from the registry.
2. **No criterion moved.** Every threshold, margin, comparator and the primary contrast are byte-identical to 1.0.0. The amendment widens the input set; it does not lower a bar.
3. **No hypothesis changed class.** HLC's endpoint is still `functional_in_vitro`. It gains H3/H4/H5 and the functional analogue of H1; it still cannot satisfy H1, and ARIANA remains the only H1-eligible dataset. Becoming retrievable does not turn an in vitro measurement into clinical pharmacokinetics.

The direction of this amendment is worth stating plainly: it **enlarges** the evidence base in which the hypotheses can fail. HLC is the cohort where the existing finding is *unfavourable* to broad panels — a targeted CYP panel beats the broad multi-CYP/POR/transporter panel in 9 of 10 endpoints — so admitting it adds opportunities for refutation, not for confirmation.
