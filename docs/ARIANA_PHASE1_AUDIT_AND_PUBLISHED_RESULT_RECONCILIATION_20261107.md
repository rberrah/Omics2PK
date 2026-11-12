# Phase 1 audit — ARIANA v5 reconciled against the source publication

**Date:** 2026-11-07
**Scope:** audit and freeze step of the five-objective validation protocol. No new model was run; this document reconciles existing claims with the published source and records what must be withdrawn or relabelled.

## 1. The ARIANA data are public, and are the data already in use

The source publication's data availability statement names a public repository:

> "Original data were uploaded to a repository accessible under the following link: https://github.com/katjajust/Solandine-Metoprolol-Correlations"

| Item | Value |
|---|---|
| File | `Solanidine_metoprolol_ratios_PGx.xlsx` |
| Size | 13,997 bytes |
| SHA-256 | `7397b13e43cd546b93f83a20b2dc84b6d49e97633173d56711c39973beedfe37` |
| Rows × columns | 47 × 16 |

The checksum is **identical** to the file already used by the local ARIANA analysis. Two consequences:

1. ARIANA is `open_public` and therefore admissible for confirmatory analysis under the public-only rule. It is the **only** public dataset carrying a clinical PK endpoint with matched omics *and* pharmacogenetics.
2. No patient-level file needs to enter this repository. The dataset is reproducibly obtainable from a URL plus a checksum, and `.gitignore` already excludes `*.xlsx`.

Columns available: `ID`, `Allel 1`, `Allel 2`, `AS`, `CNV`, `c(Metoprolol)[ng/ml]`, `Daily Dose Metoprolol`, `c(metoprolol)weighted [ng/ml]`, `lnMR OH-metoprolol/metoprolol`, `lnMR OH-solanidine/solanidine`, `lnMR SSDA/solanidine`, `Age (years)`, `male sex`, `BMI`, `weight`, `CYP2D6 inhibitor intake`.

## 2. Correction: the "biological validation" result is the publication's own headline

A previous internal report (v5) presented the recovery of the reference CYP2D6 phenotype by solanidine-derived activity as **"the strongest result of the project"**. That framing was wrong and is withdrawn.

The source publication's title is *"Very strong correlations of metoprolol- and solanidine-derived metabolic ratios in a real-world cohort with a stable metoprolol drug regimen"*, and this correlation **is its primary published finding**. Recomputing it from the public file reproduces the published values exactly:

| Correlation with `lnMR OH-metoprolol/metoprolol` | Published | Recomputed here | Match |
|---|---|---|---|
| `lnMR OH-solanidine/solanidine` | r = 0.907, R² = 0.823 | r = +0.907 [0.838, 0.947], R² = 0.823 | exact |
| `lnMR SSDA/solanidine` | r = 0.946, R² = 0.896 | r = +0.946 [0.905, 0.970], R² = 0.896 | exact |
| CYP2D6 activity score | r = 0.500, R² = 0.250 | r = +0.500 [0.249, 0.688], R² = 0.250 | exact |

The publication also already reports the univariate association with the endpoint of interest, which likewise reproduces exactly:

| Correlation with the metoprolol trough | Published | Recomputed here |
|---|---|---|
| OH-solanidine ratio | r = −0.513, R² = 0.263 | r = −0.513, R² = 0.263 |
| SSDA ratio | r = −0.533, R² = 0.284 | r = −0.533, R² = 0.284 |
| CYP2D6 activity score | r = −0.288, R² = 0.083 | r = −0.288, R² = 0.083 |

**What this changes.** The *direction* of the project's central ARIANA finding — a solanidine-derived marker tracks CYP2D6 activity far better than the genotype activity score, and correlates with the metoprolol trough about three times more strongly in R² terms — was **already established in the source publication**. It is not a discovery of this project.

**What remains this project's own contribution**, and must be stated in these terms only:

- an **out-of-fold predictive** evaluation (the publication reports in-sample correlations, not held-out prediction error);
- a **multivariable** evaluation against clinical and pharmacogenetic baselines;
- comparison against **named competing estimators** on identical folds (Ridge, XGBoost, PopPK covariate arms, reduced PBPK);
- the **mechanistic parameterisation** of the marker as enzyme activity rather than as a regression covariate.

The reproduction of the published correlation is therefore a **positive control on the implementation** — valuable, because it proves the module encodes the published biology correctly — and is relabelled as such throughout. It is not evidence for H1.

### 2.1 Why one reported number differed from the publication, verified

The v5 report gave the SSDA arm as r = 0.937 (R² = 0.878) where the publication gives 0.946 (0.896). The cause is a modelling choice, not a data difference, and it is now verified numerically:

The activity transform is `activity = exp(gamma * (lnMR - lnMR_ref))`, capped at 2.13, the top of the Rüdesheim relative-activity scale.

| Analyte | Patients at the 2.13 cap | Uncapped | Capped |
|---|---:|---|---|
| OH-solanidine | 2 / 47 | r = 0.907, R² = 0.823 | r = 0.907, R² = 0.822 |
| SSDA | **15 / 47** | r = 0.946, R² = 0.896 | r = 0.937, R² = 0.878 |

Uncapped, both reproduce the publication exactly. The cap binds for 2 patients on OH-solanidine and for **15 of 47** on SSDA, which accounts for the entire discrepancy.

**This is a finding in its own right, and it is unfavourable to the current parameterisation.** A cap that binds for a third of the cohort means the reference-centred exponential overshoots the physiological activity range for SSDA. The reference level is the training-fold mean of genotype-normal metabolisers without an inhibitor (n = 16 here). Either the exponent `gamma = 1` is too steep for SSDA, or the semi-quantitative SSDA scale (no authentic standard) is not on the same multiplicative footing as a true activity ratio. This is logged as a prespecified H5 item, not silently retuned.

### 2.2 The composite panel is not better than the single published ratio

| Activity estimate | r | R² |
|---|---:|---:|
| OH + SSDA, geometric mean, capped | 0.947 [0.906, 0.970] | 0.896 |
| SSDA alone, uncapped (= published) | 0.946 [0.905, 0.970] | 0.896 |

The two-analyte panel reported in v5 as the best estimate is **numerically indistinguishable from the single published SSDA ratio**. The apparent benefit of combining two analytes was an artefact of comparing a capped single analyte against an uncapped-equivalent combination. No multi-omic advantage is demonstrated at the level of phenotype recovery.

### 2.3 A genuine, modest new observation: the genotype arm was improved, not the omic arm

The v5 comparison (omic R² 0.896 versus genotype R² 0.470) did **not** compare like with like against the publication. The publication correlates the **raw linear activity score**. The mechanistic model correlates a **nonlinear, phenoconversion-corrected relative activity**. Decomposed:

| Genotype-based activity estimate | r | R² |
|---|---:|---:|
| Raw activity score (what the publication reports) | 0.500 | 0.250 |
| Rüdesheim nonlinear k_cat mapping of the activity score | 0.646 | 0.417 |
| Rüdesheim mapping + ΔAS = −0.54 for inhibitor intake (9/47 patients) | 0.685 | 0.470 |

So the **nonlinear activity mapping nearly doubles the explained variance of a genotype-only activity estimate** (0.250 → 0.417), and the phenoconversion correction adds a further increment (→ 0.470), on the 9 of 47 patients taking a CYP2D6 substrate or inhibitor.

This supports two of the transferred constants independently of any predictive result, and it is the one genuinely new biological observation in this reconciliation. It also means the omic-versus-genotype gap reported in v5 is **narrower than the published gap**, because the genotype baseline was strengthened. The v5 direction was conservative; only the novelty attribution was wrong.

### 2.4 Dose normalisation materially strengthens the marker-endpoint association

Not reported in the publication's main correlations, and recorded here because it affects endpoint definition:

| Correlation with the trough | Unweighted | Dose-weighted |
|---|---|---|
| OH-solanidine ratio | r = −0.513, R² = 0.263 | r = −0.625, R² = 0.391 |
| SSDA ratio | r = −0.533, R² = 0.284 | **r = −0.703, R² = 0.494** |
| CYP2D6 activity score | r = −0.288, R² = 0.083 | r = −0.375, R² = 0.141 |

Dose normalisation raises the SSDA association from R² 0.284 to 0.494. Since the cohort receives two distinct daily doses (95 and 190 mg/day), part of the unweighted residual variance is simply dose. The frozen protocol keeps the **observed trough** as the primary endpoint, because that is the quantity a clinician reads; the dose-weighted trough is retained as a prespecified secondary endpoint. Choosing between them after seeing model results is prohibited.

**Minor source note.** The publication reports p = 0.499 for the activity-score/trough correlation alongside a confidence interval of [−0.531, −0.001] that excludes zero. The recomputed value is p = 0.0499. The published figure appears to be a decimal typographical error; the correct reading is borderline significance, not absence of association.

## 3. Audit of the v5 ARIANA model inventory

Free parameters fitted on the 47 patients, and the provenance of every fixed constant:

| Arm | Free parameters fitted on ARIANA | Fixed constants and their source |
|---|---:|---|
| Ridge (any feature set) | 1 penalty (inner folds) + coefficients | none |
| XGBoost (any feature set) | tuned hyperparameters + all tree structure | none |
| M4/M5 mechanistic | **1** (baseline level, then locked) | f_m,CYP2D6 = 0.80; relative activity by activity score (Rüdesheim 2020 Table 3); ΔAS = −0.54 (external n = 88); allometric exponents 0.75/1; sex effect exp(0.572) (Eugene 2016); γ = 1; blend weight 0.5 |
| P1 PopPK published | **1** (base level) | sex 0.572; CYP2D6 effect (1−f_m)+f_m·activity with f_m = 0.80 |
| P2 PopPK estimated | 3 | allometry 0.75 |
| P3 PopPK estimated + omic | 4 | allometry 0.75 |
| Reduced PBPK | **1** (multiplicative scale factor, closed form, then locked) | full OSP drug parameter set; non-CYP2D6 intrinsic clearance solved on the published PM/NM exposure ratio of 5 |

Two anchoring disclosures stand and are not weakened:

- The mechanistic model's `f_m = 0.80` is a literature constant, not a fitted value, but it is also not validated for ARIANA.
- The reduced PBPK's non-CYP2D6 intrinsic clearance is **solved against a published exposure ratio**, because the OSP "specific clearance" parameter is not transposable to a reduced structure. The reduced PBPK is therefore **not free of empirical anchoring**. Its hepatic bioavailability remains free and emerges at 0.35 against 0.40 published, which is the module's genuine internal consistency check.

### 3.1 Items the audit confirms must not be claimed

- **No MAP-Bayesian estimation on ARIANA.** One concentration per patient makes it circular; the measured optimism factor is 3.6. Excluded by construction in the frozen protocol.
- **No NLME fit on ARIANA**, hence **Gate D3 cannot be marked compliant** from this cohort regardless of model performance.
- **The OSP observed profiles are not an independent external validation.** They were digitised by the authors of the model whose parameters are transferred, and they are group-mean profiles. The protocol labels this a parameter-transfer check.
- **The STRING/Random-Walk-with-Restart step is prioritisation, not autonomous selection.** The ranking is computed, but the retain/exclude decisions are hard-coded. H2 cannot be `SUPPORTED` on an autonomy claim until an executable rule reproduces the retained set mechanically.
- **The previous selection-dimension boundary is superseded.** It rested on too few noise realisations; the protocol now requires at least 100 repetitions before any boundary statement, and forbids asserting a universal p/n ≈ 2 threshold.

## 4. Acquisition audit: the public-only rule retires the existing queue

Every item in `docs/EXTERNAL_DATA_ACQUISITION_MATRIX_20261024.md` depends on an author request or a governance decision:

| Former priority | Access | New status |
|---|---|---|
| Zhu 2022 tacrolimus repeated TDM | author request | `ABSTAIN` |
| Phapale n = 29 tacrolimus | author request | `ABSTAIN` |
| BIOMARGIN ↔ ABIS/ISBA linkage | governed | `ABSTAIN` |
| Zhang n = 100 / Ohtsuki n = 17 hepatic | author request | `ABSTAIN` |
| Rosuvastatin n = 40 processed CSVs | author request | `ABSTAIN` |

**The project's four previously open questions cannot be answered by acquisition in this phase.** They must be answered with data already in hand, or abstained from. This is a consequence of the constraint, not a failure of the search.

### 4.1 The CYP2D6/solanidine axis cannot be replicated on public data

Checked and verified on 2026-11-07:

| Candidate | n | Access | Status | Note |
|---|---:|---|---|---|
| Müller phenotyping cocktail | 20 | no release (Europe PMC `hasData=N`, `hasSuppl=N`, no PMCID, publisher 403) | `ABSTAIN` | The only source with 15 timepoints per subject |
| Wollmann TDM cohort | 839 | author request | `ABSTAIN` | No drug concentration measured at all |
| Kiiski Finnish volunteers | 314 | author request + GDPR/consent limits | `ABSTAIN` | No metoprolol PK |
| Estonian Biobank PK recall | 114 | governed application | `ABSTAIN` | Excluded by the public-only rule |

**Consequences, recorded before any run:**

1. The H1 **replication** criterion is unattainable for this axis in this phase. An ARIANA result can be an exploratory success and nothing more. It must not be presented as confirmatory.
2. **Intra-individual temporal stability and the sampling-time effect are `ABSTAIN`.** ARIANA documents neither absolute sampling time nor repeated sampling; the one cohort with 15 timepoints releases no individual data. The H5 items that depend on them are withdrawn rather than approximated.
3. Replication must be sought as the **functional analogue** of H1 on HLC, EPA and Wegler — a different biological axis and a non-PK endpoint, which is explicitly weaker evidence and is reported in a separate column.

### 4.2 Verified local holdings

| Source | Artefact | SHA-256 | Verified |
|---|---|---|---|
| ARIANA | `Solanidine_metoprolol_ratios_PGx.xlsx` | `7397b13e…fe37` | downloaded and hashed 2026-11-07 |
| HLC | `curatedExpressionLiver.zip` | `35a80d53…e376` | recomputed, matches repository record |
| HLC | `curatedGenotype.zip` | `9add3449…1690` | recomputed, matches repository record |
| HLC | `curatedPhenotypesynapse.zip` | `141187cd…15a7` | recomputed, matches repository record |
| HLC | `human_liver_cohort.zip` | `498dd06d…0c37` | recomputed, matches repository record |
| Jarmusch | GNPS feature-based molecular networking archive | `2047ea20…c026` | hashed 2026-11-07, member list inspected |

Checksums are abbreviated to first and last four hex characters; the full 64-character values are in `registry/PUBLIC_DATA_EVIDENCE_REGISTRY.json`.

All four HLC checksums match the values recorded in `docs/HLC_ORIGINAL_400_LIVER_CYP_RNA_GENOTYPE_RESULTS_20261018.md`, so that document's provenance claims verify.

**HLC access tier requires a user ruling.** Synapse is open to anyone who registers an account and accepts data-use terms; no committee reviews the request. That is materially different from the Estonian Biobank, but it is not "retrievable without an account" either. The registry classifies it `public_registered` and **excludes it from confirmatory-admissible tiers pending an explicit ruling**. Since HLC carries most of the existing functional evidence, this ruling determines how much of the programme can proceed.

## 5. What Phase 1 changes about the claim boundary

**Withdrawn or relabelled:**

- the phenotype-recovery result, from "strongest result of the project" to "positive control reproducing the source publication's primary finding";
- the two-analyte omic panel, from "best activity estimate" to "indistinguishable from the single published ratio";
- the selection-dimension boundary, pending ≥ 100 repetitions;
- the Random-Walk-with-Restart step, from selection mechanism to prioritisation aid.

**Confirmed and strengthened:**

- the data in use are the published public data, checksum-identical;
- the published correlations reproduce exactly, so the implementation is sound;
- the nonlinear activity mapping and the phenoconversion correction each improve a genotype-only activity estimate, independently of any omic;
- the 2.13 activity cap binds for 15 of 47 patients on SSDA, which is a real calibration problem now logged rather than hidden;
- all inherited HLC provenance claims verify by recomputed checksum.

**Newly recorded as impossible in this phase:** replication of the CYP2D6/solanidine axis on public data, intra-individual temporal stability, the sampling-time effect, and any acquisition-based answer to the four previously open questions.
