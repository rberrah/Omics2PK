# Phase 1 erratum — four corrections

**Date:** 2026-11-07
**Protocol:** v1.1.0 → v1.2.0
**Registry:** v1.1.0 → v1.2.0
**Origin:** external review of the Phase 1 freeze, relayed by the repository owner.

No threshold, margin, comparator or primary contrast changed. All four corrections are to **provenance and to the status of evidence**. Each was verified against a primary source before being accepted; one of the four reviewer claims was also extended, because verification turned up a limitation the review had not identified.

---

## E1 — A second Müller study was missed, and a negative conclusion built on that miss is withdrawn

**What the Phase 1 audit said.** No public cohort can test within-subject temporal stability of the solanidine ratios. Sampling-time effect, food effect and intra-individual stability were recorded as `ABSTAIN`.

**What was wrong.** The audit checked Müller 2024, the phenotyping cocktail study (doi `10.1002/cpt.3106`, DRKS00028922), established that it releases nothing, and generalised from that single study to the whole question. It never searched for a *different* study by the same group.

**Verified.** A separate publication exists and releases its individual data openly:

> Müller JP, Friske D, Hasoumi M, Miller B, Sarömba J, Tupiec J, Koczera P, Laurentius T, Clavel T, Stingl JC, Just KS. *Solanidine metabolic ratios for CYP2D6 phenotyping after a standardized potato meal.* Front Pharmacol. 2026 Oct 2;17:1886427. doi `10.3389/fphar.2026.1886427`

Confirmed independently: the GitHub API returns the public repository `JulianPeterMueller/Solanidine-standardized-potato-meal` with two files; Crossref returns the title, journal, volume 17, article 1886427 and publication date 2026-10-02 with 11 authors. The workbook `Solanidine standardized potato meal.xlsx` (20,965 bytes, SHA-256 `68d53071…1950`) was downloaded and all four sheets parsed.

**What it contains**, verified cell by cell rather than from the abstract:

| Sheet | Content |
|---|---|
| Participant data | 20 participants, all CYP2D6 genotyped, with predicted phenotype, activity score, sex, age, BMI, eGFR, and a 3-day potato-fasting flag |
| Plasma concentration data | 3 analytes × 20 participants × 9 timepoints (0, 0.5, 1, 1.5, 2, 3, 4, 5, 29 h); **180/180 cells present for each analyte, no missing values** |
| Glycoalkaloid concentrations | The meal itself assayed — α-chaconine, α-solanine and solanidine across 4 aliquots, with mean and SD, so the **precursor dose is quantified** |
| Microbial metabolism | Two stool samples incubated with α-solanine or α-chaconine under active versus inactivated microbial community, in triplicate |

**Consequence.** The `ABSTAIN` on three H5 items is lifted: sampling-time effect, food effect (15 of 20 fasted on potato products, 5 not), and intra-individual temporal stability. Registered as `muller_potato_meal_20`, status `QUALIFIED`, objectives `["H5"]`.

**Four constraints, all recorded before any run.**

1. **It replicates nothing about pharmacokinetics.** No drug is administered. `endpoint_class` is `biomarker_only` and `no_drug_administered` is `true`. It tests whether the marker behaves reproducibly, never whether it predicts exposure. H1 remains unreplicable on public data.

2. **It is not an independent laboratory — a limitation the review did not raise.** Crossref shows **Sarömba, Stingl and Just author both this study and ARIANA**. The participants are different people (20 healthy volunteers aged 22–27 versus 47 cardiology patients), so donors are not double-counted, but shared method bias would produce agreement exactly as a true effect would. Any concordance between the two cohorts is reported with this stated. The review described this dataset as permitting "une véritable validation biologique indépendante"; that is too strong, and the registry records it as *a different population measured by the same group*.

3. **Limit-of-quantification imputation is already in the source, and is flagged only by bold font.** The workbook states its own rule: below-LOQ values were imputed as half the LOQ for solanidine (0.01 ng/ml) and 4-OH-solanidine (0.02 ng/ml), or the lowest measured value for SSDA. Extracted with openpyxl font inspection:

   | Analyte | Imputed cells | Share |
   |---|---|---|
   | solanidine (ng/ml) | 27 / 180 | 15.0 % |
   | 4-OH-solanidine (ng/ml) | 15 / 180 | 8.3 % |
   | SSDA (counts) | 9 / 180 | 5.0 % |

   **Participant 1 — the cohort's only poor metaboliser (`*4/*4x2`, activity score 0) — is below the limit of quantification at all 9 timepoints for both metabolites.** The subject the marker most needs to identify has a wholly imputed ratio. Any discriminative performance at the poor-metaboliser end rests on imputed floor values. Because imputed cells carry no data flag, a value-only pipeline would treat those floors as measurements; `scripts/qualify_muller2026_potato_meal_source.py` reads the font, and `tests/test_muller2026_potato_meal_qualification.py::test_plain_value_read_would_miss_the_imputation` pins the reason.

4. **Activity scores are not comparable with ARIANA's without alignment.** This cohort scores `*41` as **0.25**, the post-2021 CPIC convention. The Rüdesheim 2020 mapping used by the mechanistic module is built on the grid (0, 0.5, 1.25, 1.5, 2, 3) and predates that change. The mapping is not applied across the two cohorts until the convention is explicitly aligned, and that alignment is reported as a decision rather than performed silently.

**One mechanistic observation worth keeping.** The stool incubations show solanidine is formed only with an active microbial community. Solanidine supply is therefore gut-microbiota dependent, which cuts both ways and both directions are tested: it is a reason a *measured* ratio can outperform a genotype estimate, and a reason the ratio's own reproducibility may be limited.

---

## E2 — The ARIANA uniqueness claim is withdrawn

**What the protocol said.** "Exactly one public dataset carries a clinical PK endpoint with matched omics *and* pharmacogenetics — ARIANA n=47."

**What was wrong.** That phrasing implies Jarmusch has no pharmacogenetics. It does. The repository's own prior audit records Europe PMC supplement `CTS-13-972-s001.docx` as a 15-row × 11-column individual table covering **14 participants** with CYP2C19, CYP3A4, CYP3A5 and CYP2D6 genotypes and phenotypes, alongside midazolam CL/F (13 person-paired). PK and PGx are matched per individual there.

**Superseded in part by E5 below.** This section originally said the omic join was unresolved. The authors' own input `GNPS_quant.csv` is indeed absent from their repository, but the crosswalk *was* reconstructed from the GNPS processed output and the H1 comparison has already been run. See **E5**. The corrected uniqueness claim below stands unchanged, because it rests on biomarkers being *named and quantified*, which Jarmusch's ions are not:

> ARIANA n=47 is the only public cohort in which clinical exposure, pharmacogenetics and **named, quantified** biomarkers are matched at the individual level today.

The registry summary continues to report `h1_eligible_datasets == ["ariana_metoprolol_47"]`, though for the reason given in E5 rather than the one originally stated here.

---

## E5 — The Jarmusch omic join is not pending: it was done, and H1 failed there

**Added 2026-11-07, registry v1.3.0.** This corrects E2 above, and corrects the external review on two numbers.

**What is actually in hand.** The authors' input `GNPS_quant.csv` is absent from their repository, but the GNPS *processed output* `qiime2_output/qiime2_table.qza` is held locally — an HDF5 BIOM with **6,349 unique MS1 features across 1,198 analytical sample IDs**. Filtering `MS_SampleType=sample`, `blood_plasma`, day 1, time 0 min yields **13 distinct biological subjects** with predose plasma MS1, each matched to day-1 targeted midazolam oral CL/F, CYP3A4 and CYP3A5 calls, sex and measured weight. **12 of the 13** also have day-9 CL/F. Analytical replicates were verified consistent on patient/day/time before collapsing to one observation per person-day, and are never counted as independent people.

**Two corrections to the review.** It described the matched set as *10 volunteers* and the results as *"certaines améliorations numériques apparaissent, mais leurs intervalles d'incertitude recouvrent zéro"*. Against the repository's own record, both are wrong:

| | Day 1 (n=13) | Day 9 (n=12) |
|---|--:|--:|
| sex + weight | 6.875329 | 3.984289 |
| + CYP3A4/CYP3A5 | 7.343413 | 4.601193 |
| + PGx + predose MS1 PCA2 | 9.888019 | 4.748841 |
| **omic ΔMAE vs PGx** | **+2.544607 (+34.65 %, worse)** | **+0.147648 (+3.21 %, worse)** |
| descriptive 95 % interval | **[+0.4421, +5.9959] — excludes zero** | [−0.4273, +0.6899] — covers zero |
| individual errors improved | 1 of 13 | 4 of 12 |

There were **no improvements**. The omic arm was worse on both days, and on day 1 the interval excludes zero **on the harmful side**. Worse still, 25 person-to-metabolome shuffles gave a median augmented MAE of 8.0924 and **all 25 beat the true-pairing model** on day 1 — random pairings outperformed the real one. MAE values are in the source's numeric scale, as the physical CL/F units are not securely verified, so they are comparable only within a row.

**Why the negative result is conservative rather than fragile.** The GNPS/MZmine feature alignment was computed across the complete assay set **before** the leave-one-out split, so the feature universe is transductively informed by the held-out samples. That deviation can only *favour* the omic arm. The omic arm lost anyway, which makes the finding stronger, not weaker.

**Why it still cannot qualify for H1.** The frozen protocol requires preprocessing inside training folds only. A feature universe built once across all samples cannot serve as a *qualifying* H1 test, however informative it is. The entry therefore gains `H1` among its supported objectives but stays `ACQUIRED`, which keeps it out of `h1_eligible_datasets`.

**What this changes for the programme.** H1 now has **two** public clinical-PK tests rather than one, and the second has already failed:

- **ARIANA** n=47 — named, quantified biomarkers; exploratory re-analysis; not yet run under the frozen protocol.
- **Jarmusch** n=13/12 — unidentified ions, 2 global PCs, tiny n; **already negative**, with a documented preprocessing deviation that worked in the omic arm's favour.

A plausible reading, to be tested rather than assumed: 6,349 features projected onto 2 global principal components in 12–13 subjects is a weak, non-specific covariate, so the failure may be a property of *unselected* high-dimensional omics rather than of omics as such. That is precisely the H2/H4 question, and this cohort is the designated hard case for it. Four probes measured in the **same** 14 donors remain **one** biological cohort, and day 1 versus day 9 are paired occasions, not independent subjects.

---

## E3 — ARIANA conclusions are capped at exploratory

**Accepted in full.** This project has analysed ARIANA repeatedly under earlier internal versions, with model choices made while results were visible. Freezing the protocol stops *further* result-guided decisions; it cannot convert a re-analysis into an independent validation.

Recorded as `prior_analysis_exposure` with `max_claimable_status: exploratory`. No ARIANA result may be reported as confirmatory however favourable it is, or however strictly the frozen rules were followed. Note the distinction from the registry's `confirmatory_admissible` field, which describes whether a *source* may enter a confirmatory analysis, not whether a given *analysis* is confirmatory.

---

## E4 — The NLME limitation is stated precisely

**Was:** "An NLME likelihood is not identifiable from a single observation per subject."

**Now:** "One observation per subject does not make a nonlinear mixed-effects model impossible, but it leaves between-subject variability and residual error not separately identifiable, and it cannot support estimation of individual PK parameters from ARIANA alone. No NLME variability claim is made from this cohort."

The operative conclusion is unchanged — MAP-BE stays excluded by construction, and Gate D3 still cannot be marked compliant from ARIANA — but the reason is now the correct one.

---

## Not a correction: HLC

The review's arbitration on HLC — a free account with no individual authorization is admissible, a committee or author request is not — matches the ruling already made in v1.1.0, independently and on the same grounds. v1.1.0 went further than the review by typing the gate from the Synapse API rather than from documentation: the requirement on `syn4499` is a `SelfSignAccessRequirement` (id 9604913, `isCertifiedUserRequired=false`, `isValidatedProfileRequired=false`), and the ACL grants READ+DOWNLOAD to `AUTHENTICATED_USERS` and READ only to `PUBLIC`. The review's suggestion to "vérifier les conditions de chaque fichier" is satisfied: the requirement's subject is the container `syn4499`, so it governs all four files, which an unauthenticated query against the child entity `syn4629` confirms.

The owner then re-downloaded the whole folder through that gate as a single 136,114,988-byte bundle. All four inner archive SHA-256 values recompute identical to those already recorded, and the bundle's `manifest.csv` supplies Synapse's own MD5 values and metadata (`curatedPhenotype` numSamples 836, `curatedExpressionLiver` numSamples 465, Agilent Custom 44K, tissueType Liver). The documented access route is therefore reproducible, and the 2026-10-18 acquisition is independently confirmed.

---

## Validity of this erratum

The H5 items that `muller_potato_meal_20` makes testable were frozen in **protocol v1.0.0, before this dataset was known to exist** — verifiable from the v1.0.0 commit. No criterion was written or altered with this data in view, which is the strongest form the prespecification argument can take: the analysis plan predates the discovery of the data it will run on.

`tests/test_public_data_evidence_registry.py::test_amendment_records_that_the_ruling_preceded_the_runs` enforces the equivalent claim for the v1.1.0 access ruling.

## Reporting consequence

Of the three findings the Phase 1 report presented as settled, one is now wrong, one is narrowed, and one is strengthened:

| Phase 1 statement | Status after this erratum |
|---|---|
| The phenotype-recovery result is the publication's own finding, reproduced exactly; a positive control, not evidence for H1 | **Unchanged.** |
| No public cohort can test intra-individual temporal stability of the solanidine ratios | **Wrong.** Müller 2026 can, for the marker — never for PK. |
| ARIANA is the only public cohort with PK, omics and PGx | **Narrowed** to named, quantified biomarkers; Jarmusch has PK and PGx matched. |
| H1 cannot be replicated on public data | **Unchanged**, and reinforced: the newly found cohort administers no drug. |
