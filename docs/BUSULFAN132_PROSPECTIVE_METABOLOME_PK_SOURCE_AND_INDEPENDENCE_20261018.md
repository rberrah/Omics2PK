# Busulfan 132 genuine clinical predose metabolomics → measured PK clearance: independent-person qualification

**Scientific priority date: 18 October 2026**  
Original [McCune et al., *Clinical Pharmacology & Therapeutics* 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC9888309/), DOI [10.1002/cpt.2794](https://doi.org/10.1002/cpt.2794), PMID 36369996, PMCID PMC9888309. Relevant comparator studies: [Lin et al. 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5562150/), DOI 10.1007/s11306-016-1106-6, and [McCune et al. 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8214873/), DOI 10.1021/acs.jproteome.0c00599. **This is genuinely measured target-drug clearance**, not in-vitro CYP turnover, an artificial AUC, a contemporaneous drug-intensity source peak, or just a PD outcome.

## Source biological truth vs eligibility of individual original rows

| Source paper | True original biological study and paired outcome | Original numeric patient-level metadata/data actually acquired in Omics2PK |
|---|---|---|
| **2023, n132** | 132 **independent PEOPLE** with immediately predose plasma; **96 of these same people** also sampled `2-week-predose`. Total **228 original metabolomic specimens**, **841** measured endogenous plasma metabolites. Same human observed IV busulfan concentration–time observations fed Bayesian NONMEM model → individual measured plasma-anchored CL normalized to normal fat mass, original `mL/min/kg NFM` | **NOT obtained.** Article data, figure correlations, pathway/FDR tables, and source supplemental article PDFs are NOT 132×841 matrix and 132 matched PK target rows |
| **2016, n108** | 108 original humans, **147 predose plasma specimens**, 1885 untargeted ions. Split 59-discovery and 88-validation **SAMPLES**; independent PEOPLE between arms must be independently checked, not inferred from 59+88 | **NOT obtained**: original LCMS ion × original matched donor Bayesian CL table not verified in official author supplement |
| **2021, n84** | 84 HCT recipients, **230 longitudinal plasma metabolite samples**, 763 endogenous features, measured GVHD/relapse and targeted PK-guided IV busulfan | **NOT obtained**: no 84-person pre-dose multi-omic×individual busulfan CL + PD event rights-cleared quantitative table |
| **Jarmusch 2020, n14** | True original CYP probe drug plasma PK + CYP genotypes + untargeted metabolome / microbiome across 14 controlled volunteers | **NOT obtained**: original 14-person PK targeted C(t)/AUC Table S2 and molecular sample-ID crosswalk not qualified (separate previous MassIVE/PMC audit) |

**Critical study design finding:** In McCune 2023, published 13-metabolite model was fit on 132 immediately-pre-busulfan samples and *transferred* to 96 specimens collected two weeks earlier **in the same 96 original humans already represented in the 132-sample training dataset**. Their reported `R²≈0.40` in these 96 earlier samples demonstrates transfer across **sample timepoints**, NOT independent-person holdout or independently replicated external clinical predictive performance. This is a source-confirmed methodological distinction, not evidence that the published association is invalid. No authentic GSTA1/GSTM1 genotype data from these same 132 people were verified in the article’s public source tables; their discussion of PGx cannot substitute a matched genotype-versus-metabolomics comparator.

The 2023 clearance was calculated by Bayesian estimation using actually measured plasma busulfan observations and a previously developed PopPK prior of **1,610 patients / 12,380 observations**. Those 1,610 people are PRIOR population model data, **not an added 1,610 original omics sample cohort in this manuscript**. Do not count 228 plasma samples as 228 independent donors, and don't fit a patient-level molecular comparison on repeated profiles with a sample-wise split.

## Official source availability audit actually executed

Script `scripts/audit_original_busulfan_and_jarmusch_PK_releases.py`: retrieves official Europe PMC fulltext JATS and checks DOI/PMCID and supplementary **file metadata only**, never patient records. CI <https://github.com/rberrah/Omics2PK/actions/>. A DOI that validates the actual article does not establish the original 841×228 molecular source files; a publisher supplementary PDF with 13-metabolite coefficients is not the source original measured intensities × person PK.

Original 2023 article's PMC Associated Data section lists **Appendix S1, S2, S3 (PDF), Figures S1/S2 (PDF), Table S1 (PDF)**. We found no advertised original supplementary CSV/XLSX 228×841 matrix in this listing. The source may exist with the investigators and/or clinical study repository; **the absence of an attached CSV on this route is not proof it does not exist**. The other studies' source constraints were checked separately in [source catalogue](DEEP_PUBLIC_HUMAN_OMICS_PK_PD_COHORT_CATALOGUE_20261018.md).

**Original data still needed** to run this comparison:
- Per-human irreversible research pseudonym and explicit source-document crosswalk;
- authentic original continuous `busulfan_CL_per_NFM`, source `mL/min/kg_NFM`, plasma PK model method, date and dose/covariate context;
- 228 original rows of 841 pre-dose metabolite concentrations/standardized relative abundance **with feature labels, original missing-value code, batch, collection timestamps, relative dose windows and QA flags**;
- original age, sex, normal fat mass/organ function, co-prescribed phenytoin vs levetiracetam, conditioning regimen;
- **separately, true validated genotype** including GSTA1 functional haplotypes and GSTM1 copy number/deletion from **these same people**, or other mechanism-relevant genetics. If PGx was never assayed, the genotype arm must report **NOT TESTABLE**, rather than default imputing "reference" in every patient.

## Implemented original-only supervised experiment — NOT run on authentic McCune rows yet

`src/omics2pk/busulfan132_predose_pk.py` now defines explicit original-source fields, validates n132 unique donors and n96 repeat predose pairs, validates original URI source DOI, plasma PK clearance physical units, 841 measured features with missing/censor checks, timezone-qualified all assayed-before-first-dose sampling and same-person clearance/covariate concordance. It rejects duplicate donor/sample IDs and source labels, adverse donor repeat mismatch, mixed unit/period or invented baseline genotype.

The original-source-only research function `benchmark_authorized_original_busulfan` will, after legitimate authorized input:

1. Group the **96** two-week samples as test targets in 5 held-PERSON folds. In each fold, exclude **every immediate-predose profile of each test person from training**; use other 132-donor immediate-predose training persons only. Keep 36 unmatched-only immediate-predose subjects as eligible training samples. All 96 original two-week samples receive a prediction from a model never exposed to ANY of that person's molecular or PK outcome.
2. On each TRAIN dataset, `SimpleImputer → StandardScaler → univariate SelectKBest(k=13) against TRAIN PK target`. Predictive feature selection is **within folds**, unlike reusing the source authors' original 13 selected on the complete n132 (which includes test humans). No target-derived feature selection in test.
3. Fit Ridge `α=1` on log(original measured CL); backtransform to clearance in original `mL/min/kg_NFM`. Compare **clinical** `age, sex, normal fat mass` vs **clinical+13 train-only metabolites** on identical held-human samples.
4. If and ONLY if original authenticated GSTA1/GSTM1 call codes and provenance are present in all source individuals, compare `clinical + verified selected GST PGx` vs `clinical+verified selected GST PGx+13 metabolites`. A selected GST functional call/copy number is **not comprehensive PGx**. All source-verified training/test batch/ethnicity/etc. candidates should be predeclared before materialized original data.
5. Output **aggregate-only MAE and MAPE**, person counts and leakage flags. No patient ID, data row, omic concentrations, model predictions, genotype calls or source file is exported or sent to GitHub.

The source auditor cannot itself prove that a caller-labeled `source_doi` row came from the real investigators; a human/source metadata review and checksum is separately necessary. Thus no public original 132-donor numerical result is claimed.

**Tests use entirely fictional 132×841 metabolite measurements with 96 repeated FAKE people** in `tests/test_busulfan132_predose_pk.py` and verify that the fixed-time two-week heldout individuals are excluded from the donor group's training source. A synthetic test PASS is **not** biological PK efficacy. GitHub CI uses no restricted source account and never attempts to bypass data-sharing conditions.

## Actual plan for genotype-vs-phenotype

Unlike a CYP3A enzyme, busulfan elimination chiefly proceeds through **glutathione conjugation via GSTA1**, with additional GST pathways. Prioritizing this non-CYP hepatic substrate avoids extrapolating liver CYP transcriptomic n400 results to busulfan clearance. **Mechanism-matched baseline is genuine GSTA1/GSTM1 genotyping, not CYP3A5 genotype**. Do not claim that donor microsomal CYP activity is the relevant phenotype for busulfan.

If the original cohort has no GST genotype, the direct **genotype-only vs measured predose metabolome** test cannot be performed in this population: a *separate* real clinical PK cohort with genetic calls is needed. The original 132 cohort may still validate the prospective incremental clinical + predose metabolome hypothesis.

Final step would be external patient center and true PK C(t) validation, with dose/event-based NLME covariates in Omics2PK, followed by HCT GVHD/relapse/toxicity in a **separately qualified** Omics2PD outcome stage. The same-source 2021 n84 and 2023 n132 cohorts may overlap in enrollment/calendar years: verify original human ID exclusivity before calling them external to one another.

**Current state:** `REAL_PATIENT_CL_STUDY_IDENTIFIED; N132_PAIRED_SAMPLE_DESIGN_CORRECTLY_INTERPRETED; TRAIN_PERSON_EXCLUSION_ENGINE_AND_SYNTHETIC_TESTS_IMPLEMENTED; MATCHED_ORIGINAL_METABOLOME_AND_PK_NUMERIC_ROWS_NOT_ACCESSIBLE; TRUE_SAME_PATIENT_GENOTYPE_NOT_DEMONSTRATED`.

### Additional verified Jarmusch supplementary XLSX access attempt

Unlike simply relying on original article abstracts, the original PMC [Jarmusch `Table S2` link](https://pmc.ncbi.nlm.nih.gov/articles/instance/7485946/bin/CTS-13-972-s002.xlsx) was followed to its exact primary NIH/PMC attachment route and a real fresh independent GitHub Actions binary HTTP request was executed. **Observed response: HTTP 200 but `Content-Type: text/html`, only 1,817 bytes**, SHA256 `416b8c58321f2c21781ab6045a3fb7b5c1277c7844b128e3269e401f5b591d9f`, **NOT an XLSX ZIP** (`PK\x03\x04` magic absent). This is a **download interstitial / access limitation** and not proof source Table S2 lacks individual quantitative PK data. The audit deliberately **abstains**, without parsing this webpage as 14 volunteer PK values. `scripts/audit_jarmusch14_official_original_table_s2_xlsx.py` and [original source CI #39167732406](https://github.com/rberrah/Omics2PK/actions/runs/39167732406) reproduce the result.

The original author's file can potentially be downloaded via institutional Wiley access, NIH PMC browser, or obtained under appropriate authorization. Only after the legitimate original XLSX is separately provided or retrievable and exact volunteer ID/drug PK/covariate columns validated can a real genotype+omics→PK benchmark for Jarmusch be executed. Neither guessing 14 AUCs nor scraping the source figure is acceptable.

