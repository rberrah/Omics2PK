# Original RACE randomized L-carnitine study: real baseline multi-platform metabolomics → 90-day PD outcome

**Original citation:** Puskarich M.A. et al., *Clinical and Translational Science* (2021), DOI [10.1111/cts.13088](https://doi.org/10.1111/cts.13088), PMID 34216108, [author public source repository](https://github.com/UMichNMR-Metabolomics/CaPS_TimeZero), official study [Metabolomics Workbench ST001319](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST001319&StudyType=NMR&ResultType=1). Original source commit SHA **`2dcacb5d82fe07cf625ae0e4edab34e51546f3e6`** (2023 author repository).

**True original-data execution (green):** [GitHub Actions #39165799873](https://github.com/rberrah/Omics2PK/actions/runs/39165799873). This is a **new genuine original human PD study**, not generated/virtual data or borrowed statistics.

## Exact original source: three public patient-level files, matched by studyid ONLY

- Author `RACE_Demographics.csv`: n**247** human participants in analysis release; `studyid`, actual randomized arm and dose (0 placebo, 6/12/18 L-carnitine dose groups), age, sex, baseline total SOFA, survival days.
- Original `raw_nmr.csv`: n**228** individually keyed original baseline NMR assay records, 27 actual pretreatment metabolite concentrations.
- Original `raw_ac.csv`: n**236** individually keyed original **T0 PRE-treatment** LC-MS acylcarnitine assay records, 24 metabolite concentrations. All original timepoint labels in the source are T0; no post-treatment biological signals included.
- Exact original study-ID joined patients with original age/sex/SOFA and 90d outcomes on **225 different biological patients**; 119 died by day 90 according to original `survivaldays <= 90` published parent study definition, 68 placebo, 157 active-dose arm. Only this exact 225-person intersection is evaluated in **ALL models**. Do not compare source NMR n228 versus LCMS n236 separate sampling frames.

Independent recorded original files' SHA256 checksums (only metadata, **no actual people/measurements stored or printed**):

```text
RACE_Demographics.csv  eea11d13fe7431df52cc53b94cbc94a5002cf13de24092de0d83e9843f43186d
raw_nmr.csv          dbbd1f5ae37677b94a452e1399e91fdec720b65fa3341800d6290708f9236a94
raw_ac.csv           6b4990f71a4ab0affbcb11ee3f28b285029ae95d60e6b624abcb94e9ec9e8e06
```

The author repository's software LICENSE is **AGPL-3.0**. Our code is independent Python and original patient rows are fetched as a pinned original-source stream in CI, not copied to Omics2PK. Honor the primary data and software source terms separately.

## Original source-fitted models

Five original PERSON-level heldout StratifiedKFold splits, seeded `20261018`, stratified by observed 90d outcome and treatment/placebo. Fixed logistic regression C=1; `age, sex, original T0 SOFA, randomized dose, active-treatment indicator` in clinical baseline. `SimpleImputer(median) → log1p → StandardScaler → PCA(n=3)` within **training donors only** for each of the two authentic metabolomic platforms. No feature selection or outcome-informed PCA on heldout subject rows. All models use **the same 225 persons and the same exact outer folds**.

- `Clinical`: baseline clinical and assigned treatment;
- `Clinical+NMR`: + three unsupervised *baseline* PCs from the original 27 pretreatment metabolites;
- `Clinical+LCMS`: + three baseline PCs from original 24 acylcarnitines;
- `Clinical+NMR+LCMS`: both separate 3-PC assay representations;
- `Clinical+NMR+LCMS×Treatment`: exploratory treatment×6-PC interactions. No individualized causal effects inferred from this nonlinear/prognostic model.

**Both NMR and LCMS assays are metabolomics**. They are two measurement platforms within one molecular modality, **not transcriptome + proteome**, so call this "multi-view metabolomics" not true multi-omics.

## Actual measured 90-day PD outcomes, same 225 original human trial participants

Smaller heldout Brier/logloss is better; larger AUROC better.

| Prespecified source-fitted model | AUROC | Brier | log loss |
|---|---:|---:|---:|
| Clinical (age, sex, T0 SOFA, randomized original dose/treatment) | 0.724037 | 0.215019 | 0.624823 |
| Clinical + LCMS 24 acylcarnitines → 3 train-only PCs | 0.742746 | 0.206803 | 0.603147 |
| Clinical + NMR 27 metabolites → 3 train-only PCs | **0.763755** | **0.199114** | 0.591150 |
| Clinical + both NMR and LCMS platforms | 0.759553 | 0.199573 | **0.588787** |
| Clinical + both + treatment×PC interaction | 0.744094 | 0.207318 | 0.610574 |

**Original donor-heldout PAIRWISE Brier difference vs clinical, with exploratory 1,000 donor-row paired 95% percentile bootstrap of already heldout errors:**

| Comparison | ΔBrier vs same-person clinical | Descriptive paired bootstrap 95% |
|---|---:|---|
| + LCMS | −0.008216 | [−0.021677, +0.006011] |
| + NMR | **−0.015905** | **[−0.032478, +0.002000]** |
| + NMR + LCMS | −0.015446 | [−0.033879, +0.003183] |
| + both + treatment interactions | −0.007701 | [−0.027171, +0.011019] |

**Interpretation:** original baseline NMR metabolite information is numerically associated with better out-of-fold 90-day mortality risk discrimination, but the interval contains zero and is not proof of statistical/external incremental predictive value. **Adding the LCMS second metabolomics view to NMR alone does NOT reduce heldout Brier nor increase AUROC**, though combined model's log loss is slightly lower. More complexity is not necessarily better. The treatment interaction model does not establish a real patient-specific mortality *treatment benefit*: randomized treatment effect heterogeneity requires a prospectively justified causal design and external verification, not just a risk classifier with interactions.

## Critical distinctions and limitations

- This is **true clinical PD / observed event prediction**, NOT L-carnitine in-vivo plasma PK, no individual drug C(t), AUC, CL/F or mechanistic concentration effect.
- Mortality is an indirect disease/treatment endpoint strongly affected by baseline disease severity; clinical+NMR improvement is mostly **prognostic** unless true treatment-effect contrast proven. Source original `survivaldays` is used as an observable 90-day binary outcome; future survival analyses should explicitly include followup completeness, censoring, deaths as events and proportional/competing hazards assumptions as appropriate.
- Cohort is from one trial and source center(s); random splits within one cohort are **internal**, not external independent trial validation. Outcome-dependent imputation/time normalization risks controlled within our analysis, but source metabolomic quality preprocessing may already have been performed before release.
- Clinical reference does not include every comorbidity/vasopressor dose/albumin/biochemical prognostic factor. Any claim NMR independently adds to a **fully specified external prognostic gold standard** requires those data and matched n.
- Research scripts use source-only subject IDs as join keys internally then emit **aggregate metrics only**. No original source rows or identifiers committed or printed.
- The original author's paper and source method remain authoritative for interpretation of all transcript/omics data, metabolite annotations, statistical endpoints and human subjects.

**Science status:** `ACTUAL_225_PATIENT_REAL_PD_PREDOSE_BENCHMARK_REPRODUCED; TWO_METABOLOMICS_PLATFORMS_NOT_TRUE_MOLECULAR_MULTIOMICS; NMR_ADDED_PROGNOSTIC_PERFORMANCE_EXPLORATORY; TREATMENT_SPECIFIC_BENEFIT_NOT_PROVEN; CLINICAL_EXPOSURE_PK_NOT_MEASURED`.

## Reproduce responsibly

```bash
python -m pip install -e .
python -m pip install 'numpy==2.2.6' 'scikit-learn==1.7.2'
python -m scripts.validate_original_public_caps_pd
```

`src/omics2pk/caps_original_pd.py`, `tests/test_caps_original_pd.py` (fabricated source rows only), `.github/workflows/caps-original-metabolomics-pd.yml`. The original public author CSVs are downloaded into RAM ephemerally and never saved, uploaded or distributed as public outputs.
