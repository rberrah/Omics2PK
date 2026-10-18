# First original human clinical E3 benchmark: metformin predose GC-TOF → AUC0–24 / CLr (2026-10-18)

**Primary original scientific source:** Rotroff DM et al., *Frontiers in Pharmacology* 2016 DOI https://doi.org/10.3389/fphar.2016.00135, PMC4906013. The official Europe PMC public original paper supplementary archive is https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4906013/supplementaryFiles .

**Actual original-data run with numerical results:** https://github.com/rberrah/Omics2PK/actions/runs/39145479441 . Source-identity/crosswalk run: https://github.com/rberrah/Omics2PK/actions/runs/39144998287 . Original source schema and synthetic privacy tests run successfully. **This is a real original-source quantitative human in-vivo PK evaluation, not a simulated fixture.** All patient-level clinical and omics observations are downloaded ephemerally, with only aggregated evaluation and source digests printed.

## Source and exact clinical specimen join

The study dosed **33 humans** with metformin and acquired original untargeted **GC-TOF molecular serum/plasma samples at A (before metformin dosing), B and C (after dosing)**. The original supplementary contains 5 Excel workbooks and 16 worksheets. The complete 118-row untargeted metabolomic source worksheet includes **99 authentic labelled individual samples** of form FOUR-LETTER-SOURCE-PREFIX + FOUR-DIGIT-SUBJECT + -A/-B/-C (33 complete A/B/C triplets), plus 18 source technical/sample-QC rows which are NEVER considered biological people. Each original sample A contains **384 numeric metabolomic feature measurements**.

Exact original clinical source uses 33 unique 4-digit person identifiers. The four-digit component of each original molecular sample ID matches **all 33 of the 33 original patient identifiers**, uniquely, without sorting subjects or guessing row order. Two original primary xlsx binaries are individually SHA256 pinned:
- **Clinical PK source:** 76ef6d62e1809f37b7fcf0d9a07e38eb55a67331b23bd058890911547fa9681a
- **Predose metabolomic source:** 1debb35eac7c22383031d4aaeaa47c0362d57609956f16fa4420626efd4a8fe6

The original metformin-PK subject-level worksheet reports:
- **AUClast (0-24)**, 32 valid independently measured subject-level records (actual in-vivo drug exposure AUC, distinct from the separately reported glucose-response AUC).
- **CLr_renal clearance**, 32 valid independently measured subject-level records.
- Age, Gender, Weight and other PK explanatory variables in the same patient-numbered table.

The original clinical header omits the **physical AUClast measurement unit**. This report keeps errors in **numeric original source table units**, and DOES NOT arbitrarily label AUC MAE ng·h/mL or renal CLr MAE mL/min until article methodology/codebook confirms physical units. One original clinical subject lacks both endpoints; it is never imputed.

## Fixed benchmark (exploratory E3)

- **Independence**: 32 original unique humans, one *A PRE-DOSE* molecular specimen per subject, two measured clinical drug outcomes; 32 leave-one-person-out outer tests for each outcome.
- **Reference**: age + sex + weight, fixed Ridge(alpha=1), log1p observed drug endpoint, original-units backtransformation.
- **Omics-augmented**: reference + 3 *unsupervised principal components* computed from **384 predose features**, with sample missing-value imputation, z-score scaling and PCA fitted **strictly on the outer training donors**. No supervised feature selection outside folds, no use of B or C samples. Model and PCA component count fixed without selecting on these original outcome results.
- **Secondary**: 3 pre-dose omics PCs only, fixed Ridge(alpha=1), same outer folds.
- **Uncertainty**: descriptive 1,999 paired bootstrap re-samples of the already heldout per-person absolute-error differences. This is NOT a nested/refit bootstrap, nor independent external validation.
- **Environment**: Python 3.11, NumPy 2.2.6, scikit-learn 1.7.2 pinned in the E3 original-source CI. Synthetic-only tests separately verify output confidentiality, 32-person requirement and full fold-local PCA/normalization.

### Actual source-heldout measured clinical PK results

All errors below are in **original numerical source table units** (physical units still require verification).

| Original observed clinical endpoint | n | Model | MAE | RMSE | Bias | ±20% predictions |
|---|---:|---|---:|---:|---:|---:|
| Metformin AUClast 0–24 | 32 | Clinical age/sex/weight | **2826.0472** | 3264.2271 | −285.7454 | **56.25%** |
| Metformin AUClast 0–24 | 32 | Clinical + 3 predose molecular PCs | 2825.0240 | 3424.2538 | −48.6526 | 50.00% |
| Metformin AUClast 0–24 | 32 | Predose molecular PCs only | 2673.1230 | 3111.4288 | −162.4417 | 53.125% |
| Metformin CLr | 32 | Clinical age/sex/weight | **115.1113** | **131.2252** | −17.7205 | 37.50% |
| Metformin CLr | 32 | Clinical + 3 predose molecular PCs | 118.1449 | 137.7472 | −15.0257 | 37.50% |
| Metformin CLr | 32 | Predose molecular PCs only | 121.3379 | 159.0275 | −23.4347 | 46.875% |

**Prespecified primary incremental test (clinical + omics vs clinical only):**
- Metformin AUClast 0–24: ΔMAE = **−1.023245** numeric source units (≈ **0.04%** MAE reduction, negligible), exploratory paired 95% bootstrap interval **[−592.810886 ; +681.066112]**; n18 donor predictions with improved absolute error vs 14 worse, and a worse aggregate RMSE and P20.
- Metformin renal CLr: ΔMAE = **+3.033606** (≈ **2.64% worse** vs clinical), exploratory interval **[−11.968748 ; +21.604118]**; 20 donor improvements vs 12 worse, but larger losses in some donors caused the worse mean MAE.

Do not claim "omics improves in-vivo clinical PK" from these primary comparisons. That would be unsupported. The favorable omics-only AUC MAE is an exploratory model comparison involving a different predictor set, **not** the prespecified incremental improvement over clinical baseline, and is not a reason to re-select the primary endpoint/model.

## Scientific interpretation and restrictions

**Achievement:** the project has a genuine, reproducible, completely public original 33-person clinical OMICS↔patient PK dose/exposure source, exact original human identifiers and **32 original measured drug AUC / renal clearance events**, and successfully executes train-only donor-heldout E3 models. No author-controlled closed Google Drive is required for this n33 cohort.

**Limitations:**
- Only **32** complete individuals, 384 features and three PCA components; these high-dimensional small-n generalization estimates are highly uncertain.
- This is **internal** leave-one-human-out, not independent biological external validation, no leave-centre-out replication, no bootstrapped refit.
- Genotypes, transporters, formulation/comorbidity and kidney function baseline comparators are not yet fully explored. AUC / CLr source physical units need direct assay methods verification before communicating absolute magnitude.
- The shared study's predose endogenous metabolites are valid prospective candidate predictors; causality cannot be inferred and the subsequent B/C profiles are excluded from prediction.
- No SWEET, XGBoost, M2 combined omics, M3 mechanistic PopPK or M4 network has been run in this benchmark. This is a **controlled preliminary E3 reference** from which those predeclared methods may be developed WITHOUT tuning conclusions using the same heldout benchmark as a final external test.
- The 2016 original cohort differs from the 2022 metformin n17 and from rosuvastatin n40 and Wegler n54. Neither metformin cohort may be silently pooled, and original drug AUC cannot be substituted with post-metformin glucose AUC.

## Reproduce from the immutable original data

- Package: src/omics2pk/rotroff33_clinical_e3.py
- Real original data runner: scripts/validate_rotroff33_original_metformin_clinical_e3.py
- Donor-key contract: scripts/qualify_rotroff33_original_metformin_patient_keys.py
- Synthetic tests: tests/test_rotroff33_original_patient_keys.py and tests/test_rotroff33_clinical_e3.py
- CI: .github/workflows/rotroff33-e3-original-pk-benchmark.yml

Execute from checked-out Omics2PK (requires internet access to Europe PMC):

```bash
python -m pip install -e .
python -m pip install 'numpy==2.2.6' 'scikit-learn==1.7.2'
python -m scripts.validate_rotroff33_original_metformin_clinical_e3
```

The command emits only original source digests and **aggregate** metrics. Individual original person IDs, metformin AUC/CLr and metabolomic observations are not serialized in repository source, GitHub Actions logs or artifacts.

**Result evidence status: `AUTHENTIC_MEASURED_CLINICAL_E3_INTERNAL_N32_CONFIRMED; MOLECULAR_ADDED_VALUE_NOT_CONFIRMED; EXTERNAL_VALIDATION_PENDING`.**
