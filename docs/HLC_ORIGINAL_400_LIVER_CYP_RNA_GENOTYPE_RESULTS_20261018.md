# Real original matched human liver cohort: CYP RNA activity, genome SNP and whole-microarray experiments (18 October 2026)

**Evidence:** user-provided genuine original Sage Synapse HLC source files, analyzed **privately, without putting original ZIPs, participant IDs, genotypes, phenotype observations or expression matrices in GitHub**. Reproducible *code only* is included in this PR; all output below is aggregate. Citation to primary original HLC data generators: Schadt et al. 2008 *PLoS Biology* (PMID 18462017); Yang et al. 2010 *Genome Research* DOI [10.1101/gr.103341.109](https://doi.org/10.1101/gr.103341.109), and original Sage Bionetworks Human Liver Cohort. Data-sharing rights and the source terms explicitly prohibit attempts to identify donors.

## Source-proven original biological crosswalk (QUALIFIED)

Four supplied authorized original archives, never committed or redistributed:

| Original source role | Structure actually inspected | SHA-256 of original archive |
|---|---|---|
| curatedExpressionLiver | **40,638 liver Agilent microarray probes × 466 individual original donor sample ID columns**; probe annotations include CYP1A2, CYP2A6, CYP2B6, CYP2C8, CYP2C9, CYP2C19, CYP2D6, CYP2E1, CYP3A4/5, POR and other transporters | `35a80d532b967fdfbc73eaff571ddade17b90b30121ebe48465870245fe5e376` |
| curatedPhenotype | **18 original phenotype/activity fields × 836 unique liver individual IDs**; CYP activities and age, sex, anthropometrics, originally measured in HLM | `141187cd0e226ba43389bd10ab6973661ad935fa8bad51f9984f1833a2cf15a7` |
| curatedGenotype | original array **195 genotyped people**, genotype table 449,699 SNP rows (500,568 SNP feature annotations); only two *preselected* pharmacogenetic rsIDs verified available: rs776746, rs1799853 | `9add3449638a64b5db7ddb01ca85d8ae03d441679447cc133b9233dd91981690` |
| Historical human_liver_cohort | 836-row original phenotype table, 466-expression sample original 40k array, source README and Sage sharing terms, precomputed liver networks | `498dd06d106779adb4cc833a51a49be8f4e00cad1d07791cabc06d6c1b090c37` |

**Exact joins proved numerically:** all **466** expression sample IDs occur uniquely in the 836-donor phenotype data; both original phenotype source formats have **836/836 identical person IDs** and identical measured numeric CYP readouts on original common quantified samples (CYP1A2 and CYP2C9 original source 488 assays each, original source CYP3A4M/T n487, all matched exactly in source copy). The old historical and curated expression sources share the same **466** subject columns and **40,638 source probes**; 3,189 independent spot-checked numeric source sample/probe cells matched exactly. The genotyped source has **195** IDs in its original table, but only **178** overlap the original liver expression donors; do not call it 195 paired people.

**Original phenotype endpoint interpretation:** measured donor CYP probe substrate turnover assays from human liver microsomes, **not actual clinical PK plasma concentrations or human AUC**. The separate original README describes original activity assay units as nmol/min/mg protein; source raw numeric assay units, outliers and normalizations should be verified against Yang 2010 assay-specific supplement before interpreting absolute MAE magnitudes physically. Report primarily within-endpoint relative MAE.

## Scientific protocol: original n387–397 actual metabolic assay outcomes

Original donor ID split: fivefold shuffled KFold(random_state=20261018), **unique original biological donor ID within each group**. Fixed Ridge α=1 on log1p(CYP source activity), predictions backtransformed to original source scale; source-independent standardization and median imputation trained **ONLY on training donors** within each fold. Prespecified model arms for every endpoint:

1. age+sex baseline (bodyweight missing in most HLC humans, so not inserted or imputed as though observed);
2. the gene-specific **target CYP transcript microarray probes alone** (2–4 original annotated array probes);
3. age+sex+gene-specific CYP probes (**primary incremental contrast**);
4. age+sex+all measured CYP1A2/2A6/2B6/2C8/2C9/2C19/2D6/2E1/3A4/3A5, POR/ABCB1/SLC22A2 family panel of source-annotated probes (**secondary broader molecular hypothesis**).

No search for winning CYP panel size or result-guided gene combinations. Original platform probe-gene annotations are supplied by source; a microarray gene symbol does not establish orthogonal qPCR/protein quantities or whole-transcriptome measurement accuracy.

### Original actual donor-heldout results, measured CYP activity

MAE numerical original phenotype source units, smaller better.

| Observed probe CYP activity | n independent measured livers | Clinical age+sex MAE | Clinical + matching CYP transcript MAE | Relative gain | Clinical + broad CYP/POR/transport RNA MAE |
|---|---:|---:|---:|---:|---:|
| CYP1A2 | 397 | 477.107 | **399.395** | **16.29%** | 554.614 |
| CYP2A6 | 396 | 627.885 | **557.330** | 11.24% | 584.425 |
| CYP2B6 | 397 | 234.057 | **206.682** | 11.70% | 231.987 |
| CYP2C8 | 387 | 58.848 | **49.797** | 15.38% | 51.936 |
| CYP2C9 | 397 | 272.456 | **244.353** | 10.32% | 250.047 |
| CYP2C19 | 397 | 46.342 | 43.838 | 5.40% | **41.937** |
| CYP2D6 | 396 | 115.655 | **98.380** | 14.94% | 103.049 |
| CYP2E1 | 397 | 1281.911 | **1234.069** | 3.73% | 2899.692 |
| CYP3A4 midazolam activity | 396 | 851.109 | **719.247** | 15.49% | 733.487 |
| CYP3A4 testosterone activity | 396 | 3571.417 | **3144.195** | 11.96% | 3206.573 |

**Finding:** Target CYP RNA improves MAE over the same-cohort age/sex baseline for all 10 original measured enzyme activity assays. **In 9/10 endpoints the broad RNA multi-CYP/POR/transporter panel DOES NOT beat the smaller same-gene targeted RNA panel**. Only CYP2C19 benefits numerically from the broad panel in this fixed Ridge configuration; this is an exploratory observation, not proof other genes are biologically unimportant. The CYP2E1 panel is especially worse, consistent with overfitting/multicollinearity or outlier sensitivity, not proof a particular CYP pathway is irrelevant.

### Independent incremental RNA negative-control inference (10 original outcomes)

The original measured clinical age/sex and actual CYP activity are kept fixed. For each CYP we do **199 within-age-quartile × sex** random donor mismatches of the CYP RNA probes as a joint vector, refit on the same outer splits, and compute the empirical one-sided null p `(1+ # null ΔMAE <= observed ΔMAE)/(199+1)`. This approximates conditional exchangeability, is not a causal or external-validation test. Empirical resolution is 0.005 with 199 nulls, and Bonferroni nominal α=0.05/10=**0.005**. **All ten** return p=0.005, at the minimal attainable p; this result is not a guarantee of independent generalization.

Paired original heldout absolute-error bootstrap 1999 draws, **exploratory 99.5%** (Bonferroni-width) percentile loss intervals, without refitting every bootstrap model:

| Actual measured CYP endpoint | Original target RNA minus baseline ΔMAE | Descriptive 99.5% paired loss bootstrap interval |
|---|---:|---|
| CYP1A2 | −77.712 | [−113.245, −44.562] |
| CYP2A6 | −70.555 | [−110.415, −28.594] |
| CYP2B6 | −27.375 | [−42.034, −12.299] |
| CYP2C8 | −9.052 | [−12.301, −5.417] |
| CYP2C9 | −28.104 | [−44.931, −12.455] |
| CYP2C19 | −2.505 | [−4.462, −0.342] |
| CYP2D6 | −17.275 | [−28.123, −6.933] |
| CYP2E1 | −47.842 | [−99.409, +4.543] |
| CYP3A4 midazolam | −131.862 | [−188.052, −80.654] |
| CYP3A4 testosterone | −427.222 | [−625.198, −241.424] |

**Interpretation:** 9/10 descriptive 99.5% intervals exclude zero; CYP2E1 does **not**, despite its donor-RNA-shuffle empirical p reaching 0.005. They are distinct null targets and inferential assumptions. No external donor centre, randomized prospective trial, mechanistic causality or actual full omics RNA→protein test was done.

## Entire source liver transcriptome (40,638 probes) versus the matching CYP alone

On five **preselected** CYP endpoints we performed train-only `SimpleImputer(median) → StandardScaler → PCA(randomized, 8 components, fixed seed)` over the **complete original 40,638 array probes** inside each training fold (never fit on heldout individuals). Same original original fivefold splits and Ridge α1 fixed. Tested age+sex + CYP-target probes, age+sex + genome-wide 8 PCs, and age+sex + both.

| Actual source activity | n | Clinical + matching CYP RNA | Clinical + eight genome-wide PCs | Clinical + matching CYP RNA + eight PCs |
|---|---:|---:|---:|---:|
| CYP1A2 | 397 | 399.395 | 430.359 | 399.753 |
| CYP2C9 | 397 | 244.353 | 245.010 | **242.664** |
| CYP2D6 | 396 | 98.380 | 113.015 | **96.098** |
| CYP3A4 midazolam | 396 | 719.247 | 814.962 | **705.895** |
| CYP3A4 testosterone | 396 | 3144.195 | 3466.780 | **3057.984** |

**Finding:** global unsupervised transcriptome PCs **alone** do not beat CYP-target probes on these five primary phenotype endpoints. Adding PCs alongside the targeted CYP signal yields small improvements in 4/5, which must be tested for stability on genuinely new original human donors and with alternative predeclared dimension/hyperparameters; this one-run result is not external validation or proof of a particular expression network.

## Genuine two-layer molecular experiment: liver RNA + original germline SNPs

Original genotyping source: 195 people total, **178 share an exact original liver expression ID**, with **171–172 also having the CYP activity and age/sex baseline**. Original SNP features use numeric dbSNP identifiers (not the string `rs`). The original available variant labels were confirmed: `rs776746` for CYP3A5 and `rs1799853` for CYP2C9. Other commonly cited CYP variants were not represented by those specific rsIDs, so no gene-based PGx proxy alleles were imputed or fabricated. Source calls were CC/CT/TT (with some NA), encoded **dosage of T per observed array**, with no unverified claim of reference allele orientation or CYP3A5 expressor phenotype.

Prespecified identical-subset comparison: age+sex vs age+sex + CYP3A4 and CYP3A5 RNA for CYP3A4 phenotypes; CYP2C9 RNA for CYP2C9 outcome; age+sex + SNP; age+sex + RNA + corresponding SNP. Three fixed fivefold shuffled partitions (seeds 20261018/19/20) and mean heldout MAE across runs, with TRAIN-ONLY missing SNP imputation.

| Enzyme activity (original source units) | Original common n | Age+sex MAE | + RNA MAE | + SNP MAE | + RNA+SNP MAE |
|---|---:|---:|---:|---:|---:|
| CYP3A4 midazolam | 171 | 706.771 | 641.842 | 697.803 | **641.836** |
| CYP3A4 testosterone | 171 | 3452.323 | 3188.553 | 3453.227 | **3176.100** |
| CYP2C9 | 172 | 182.367 | **172.270** | 183.260 | 173.455 |

**Finding:** these particular genotype calls show negligible/inconsistent incremental predictive gain after target RNA. This **does not imply** that population pharmacogenetics or comprehensive germline genomic data have no value in human CYP clearance. Variant set, orientation, HLC ancestry composition, allele frequency, substrate, sample size, high-dimensional signal and study context all constrain the claim. Source does not contain matched CYP quantitative protein levels, and the germline SNP dosage is genuine genomics but NOT proteomics.

## Nonlinear validation: same original HLC liver people, fixed XGBoost M1-type baseline

An additional **actual original source** study test ran one prespecified nonlinearity comparison using XGBoost 3.1.3, exactly the original person IDs and 5 identical donor-heldout folds from the Ridge experiment. No original rows were transferred to GitHub. Hyperparameters were **fixed** before execution (150 estimators, max_depth=2, learning_rate=.035, min_child_weight=8, subsample=.85, colsample_bytree=.9, L2 reg_lambda=10, histogram trees; log1p endpoint). We did not tune models on outside CV results; this is an exploratory **fixed XGB** test, not the full nested/hyperparameter-optimized M1 XGBoost or official SWEET.

| True observed CYP substrate activity | Original n | XGB age/sex | XGB age/sex + same CYP RNA | XGB age/sex + broad CYP/POR/transport RNA |
|---|---:|---:|---:|---:|
| CYP1A2 | 397 | 477.582 | **395.501** | 402.832 |
| CYP2A6 | 396 | 634.484 | **554.625** | 564.629 |
| CYP2B6 | 397 | 233.002 | 213.797 | **213.218** |
| CYP2C8 | 387 | 59.112 | 50.323 | **50.112** |
| CYP2C9 | 397 | 276.338 | **247.648** | 248.800 |
| CYP2C19 | 397 | 46.654 | 44.067 | **42.519** |
| CYP2D6 | 396 | 115.590 | **95.172** | 96.064 |
| CYP2E1 | 397 | 1270.449 | **1228.312** | 1228.512 |
| CYP3A4 midazolam | 396 | 849.495 | **703.854** | 707.335 |
| CYP3A4 testosterone | 396 | 3549.908 | **3166.839** | 3262.347 |

Original source numeric MAE values, smaller better. Fixed XGBoost finds a gain from target CYP RNA over age+sex in **all 10**, but the broad multi-CYP panel yields lower MAE than the small CYP-specific RNA arm in only **3/10** activities. Relative to fixed Ridge on the exact same donors/partitions, target RNA XGBoost is better in **5/10** and worse in 5/10. This supports a more robust conclusion that adding features or model complexity is **not consistently superior** on this source. Still no measured same-human CYP protein concentrations and no independent external HLM cohort in this analysis.

XGBoost source-private model code: `src/omics2pk/hlc_real_xgboost_validation.py`; main private runner option `--xgboost`. This comparison was actually executed on the original uploaded files privately, and the GitHub CI never accesses those data.

## Scientific limitations and next required experiment

- **E2 in-vitro** measured CYP probe substrate turnover in microsomes, not a popPK model, human clinical CL/F, midazolam AUC or patient plasma samples. CYP3A4 metabolite formation is not necessarily CYP3A4-only in every HLM due CYP3A5 activity.
- Same full original study donor cohort; no independent external clinical/liver centre, CV over source-processed arrays may have batch/normalization learned across original samples before Omics2PK splits.
- Original dataset has 466 array profiles and 836 phenotype keys, but activity complete-case subsets are n387–397 with age/sex: **do not claim all 466 have observed CYP outcome** or 195 donors have paired genetics.
- AGE and SEX alone are an intentionally simple benchmark; original disease status, genotype uncertainty, induction, comedications and measured protein expression not fully controlled. No population-size-independent inference about multi-omic superiority.
- Source author's genotype sharing terms and participant privacy remain in force; derived aggregate metrics and generic code are the only material permitted in public repository in this project.
- Before claiming that **transcriptomics is sufficient** and multiomics cannot help, locate the ORIGINAL SAME HUMAN quantitative CYP proteins plus RNA and these exact assays; compare donor-heldout and independent external performance on EXACT identical donors. Without proteins, the RNA-versus-RNA+proteomics hypothesis remains **unanswered**.

## Reproduction and confidentiality

Generic source-only reproducible code in:
- `src/omics2pk/hlc_real_cyp_validation.py`
- `src/omics2pk/hlc_real_fullarray_validation.py`
- `src/omics2pk/hlc_real_genomics_validation.py`
- `scripts/benchmark_authorized_hlc_original.py`
- `tests/test_hlc_real_cyp_validation_contract.py` (ONLY fabricated sample matrices)
- `.github/workflows/hlc-actual-cyp-research-contract.yml` (ONLY fabricated sample arrays; **never downloads or stores original original files**).

Local or authorized research-workspace reproduction can run `python -m scripts.benchmark_authorized_hlc_original --expression <local-private-expression.zip> --phenotype <local-private-phenotypes.zip> --genotype <local-private-genotype.zip> --fullarray`; emits aggregate JSON only. The uploaded ZIPs have been analyzed directly in the user conversation's transient execution environment; **NONE has been committed to GitHub, and no individual row/participant ID is reproduced in public artifacts**.

**Status: `ORIGINAL_N400_HLC_RNA_TO_ACTIVITY_E2_QUANTITATIVE_VALIDATED; PAIRED_SOURCE_GENOMICS_AND_WHOLE_RNA_ABLATIONS_EXECUTED; SMALL_VARIABLE_LEVEL_INFORMATION_ADVANTAGE_FAVORS_TARGET_CYP_RNA; ORIGINAL_SAME_DONOR_CYP_PROTEOMIC_LAYER_ABSENT`.**
