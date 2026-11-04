# Clinical real PK vs tissue-qualified predose metabolomics: named, unknown and pathway panels
**2026-11-04 · SOURCE EXECUTED · exploratory E3, not a new clinical validation.**
**Reproducible original data action: [Rotroff n32 real clinical pathway LOOCV #41902687694](https://github.com/rberrah/Omics2PK/actions/runs/41902687694) — SUCCESS.**
**Original source chemistry-label gate: [Rotroff feature-annotation #41902196861](https://github.com/rberrah/Omics2PK/actions/runs/41902196861) — SUCCESS.**

### Original experimental source and biological location
- **Rotroff DM et al. 2016**, *Front Pharmacol* DOI [10.3389/fphar.2016.00135](https://doi.org/10.3389/fphar.2016.00135), [PMC4906013](https://pmc.ncbi.nlm.nih.gov/articles/PMC4906013/). Plasma GC-TOF taken from 33 **human volunteers before first metformin administration** (original sample A). Post-exposure B/C samples excluded.
- Source original publisher workbook binary checksums required: clinical `76ef6d62e1809f37b7fcf0d9a07e38eb55a67331b23bd058890911547fa9681a`, metabolomics `1debb35eac7c22383031d4aaeaa47c0362d57609956f16fa4420626efd4a8fe6`. Source IDs verified 33:33; **32 original humans with actual PK AUC or renal clearance measurements**. Same origin & drug, not matching different studies by demographic profile.
- The original measured matrix explicitly includes **384 independent GC-TOF features, 155 source chemically named and 229 original numbered/unknown feature labels**, as confirmed independently by first-row header audit, exactly matching the paper's 155 structurally identified/384 quantified report. The source XLSX names alone do **not** establish authentic standard/MSI level 1 identity for each compound.
- **Not liver expression or organ function:** This is cell-free human circulating **plasma metabolomics**. No source genotype or actual PopPK structural estimate fitted here. Mechanistic metformin renal secretion by transporters is not directly assayed via the plasma metabolome.

### Seven frozen analysis arms on identical original subject folds

**Outcome**: separately observed metformin **drug AUClast 0–24** and measured **renal clearance CLr**, not glucose response/PD. Both in original clinical table numeric units, whose physical units have not been verified from primary source. Lower MAE is better.

**Method**: internal leave-one-human-out n=32, exact same heldout original patients in every arm; clinical baseline age + sex + weight standardized **training fold**; Ridge(alpha=1) on log1p PK outcome; molecular median imputer, standardization and PCA **training fold only**. All 384, named 155 and unknown 229 each use 3 PCs; the two path panels and their union use 1 PC. No drug outcome-guided metabolite selection or use of B/C postdose samples. **Original previously validated baseline and all-384 results reproduced numerically** to within a 0.02 tolerance (in practice 1e-4 or better).

| Predictors | n original predose MS features | Measured metformin AUC MAE | Measured metformin CLr MAE |
|---|---:|---:|---:|
| **Clinical reference only** | 0 | **2826.047211** | **115.111320** |
| Clinical + all source signals | 384 | 2825.023966 | 118.144925 |
| Clinical + original chemically named | 155 | 3015.952253 | 123.891019 |
| Clinical + original unknown/numbered | 229 | **2709.465855** | 116.999330 |
| Clinical + partial urea-cycle | 3 | 2893.203810 | 119.551999 |
| Clinical + partial purine degradation | 2 | 3074.020534 | 117.471475 |
| Clinical + union of those partial panels | 5 | 2907.014025 | 118.647478 |

**Exact source pathway coverage** (partial, not verified Reactome reaction graph): urea group = citrulline, ornithine, urea (no labeled arginine); purine degradation = inosine and hypoxanthine (no labeled xanthine). Uridine is a pyrimidine and **not** put in the purine-only set.

Paired descriptive bootstrap 95% intervals on **already heldout** absolute-error differences (n32 human units, 1999 resamples; **not** a refit/bootstrap confidence interval and no familywise multiplicity correction; positive delta is worse than clinical-only):

| Candidate minus clinical reference | ΔMAE AUC [CI] | ΔMAE renal CLr [CI] |
|---|---|---|
| All 384 | −1.023 [−592.811, +681.066] | +3.034 [−11.969, +21.604] |
| Named 155 | +189.905 [−134.597, +565.950] | +8.780 [−3.959, +22.845] |
| Unknown 229 | **−116.581 [−621.830, +378.011]** | +1.888 [−9.880, +15.966] |
| Urea panel | +67.157 [−104.422, +276.031] | +4.441 [+2.471, +6.895] |
| Purines panel | +247.973 [+65.033, +474.097] | +2.360 [−1.007, +5.950] |
| Union | +80.967 [−22.321, +203.223] | +3.536 [+1.325, +6.229] |

**Clinically meaningful checks:** Original AUC relative accuracy within 20% was 56.25% for clinical-only, 50.0% for full/named, and **46.875% for unknown-only** despite unknown's numerically lowest aggregate MAE. Thus it is **not valid** to assert that unknown metabolites clinically improve AUC. None improves both PK outcomes reproducibly.

### Interpretation in biological cascade terms

The negative result does **not** prove that blood metabolomics is useless or that genotype is sufficient. The source panels were selected as partial biochemical reaction-neighborhood proxies based on original author article mechanisms and are **not complete pathway activity scores**. The original article focused significantly on **post-dose response and glucose PD**, not a prospective PK mechanism. Our test is exploratory reuse of those hypotheses, not independent a priori pathway validation. A plasma metabolite concentration may reflect production, uptake, diet, gut biology, renal clearance or compartment transport, rather than specific target tissue CYP activity or drug clearance. Renal clearance of metformin in this cohort is an important target but the current baseline lacks independently source-verified predose kidney function and individual real pharmacogenetics.

**Appropriate short result for Pierre Marquet:** in the first source-matched, true *human* PK dataset (32 subjects), adding untargeted plasma metabolites, chemically named features or small incomplete urea/purine panels did **not** clearly improve prediction over a basic clinical reference. Some unassigned features produced numerically lower AUC MAE but also worse ±20% accuracy, with a descriptive interval spanning no benefit. The data are too small and do not include the stronger clinical + PGx + PopPK reference to answer the thesis's central hypothesis.

### Next independent steps with no wet-lab budget

1. Audit original predose kidney function + creatinine/eGFR metadata and fit a **stronger genuine source covariate baseline** before any claims (do not use observed metformin CLr to predict its own AUC as baseline).
2. Obtain original rosuvastatin n40 AUC/Cmax and linked predose intensities (MTBLS3128; [source access verdict](ROSUVASTATIN40_GITHUB_VS_EBI_FIVE_PROCESSED_CSV_VERDICT_20261018.md)). No invented matching or raw spectral-to-PK workaround.
3. Pair with independent tissue-mechanistic human liver CYP activity from HLC/EPA and GTEx liver-vs-whole-blood detection floor, without implying blood enzyme functionality.
4. For future source, separately qualify measured **actual genotype**, cellular specimen origin, source analytical assay confidence, functional tissue and **real PopPK baseline** before comparing omic modalities.

### Reproducibility / privacy
- `src/omics2pk/rotroff33_pathway_e3.py`: original name-to-numeric feature alignment verified one-to-one against original predose values, train-only LOOCV fixed methods.
- `scripts/validate_rotroff33_pathway_e3.py`: downloads original source once ephemerally and prints only aggregate metrics.
- `tests/test_rotroff33_pathway_e3.py`: synthetic-only negative schema tests.
- `.github/workflows/rotroff33-original-pathway-pk.yml` validated actual published source and original clinical baseline replay. **No individual source identifiers, profiles, concentrations, drug PK outcomes or predictions uploaded.**
