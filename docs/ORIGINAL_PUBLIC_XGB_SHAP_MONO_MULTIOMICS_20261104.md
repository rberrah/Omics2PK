# Original public omics → observed pharmacology: XGBoost, exact Tree SHAP and source-aware white-box attempts
**4 November 2026 · original datasets actually executed. Research exploration, NO claim of independent clinical validation.**

## Conclusions first
1. **Actual human PK E3:** on 32 published metformin-treated humans, predose plasma GC-TOF all 384 peaks **does not** robustly improve source-clinical XGBoost (age, sex, weight, serum creatinine). Unknown 229 outperform the baseline numerically on AUC and renal CLr, but neither descriptive paired confidence interval excludes zero. High SHAP contribution is not verified chemical identity, enzyme function, causal transport or better clinical accuracy.
2. **True 3-omics test, but PD E1:** on 47 original NCI-60 cell lines, 17,987 RNA + 3,171 SWATH proteins + 280 source small molecules, XGBoost 3-omics gives best numerical doxorubicin **GI50** MAE = 0.3230. Relative to RNA-only 0.3336, delta −0.01062 CI [−0.02841,+0.00720], uncertain. Different lab/sample batches, **not PK or human patients**. A prior *Ridge-PCA* 3-omics source replay was **worse** than RNA alone; true conclusion is model dependence, not guaranteed multi-omics benefit.
3. **Partial biology:** AUC high-importance unknown IDs `222169` and `221574`, and named 1-monostearin; CLr unknown `218765` are targets for spectra acquisition. Narrow a priori biochemical panels are not useful in the present n32 PK data. SHAP-top5-only Ridge can fail catastrophically on a human held-out CLr: **safety-abstention**, not a hidden valid model.
4. **New original molecular source:** the official public BinDiscover GitHub contains raw experimental reference **MSP** source files for **15,818 unknown** and **763 curated known** GC-EI entries, independently downloaded in the same CI. An exact-ID screen over all 229 original Rotroff anonymous headers finds **four accession overlaps**, with raw reference spectra present in the public MSP, but none of the high-SHAP three; no source-specific RI/derivatization crosswalk was established. These are **four potential reference accession associations**, zero verified identities.

## 1 — Observed individual human PK (Metformin: Rotroff et al. 2016)
- Primary reference [Rotroff 2016](https://doi.org/10.3389/fphar.2016.00135), original authentic clinical & metabolomic Excel bytes from official publisher bundle, original hashes `76ef6d62e1809f37b7fcf0d9a07e38eb55a67331b23bd058890911547fa9681a` and `1debb35eac7c22383031d4aaeaa47c0362d57609956f16fa4420626efd4a8fe6`. Exact same 33 people original predose A samples, 32 original observed AUClast 0–24 and renal CLr measurements. 384 GC-TOF original columns = 155 source-named chemical candidates and 229 unknown numbered peaks. No RNA/proteomics, PGx or authentic PopPK comparator in this study.
- **Actually run, source replay SUCCESS:** [workflow #41912428806](https://github.com/rberrah/Omics2PK/actions/runs/41912428806), with explicit later safety gate [#41915742531](https://github.com/rberrah/Omics2PK/actions/runs/41915742531).
- **Model:** same source patient heldout (LOO 32), separate AUC/CLr true drug targets, base clinical age/sex/weight and **original screening serum creatinine** (do not include post-dose urine or CLr in feature set), XGB=100 depth2 learning .045 min_child_weight3 L1/2 regularized, model hyperparameters fixed independently. Native **exact Tree SHAP** via `Booster.predict(DMatrix(),pred_contribs=True,approx_contribs=False)`, additivity checked against model's raw log1p(PK) score. Training-only median imputation, training-only SHAP for feature selection. No hyperparameter search using heldout data, no molecular samples after administration.
- **No physical unit labels verified for original drug AUC**, so MAE below in published table numerical units; not `mg*h/L` unless confirmed.

| XGB plus published clinical covariates | Observed drug AUClast 0–24 MAE | Observed renal CLr MAE |
|---|---:|---:|
| Clinical only | **2686.78361** | 114.47331 |
| + all 384 GC-TOF | 2692.81140 | 115.68338 |
| + 155 named candidates | 2834.35069 | 125.84632 |
| + 229 original unknown bins | **2619.88862** | **107.06987** |
| + predefined incomplete urea/purine markers, 5 signals | 2757.05855 | 113.59000 |
| 5 predefined biochemical markers + Ridge white-box | 3124.03261 | 121.94153 |
| Nested train SHAP top5 + Ridge (interpretability attempt) | 3628.65259 | **ABSTAIN: 1 fold unbounded extrapolation** |

**Descriptive** paired 1,999-resample bootstrap of heldout error difference from XGB clinical-only (not a model-refit CI, uncorrected for multiplicity, positive is worse):

| Contrast | AUC ΔMAE [95% descriptive interval] | CLr ΔMAE [95% descriptive interval] |
|---|---|---|
| +384 | +6.0278 [−292.3385,+314.5157] | +1.2101 [−12.9982,+14.3640] |
| +155 named | +147.5671 [−84.2108,+371.9750] | +11.3730 [+0.1658,+23.0431] |
| +229 unknown | **−66.8950 [−374.2792,+254.9591]** | **−7.4034 [−21.9746,+6.9207]** |

**Source public features explaining this *model*, not the causal PK physiology:**
- AUC all384: among metabolite abs SHAP (excluding 4 clinical variables), **86.056% unknown** and 13.944% named; source original numeric bin `222169` mean |SHAP| 0.03867 log-output units selected train-SHAP top5 in **32/32 LOO folds**, `221574` 0.01991 top5 **27/32**, source-named `1-monostearin` 0.00802 top5 **22/32**, then `223505` etc; original source `adenosine` and `hippuric acid` lower.
- CLr all384: unknowns represent **76.084% abs omic SHAP**, top bin `218765` 0.06538 top5 **32/32**; these are not mass values and not transporter protein identifiers.
- These percentages are only relative absolute attribution across molecular columns for a particular XGB model; their base features are more numerous than named. They are **not** variance explained or identification posterior probabilities. Correlated ion signals cause SHAP credit splitting; train-only stability is useful but n32 is tiny.
- **White-box failure:** the small 5-selected-bin Ridge regression (after train-only SHAP, train-only median + scaler) generated an unbounded clinically absurd CLr prediction in 1 heldout fold. Without using heldout truth, a **posthoc safety diagnostic** flags predictions exceeding 2x largest observed training endpoint. Clinical efficacy metrics for this arm are WITHHELD (`ABSTAIN_UNSTABLE_POSTHOC_TRAIN_RANGE_EXTRAPOLATION`), not clipped or re-scored. Its AUC MAE=3628.7 also underperforms baseline. Simply selecting SHAP highest peaks and fitting a linear model is not a successful functional phenotype. More robust white-box should prespecify ridge regularization, biology-informed units and stability criteria on a larger independent cohort.

## 2 — Actual three-layer original cell-line source (NCI60; **not PK**)
Original mixOmics NCI60 doxorubicin GI50; IntLIM 17,987 RNA genes /280 metabolite signals; Guo et al. 2019 [SWATH protein source](https://doi.org/10.1016/j.isci.2019.10.059), 3,171 protein accessions, exactly 47 matched cell lines, different originating laboratories/cell passage. **Original source XGB 3-repeats × 5 independent line-disjoint CV, 141 correlated scores across 47 biological lines**, train-only unsupervised variance selects RNA 300, protein 300, metabolite 120, exact SHAP on true heldout test predictions. Source real **SUCCESS** [workflow #41915019343](https://github.com/rberrah/Omics2PK/actions/runs/41915019343).

| XGB + tissue-lineage class | Heldout original GI50 score MAE |
|---|---:|
| Type tumour only | 0.349593 |
| RNA only | 0.333638 |
| Protein only | 0.338472 |
| Metabolites only | 0.327104 |
| RNA + protein | 0.328180 |
| RNA + metabolites | 0.328804 |
| Protein + metabolites | 0.330808 |
| **RNA+protein+metabolites** | **0.323019** |
| Predefined TOP2/ABCB1/ABCC/ABCG/CBR candidate gene RNA, Ridge | 0.365349 |
| Predefined candidate gene RNA+proteins, Ridge | 0.377913 |

Descriptive paired threeomics−RNA ΔMAE **−0.010620** [**−0.028414,+0.007196**], 29/47 cell lines showing improved individual mean-heldout error. Threeomics−metabolites **−0.004086** [−0.040780,+0.032173], 23/47 cell lines improve. All comparisons include zero (without multiple comparison adjustment). **This is not a demonstrated clinical or causal advantage** and absolute GI50 MAE is unrelated to metformin PK MAE.

Threeomics absolute test SHAP attribution within molecular layers (NOT explained variance, dependent on input representation and number of features):
- RNA **52.967%**; source genes `TMEM139`, `TNFRSF21`, `GOLT1A`, `ATP8B1` and `CLEC2B` highest in this fit.
- Proteome **35.573%**; source gene-symbol labels `IL18`, `PLXNB2`, `LGALS3`, `BAIAP2L1`, `JUP`; gene symbol is not a specific peptide/active protein/PTM.
- Metabolome **11.459%**; source unidentified `X-4046`, source candidate `uracil`, `3-phospho-d-glycerate`, `N-acetylneuraminic acid`, `uridine`.
- Predeclared **doxorubicin** TOP2A/TOP2B, ABCB1, ABCC1, ABCG2, CBR1, CBR3 targeted white-box RNA or RNA/proteins was worse in this n47 fixed Ridge comparison; this does NOT show that these pathways are not involved. Their **actual activity was not measured**.

**Crucial model sensitivity:** original prior fixed PCA+Ridge benchmark with the same 47 source lines had cancer type+RNA MAE **0.334666** and cancer type+3omics **0.358124** (worse), versus XGB 0.333638→0.323019 (better). A graph/multiomic signature must be validated by different methods and on independent tissue batches, not selected because one model score looks favorable.

## 3 — Exact study bin IDs vs publicly accessible spectra: a new source milestone
The [BinDiscover source repository](https://github.com/metabolomics-us/bindiscover/tree/045e7af0861ca04b97be902f4f557e8abb2d4940/api/newer_datasets) contains actual GC-EI MSP:
- `GCBinbase_unknowns.msp` original **39,902,584 bytes, 15,818 records**, source hash `986c73af77f64763bf88ba11dad65a4754f3b8177b9af32673fb075ba92374f4`.
- `GCBinbase_knowns_curated.msp` original **1,870,093 bytes, 763 records**, source hash `329a63da6e58644ba0fb98e11482e9ce6b0e89638fc4f700faf03d80beb2c349`.
- **Source-executed success [workflow #41915819310](https://github.com/rberrah/Omics2PK/actions/runs/41915819310)**: exact original unknown bin IDs 229 aligned against **identity header fields** `BinId` and `Name` (not against random fragment m/z or RI). **4 anonymous Rotroff bins share reference unknown MSP IDs**: `210286` (105 declared fragment ions), `210403` (11), `213182` (92), `356985` (181). Source instrument **Leco Pegasus IV**, mode positive and original library RI respectively **661999**, **463407**, **298062**, **229302** in source RI numerical scale; no cross-study same-method calibration proven. **225** original anonymous bins do not show exact header ID matches in these particular public files.
- The top SHAP bins `222169` (AUC), `221574` (AUC) and `218765` (CLr) have **no exact ID** in the downloaded pinned public MSP headers. The historical public BinDiscover API returned no usable spectra for these three ([feasibility run #41913618382](https://github.com/rberrah/Omics2PK/actions/runs/41913618382)); an API failure is **not** proof no data elsewhere.
- The source numbers are **bin IDs NOT molecular mass or specific PTM**. A matching BinId across sources is a hypothesis for a common library entry, not an experimentally validated exact spectrum/derivatization/crosswalk, nor a transported substrate. **0 authentic new chemical structure/transport relations are assigned**. No raw spectra copied to repo, no patient-level source data saved. For the four candidate overlaps, the next gate is actual deconvoluted **study Rotroff EI spectra/RI** from original lab or vendor to compare with public MSP, not just same numeric ID.

## 4 — Reasonable next research steps without a perfect cohort

1. **Disentangle chemical identification from feature predictive value**: Preserve all 384 original GC-TOF features; no arbitrary removal of unknowns. Validate the BinBase accession namespace, documented derivatization, retention standard and original source link for four overlaps. Ask data producer for EI spectra/RI of top SHAP unknowns, starting `222169`, `221574`, `218765`. Once raw study ions exist, do true matched-ion EI cosine + retention tolerance (and standard/MSI certification), not a database-name guess.
2. **Mechanism:** For source-identified structures, use ChEBI, Rhea, UniProt and curated OCT2/MATE/OAT and CYP/P-gp literature, with tissue, mechanism evidence and negative controls. P-gp is a legitimate **tacrolimus** transporter hypothesis, not automatically metformin. Differentiate abundance of HSA protein/proteoform, PTM and **experimentally measured** free fraction/transport/clearance.
3. **Replicate clinical PK:** Acquire rosuvastatin MTBLS3128 individual predose spectra + original AUC/Cmax targets or target drug tacrolimus clinical PGx+PK; only then test RNA vs proteomics vs metabolites vs all on **same patients**. Where multi-omics absent, evaluate each modality on a qualified matched donor in vitro, not by merging unrelated donors into fictional patients.
4. **Improve whitebox:** choose mechanisms pre-specified by external reactions/transporters and enforce stability/credible prior contrasts in nested training; a numerical SHAP top5 is not automatically a biochemically identifiable score.

**Code**:
- `src/omics2pk/rotroff33_xgb_shap_whitebox.py` + `scripts/validate_rotroff33_xgb_shap_whitebox.py`, CI `.github/workflows/rotroff33-real-pk-xgb-shap.yml`.
- `scripts/benchmark_original_nci60_threeomics_xgb_shap.py`, CI `.github/workflows/nci60-original-threeomics-xgb-shap.yml`.
- `scripts/audit_rotroff_vs_bindiscover_public_msp.py`, CI `.github/workflows/rotroff-bindiscover-public-msp-crosswalk.yml`, immutable public MSP sources downloaded ephemerally.
- Source authentication and n32/47 biological identity guard, train-only steps and **no raw source data committed**.
