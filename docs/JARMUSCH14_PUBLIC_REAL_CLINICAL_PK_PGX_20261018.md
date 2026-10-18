# Real original Jarmusch human CYP PK + germline PGx versus PRIOR functional phenotyping (18 October 2026)

**Original data provenance:** Jarmusch et al. 2020, DOI [10.1111/cts.12785](https://doi.org/10.1111/cts.12785), PMCID [PMC7485946](https://pmc.ncbi.nlm.nih.gov/articles/PMC7485946/), [author's original research repository](https://github.com/alan-jarmusch/UntargetedMS-DrugMetabolism) pinned commit `95586620ead4ac5241dc1e99a98e047c1f99b17b`. Actual per-human clinical target-PK and CYP genetic information was NOT in the separately uploaded Wiley Table S2; it was discovered in the ORIGINAL author `Materials/Metadata.txt`. That public source file has **1423 analytical MS-sample records**, **215 metadata columns**, and only **14 independent biological volunteers**. Original source sha256 `152dd683a1c5904109e8a360a1d67d7631f4be39213c31a8393c652e665fc628`. Actual original rows are retrieved ephemerally in CI only, never saved/committed or printed.

## User-uploaded Wiley Table S2 qualification

User-supplied original `CTS-13-972-s002.xlsx` verified in this private chat; SHA256 `4d84d0ae3a26c6d9f0f46a4804fd87e828d80d9c1a3aabf9d65ff59deafe4ace`, **ONE** worksheet `Serene_FDR_sOTU_annot.txt`, **29 differential 16S sOTUs**, source 7 fields: `16S_sequence`, `Statistic`, `raw pvalue`, `Taxon`, `day_1_mean`, `day_9_mean`, `Enriched_Day`. Observed 24 sOTUs enriched day 1 and 5 day 9, with p-values 0.000999001–0.00999001. These are **group mean relative abundances, NOT 29 independent humans and NOT an original patient×microbiome matrix**. This workbook has NO subject keys, no quantified targeted plasma CL/F, no omeprazole AUC ratio and no person/day matched omic vectors. The original uploaded source file is **NOT committed**.

## Original human study: TRUE systemic drug PK endpoints and time ordering

The author `Metadata.txt` includes original sample-type and actual individual patient study day, sex, body weight, source-genotyped CYP2C19, CYP3A5 and CYP3A4 calls, and independent original clinical noncompartmental PK for same 14 persons over two doses of four drug probes.

- **Day 1**: administered oral caffeine (CYP1A2), omeprazole (CYP2C19), dextromethorphan (CYP2D6), midazolam (CYP3A); original targeted clinical PK was measured *after* day 1 administration, before cefprozil.
- **Days 2–8**: cefprozil 500 mg twice daily.
- **Day 9**: same oral cocktail and new original targeted PK after intervention. Day1 functional probe PK can be a strictly **PRIOR observation** for predicting future day9 PK, but is NOT a non-invasive biological Omics2PK covariate and cannot be used to predict the day1 PK as if measured before the first dose.
- CYP3A oral midazolam `CL/F`: source **14 day1 and 13 day9 quantified**, **13 genuine person-paired** results.
- CYP2C19 source `OME/5-OH OME` AUC ratio: **14 day1 /14 day9**, 14 genuine pairs. **It is a substrate/metabolite AUC ratio, NOT CL/F**.
- CYP1A2 caffeine oral `CL/F`: **14 day1 /13 day9**, 13 genuine pairs. Original data have NO directly assayed CYP1A2 genotype, so we do NOT invent or substitute a CYP genotype arm.
- Source analytical mass-spec rows share repeated values for same individual/day; all original matched copies were checked for consistent pharmacokinetic and genetic original donor identity before collapsing to **one independent person/day PK observation**.
- Published original Table 1 quotes midazolam/caffeine mL/hour, but numerical magnitude is more consistent with L/hour; **original physical clearance unit remains independently unverified**. Report source numerical CL/F scale only until provenance confirmed.

**Important original source metadata error trapped:** field `Enrollment_Subject_Age` literally contains **`yes`** for eligible volunteers, NOT numerical years. We deliberately omitted biological age instead of fabricating or coercing it. Thus source-verified clinical baseline consists of **sex + measured baseline body weight**. The author processed original plasma *untargeted LCMS sample runs* are present for only **13 volunteers on day1 and 12 on day9**, despite the **14 /13/14 original targeted PK measurements** represented in cross-tissue metadata. Never silently assume 14 source plasma MS profiles are available to join all original PK subjects.

## ACTUAL original leave-one-person-out PK comparison, no pseudoreplication

Every model has the **same original complete-paired subjects for that same CYP endpoint**, outer `LeaveOneOut` by person (13 for midazolam/caffeine, 14 for omeprazole). No drug PK from a test human/day1 or day9 is in training of another test-person's response; day1 clinical phenotype is legitimately available to predict later day9 in that person's **test input**. Train-only genotype one-hot encoding `handle_unknown=ignore`, sex/weight and log day1 PK normalization, fixed Ridge α1 fitted to `log(day9 original PK)` and backtransformed. No model/hyperparameter selection by heldout error.

| Original observed day9 NCA outcome | independent paired volunteers | Clinical sex+weight MAE | + authentic CYP genotype | + independently measured prior day1 probe PK | + genotype and prior PK |
|---|---:|---:|---:|---:|---:|
| **CYP3A midazolam oral CL/F** | 13 | 3.76108 | 4.19763 | **3.09465** | 3.47691 |
| **CYP2C19 OME/5-OH OME AUC ratio** | 14 | 0.63293 | **0.60639** | 0.63619 | 0.60059 |
| **CYP1A2 caffeine oral CL/F** | 13 | 16.27584 | Not assayed | 14.84829 | Not assayed |

All MAE in **source original PK numeric scale**, not automatically physically verified original CL unit. The CYP2C19 genotype+prior numerical MAE **0.60059** is slightly lower than genotype alone (0.60639), but other metrics and model complexity matter and n14 is too small for broad claims.

Original **CYP3A**: prior functional phenotype improves relative to genetic baseline (4.19763 → 3.09465, about **26.28% lower MAE**); adding genotype on top of prior PK worsens MAE (3.09465 → 3.47691), while adding the functional phenotype on top of genetic reference improves genotype (4.19763 → 3.47691). Original **CYP2C19**: genotype alone 0.60639 is numerically better than prior functional probe 0.63619; combined 0.60059 is only marginally better than genotype alone. **CYP1A2** cannot answer genotype superiority at all without actual authentic CYP1A2 variants.

### Descriptive uncertainty of person-heldout errors (1999 paired bootstrap resamples, no refitting)

Contrasts are prior phenotype minus genotype, so negative values favor the prior functional phenotype.

- CYP3A midazolam: `clinical+prior PK` versus `clinical+genotype`, ΔMAE = **−1.10299**, *descriptive* 95% **[−2.08315, −0.23599]**.
- CYP3A midazolam: `clinical+genotype+prior PK` versus `clinical+genotype`, ΔMAE = **−0.72072**, descriptive 95% **[−1.43014, −0.00399]**.
- CYP2C19 omeprazole: prior alone vs genotype, ΔMAE = **+0.02980**, descriptive 95% **[−0.14137, +0.18760]**.
- CYP2C19 omeprazole: genotype+prior vs genotype, ΔMAE = **−0.00580**, descriptive 95% **[−0.10836, +0.09016]**.

These intervals only bootstrap the errors from existing LOPO predictions: they do **NOT** refit each model on each bootstrap sample, control for model/endpoint multiplicity, account for external cohort heterogeneity, or prove statistical clinical usefulness with n13. The CYP3A interval being below zero is a descriptive promising finding only, not proof of pharmacogenomics-inferior omics. The original source genotypes are a narrow assayed subset and not independently phased full CYP diplotypes (CYP3A5 *1/*3 and CYP3A4 selected source calls). A prior administered-probe PK parameter incorporates exposure, absorption, CYP3A-mediated first pass and other physiology, rather than being a pure independent enzyme assay.

## GENUINE Omics2PK missing data and what to obtain next

This source is the FIRST FOUND study whose original independent volunteer PK and CYP genotypes are publicly available without a university subscription. To make it TRUE molecular Omics2PK, add genuine measured **PRE-DOSE** metabolomic ion intensity vectors from the **same person's** LC-MS study-day-1 before the probe administration, with sample`filename` matching author metadata. The original author R analysis [Code/20191030_UntargetedMS_DrugMetabolism.R](https://github.com/alan-jarmusch/UntargetedMS-DrugMetabolism/blob/master/Code/20191030_UntargetedMS_DrugMetabolism.R) refers to a processed feature-by-sample input called `GNPS_quant.csv`, **which is not in the author's GitHub repository**. Its raw source is [MassIVE MSV000082493](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=b8ffc66bab8d435aa863a0fed8fff6fa), with potentially corresponding [MSV000084008](https://www.omicsdi.org/dataset/gnps/MSV000084008). The original 29-sOTU Table S2 is GROUP summaries, cannot substitute for donor-level BIOM/16S abundances. Public metadata show only **13 source individuals with plasma untargeted MS on day1 (including T0 predose)** and 12 on day9; the true cross-section can be smaller still after requiring actual PK values. These are not sufficient for high-dimensional omics performance claims without shrinkage and independent replication.

Required benchmark after actual original `GNPS_quant.csv` is retrieved and author sample filename crosswalk source validated:
`sex+weight+CYP PGx` vs `+day1 T0 ONLY baseline metabolite scores` to predict real **day9 PK**, with `+day1 measured probe PK` separately labeled as a reference functional test and **not** as omics. Train-only omics PCA/shrinkage and original complete subject intersection; no choosing omics features using treatment outcome, no comparing different patient sample sets. Biological PK variability after cefprozil may depend on altered gut microbiome, but original group sOTU changes alone cannot be used as a patient-level predictor.

**Do not place any original user-uploaded XLSX, 2.8MB original volunteer metadata, individual MS signals, subject IDs, clinical PK measurements or genotype calls in public Omics2PK GitHub.** Public repo gets code, fake fixture tests, provenance digests and aggregate original cohort metrics only.

**Current science:** `JARMUSCH_REAL_CLINICAL_TARGET_PK_AND_MATCHED_CYP_GENOTYPE_READY_AND_RUN; N13_CYP3A_PREVIOUS_FUNCTIONAL_PROBE_DEMONSTRATION; METABOLOMIC_OR_MICROBIOME_INDIVIDUAL_FEATURE_MATRIX_NOT_QUALIFIED; NO_OMICS_SUPERIORITY_CLAIM`.
