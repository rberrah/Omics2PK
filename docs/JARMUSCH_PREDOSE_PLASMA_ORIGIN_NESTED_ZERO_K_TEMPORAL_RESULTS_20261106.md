# Omics2PK decision record — sampling tissue, molecular provenance and 0/1/3-feature choices

**Original-source-qualified update: 2026-11-06.**
**Research-only:** actual same-volunteer CYP3A4/5 PGx, untargeted predose
plasma MS1 and oral midazolam CL/F (Jarmusch 2020, PMID 32216086,
DOI 10.1111/cts.12785).
No claim of validated tacrolimus prediction.

## 1. Which organ actually supplied each biological modality?

| Cohort | Predictor analyzed in Omics2PK | Anatomical specimen source | Functional outcome and source | What it cannot establish |
|---|---|---|---|---|
| Jarmusch | Untargeted MS1, source-annotated MS2, selected at **day1 t0** or **day9 t0** | **Peripheral blood-derived plasma**, not liver/intestine tissue | True measured systemic midazolam oral CL/F, plasma concentration profiles at days 1 and 9 | Hepatic CYP3A abundance or isolated intrinsic CL; genomic DNA collection tissue not verified from article |
| Jarmusch (currently unused) | Urine MS1 | Urine, dynamic urinary excretion | Original urine drug and endogenous metabolite profiles | CYP3A-specific hepatic phenotyping without renal and timing adjustment |
| Jarmusch (currently unused) | Fecal MS1 and stool 16S | Stool/fecal gut lumen community, **not intestinal enterocytes** | Gut microbial composition / fecal metabolome | Human enterocyte CYP3A mRNA or activity |
| Jarmusch (unused) | \`skin\` labeled MS source assays | GNPS metadata contains analytical source label \`skin\`; detailed collection provenance not established by paper here | None in the current model | Liver-specific function |
| HLC / Yang 2010 | SNP (original rs776746 or rs1799853), CYP-target RNA microarray | **Human liver tissues**: liver DNA and liver RNA | Liver microsomal CYP activity from same donor's liver | Patient in-vivo drug AUC, or blood RNA as surrogate for liver |
| Wegler n15 | Corresponding CYP protein or 5-CYP panel | **Human liver**, isolated hepatocytes and microsomes | Measured *in-vitro* substrate CLint in each liver compartment | A noninvasive clinical plasma marker |
| US EPA n100 | CYP qPCR assays | Human liver | CYP activity in microsomal preparation | Original patient PK or paired germline genotype |

**Source gate:** [Original Jarmusch specimen-time audit #42131046154](https://github.com/rberrah/Omics2PK/actions/runs/42131046154).
Original GNPS source has 1,440 technical metadata rows, 1,202 unique
analytical sample names, **13** biological people with day1 *predose plasma*,
**12** with day9 *predose plasma*, and **11** having both. The separate strict
PK complete-case join finds **10 of those 11** with an *observed numerical
day9 midazolam CL/F*. Preserve this effective n10, do not impute the missing
outcome or count analytical rows as patients. Blood was taken pre-dose and
after dosing 5 min, 0.5, 1, 2, 4, 5, 6 and 8 h on both days (Jarmusch 2020).
Cefprozil was taken orally twice daily on days 2–8.

**The two plasma time clocks correspond to different questions:**
- Day1 t0 -> day1 CL/F: contemporaneous pre-first-dose metabolic signal.
- Day1 t0 -> day9 CL/F: **eight-day baseline prognostic marker** before
  antibiotic intervention, not an assay of post-antibiotic CYP3A.
- Day9 t0 -> day9 CL/F: contemporaneous pre-second-dose circulating phenotype
  after antibiotics. Day9 t0 plasma is not day9 liver or intestinal biopsy.
The latter two comparisons are source-qualified on the **same n10 donors**.

## 2. What exactly justifies 0 versus 1 or 3 additional molecular predictors?

**Not an observed number of 'true causal metabolites'.**
The pipeline chooses **model capacity for each independently trained
outer fold**. It is not a rule switching biomarker assays for a
particular heldout patient.

Current reproducible source code:
\`scripts/benchmark_original_jarmusch_nested_pk_germline.py\`:
- Outer **leave one biological person out**; all train-fitted operations
  exclude that person's observed plasma MS1 and CL/F target.
- **4-fold inner model-selection CV**, repeated from scratch inside
  each outer training set. 18 *putative* steroid/bile GNPS library
  matches versus 207 other library candidates, versus full 6349 MS1
  analyte-ion space; source is mixed-sample GNPS project, but predictor
  matrix uses only **plasma t0 rows**.
- k in {**0,1,3**}, Ridge alpha in {**1,20,100**}. The \(k=0\)
  candidate is always available even within each molecular candidate
  group, ensuring the model is permitted to reject all omics.
- Train-only candidate ion prevalence >=50%, spectral feature
  standardization, genotype onehot and supervised ion ranking.
  Current ranking uses **marginal correlation with log(CL/F)**,
  *not residual association conditional on PGx*. This is a
  limitation of biological interpretation and potentially
  suboptimal additional-information selection.
- Select minimum *inner-fold heldout absolute **log(CL/F)** error*.
  If exactly tied, prefer **fewer** omic ions, then **higher**
  regularization alpha. The reported outer MAE is calculated
  in original CL/F numeric units, so inner and outer metrics
  are not exactly the same loss.
- Diagnostic \(\Delta_{\rm inner}=\min_{\alpha,k>0}
  {\rm MAE}_{\log} - \min_{\alpha,k=0}{\rm MAE}_{\log}\).
  Positive -> zero molecular features preferred on training data;
  negative -> one or three molecular predictors preferred on
  training data. Code reports aggregated mean/median/range,
  k counts, and paired heldout PK error without person values.
- k1 or k3 does **not** establish molecule identity, biological
  enzyme specificity, causal mediation or an independently
  reproducible clinical improvement.

### Actual original-sample k justification: day1 oral midazolam n13

| Molecular candidate source | Decision k0 | Decision k1 | Decision k3 | Median inner logMAE (best k>0 minus best k0) | Outer MAE |
|---|---:|---:|---:|---:|---:|
| Clinical weight+sex+PGx | 13 | 0 | 0 | Not applicable | **6.195090** |
| +18 *putative* steroid/bile ions | **13** | 0 | 0 | **+0.016277** | **6.195090** |
| +207 MS2-annotated ions | **13** | 0 | 0 | **+0.016325** | **6.195090** |
| +6349 raw untargeted MS1 ions | 5 | 3 | 5 | **−0.008383** | 6.612286 |

For the two restricted panels, **all 13 independent
outer-training experiments** found that the best nonzero
omics configuration had a worse inner heldout logMAE than
some k0 configuration. There was no biological 'feature absence':
all 18 putative steroid/bile ions are detected in the 13
original day1 predose plasma samples at least once, and all
occur in >=50% of those samples (audit original data).
The selected zero therefore means "insufficient conditional
generalization improvement under *this* estimator and tiny cohort".
For full 6349 MS1, k>0 was preferred in **8/13** inner train
experiments, yet actual outer MAE **worsened**: inner selection
does not guarantee external predictive usefulness.

### Actual original-sample day9 CL/F, n12; molecules still from day1

| Molecular group | k0 | k1 | k3 | Outer MAE |
|---|---:|---:|---:|---:|
| Weight+sex+PGx | 12 | 0 | 0 | 3.310359 |
| 18 putative steroid/bile ions, baseline J1 | 9 | 3 | 0 | **3.252587** |
| 207 annotated source ions, baseline J1 | 7 | 3 | 2 | 3.780407 |
| All 6349 baseline J1 MS1 ions | 4 | 4 | 4 | 3.675309 |

The 18-ion gain on n12 is **only 0.057772 absolute MAE**
(~1.75%), descriptive paired bootstrap interval
[-0.116065,0] (not model-refitted/confirmatory, multiple
arms explored). It is **not** post-antibiotic omic
measurement and cannot be promoted as evidence of an
acute hepatic CYP3A functional biomarker.

## 3. NEW real same-n10 matched timepoint test

Source run: [GitHub Actions #42132075282 (SUCCESS)](https://github.com/rberrah/Omics2PK/actions/runs/42132075282).
Same original 10 people for both predictor clocks; same original observed
day9 CL/F and original weight/sex/CYP3A4/5 PGx on each arm. Both
plasma timepoints precede administration of the oral test probe
*at their respective study visits*.

**Heldout original day9 CL/F MAE:**

| Prediction model, same 10 original people | J1 baseline plasma omics | J9 predose plasma omics |
|---|---:|---:|
| Weight + sex + original CYP3A4/5 PGx | 3.356567 | 3.356567 |
| +18 putative steroid/bile ions | 3.535090 | 3.611280 |
| +207 MS2 source ions | 3.991813 | 3.450288 |
| +6349 untargeted MS1 ions | **3.081850** | 3.705224 |

The baseline J1 6349-ion arm improves numerical original
day9 CL/F MAE by **0.274717 = 8.19%**; however,
it chooses k0 in 6/10, k1 in 2/10, k3 in 2/10,
and the descriptive 95% paired patient loss interval
**[-0.790684,+0.141669] includes zero**.
This is a **post-hoc exploratory comparison after several
candidate panels / sample clocks have been tested**.
No external cohort, correction for multiple hypotheses,
verified MS1 ion identity or stable selection across
sites; no independent positive validation can be claimed.
At **contemporaneous J9 plasma**, every molecular arm is
numerically worse than weight+sex+PGx in the same ten people.
Do **not** interpret this as proof J9 omics intrinsically
harmful or J1 omics intrinsically superior.

## 4. Implications for genotype and sampling design in tacrolimus

The relevant organ depends on the transplanted organ and the PK pathway:
- Kidney transplant recipients: recipient liver and intestinal CYP3A5
  genotypes dominate the *systemic oral PK* comparison, not donor
  kidney CYP3A5 genotype (donor kidney genotype can matter for
  intrarenal processes). Source https://pmc.ncbi.nlm.nih.gov/articles/PMC5846312/
- Liver transplant: donor liver **CYP3A5** genotype and recipient
  intestinal CYP3A5 genotype can affect tacrolimus PK differently.
  Source https://pmc.ncbi.nlm.nih.gov/articles/PMC6709420/
- Direct liver RNA/protein is mechanistically close to hepatic
  enzyme activity, but inaccessible as a routine pre-dose
  blood test. PBMC transcript abundance does NOT automatically
  approximate hepatocyte CYP3A4 expression.
- Peripheral **plasma** candidate endogenous metabolites may
  reflect hepatic enzyme function but are affected by gut
  production, transport, chronic vs acute sampling timing,
  protein binding, and renal/biliary elimination. A strong
  circulating association is NOT a tissue-specific CYP assay.
- **Negative clinical precedent:** 4beta-hydroxycholesterol
  was measured alongside tacrolimus PK in 79 kidney
  transplant recipients (625 marker assays,
  1999 tacrolimus whole-blood concentrations), without
  credible incremental prediction early after transplant
  (PMC5465334); another n279 with full
  CYP3A4*22/CYP3A5 PGx and pretransplant 4beta-OHC
  also found no early dose-prediction value (PMC5651327).
  These are **literature evidence**, not independently
  re-fitted original Omics2PK data. This particularly
  cautions against treating one serum enzyme proxy as a
  general genotype-to-PK solution.

**Scientific approval gate before choosing k>0 in a real clinical
trial:** verified specimen (organ, matrix, tissue,
aliquot processing), verified **sample time before prediction**,
known molecule/assay, explanatory mechanism or
biologically testable pathway, true incremental
**patient-heldout** performance vs a robust
clinical+full PGx comparator, stability and uncertainty,
cost and replicate-site transportability. A model
selecting k0 is a scientifically useful *negative
finding* only within its tested candidate space.
k1 or k3 must **earn a meaningful heldout gain**;
CV inner selection alone is insufficient.

**Development protocol**: the next independently frozen,
preferably larger cohort should compare marginal vs
PGx-residualized ion ranking, a predefined 1-standard-error
parsimony rule, validated targeted metabolites,
actual clinical effects and external center holdout.
Treat any further re-tuning on these n10/n13 source
outcomes as **exploratory model development**.
