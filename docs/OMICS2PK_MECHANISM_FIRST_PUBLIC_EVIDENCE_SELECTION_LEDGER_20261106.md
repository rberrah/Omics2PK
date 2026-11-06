# Omics2PK — Mechanism-first scientific charter, original public evidence and variable-selection ledger
**2026-11-06 · exploratory research charter · original human data only · clinical tacrolimus gain unproven**

## 1. Question & the two VALID ways to answer it

**Mechanistic / hypothesis-driven:** A named variant or physiological process predicts a specific molecular consequence (splicing, abundance, activity, transporter function, signaling), localized to a specific **organ** and a specific **specimen / time**. Define a source-supported restricted cascade and assay in advance. Test: (A) mutation -> intermediate signal, (B) signal -> separately measured enzyme activity/PK, and (C) **incremental** patient-heldout performance of that signal beyond authentic full clinical+PGx. A negative result can reject one specific link without disproving the entire mechanism.

**Discovery / data-driven:** on training-only data, screen broader omic signals (XGB/SWEET, network inference, SHAP/permutation importance, pathways); **SHAP is descriptive model attribution, not a causal mechanism**. Generate mechanistic hypotheses and use independent assays/cohorts to test them. Freeze any new biomarker panel before formal evaluation. Never reuse the original cohort as if it were external validation.

These axes can meet: evidence/prior pathway restrictions -> omic selection -> functional assay -> model AUC/CL/F, with a null/no-omics model always available. Do not make a negative feature-set selection by itself equivalent to a negative biological experiment.

## 2. Variant *mechanism* -> layer and specimen decisions (CYP3A example)

**Important correction:** CYP3A5*3/*3 is **NOT missense**. CYP3A5*3 (rs776746, c.219-237A>G) is intronic, activates a cryptic intron-3 splice acceptor, may yield pseudo-exon insertion, shifted frame, premature stop and degradation of aberrant RNA by **NMD**. CYP3A5*6 (rs10264272) is a *synonymous coding substitution that disrupts splicing*; CYP3A5*7 (rs41303343) is a frameshift. CYP3A4*22 (rs35599367) is intronic and increases an alternative nonfunctional splice transcript with **liver-vs-intestine tissue dependence**.

**Authoritative sources:**
- Joint CYP3A4/5 genotyping consensus https://pmc.ncbi.nlm.nih.gov/articles/PMC10565868/
- CYP3A5*3 splicing / NMD https://pmc.ncbi.nlm.nih.gov/articles/PMC7907328/
- Liver-specific CYP3A4*22 splice phenotype https://pmc.ncbi.nlm.nih.gov/articles/PMC4674354/

| Actual sequence mechanism | Prior causal expectation | PRIORITY experimental layer | Specimen that can test it | What **not** to infer |
|---|---|---|---|---|
| CYP3A5*3 splice / frameshift / PTC | Incorrectly spliced transcript, lower full-length functional CYP3A5 | Junction-/isoform-specific RNA, functional CYP3A5 protein, activity; downstream targeted endogenous substrate-product as indirect proxy | Liver / intestinal CYP3A5-expressing cells for splice and protein; peripheral plasma/urine for indirect circulating functional proxies only | Total bulk RNA count = active protein; CYP3A5*3 = missense; NMD always complete |
| CYP3A5*6 splice (synonymous) | Exon skip and reduced active protein | Isoform-specific RNA and protein/activity | Relevant CYP3A5-expressing tissue | Synonymous sequence = functionally neutral |
| CYP3A5*7 frameshift | Truncated or absent protein | Protein/peptide identity and activity, RNA only if NMD mechanism plausible | Relevant organ tissue | Gene-specific bulk RNA = active full-length protein |
| CYP3A4*22 intronic | Tissue-specific splicing and CYP3A4 expression loss | Alternative-splice RNA ratio + CYP3A4 protein/activity | Liver; intestinal tissue a distinct biological control | A blood/PBMC gene expression level necessarily mirrors hepatic response |
| **Hypothetical true missense** | Active-site change with transcript/protein still present | Protein abundance + function, substrate specificity, direct activity | Human liver/intestine or controlled recombinant system | Normal RNA/protein abundance excludes abnormal activity |
| **Hypothetical late stop** | May escape NMD depending on exon junction context; altered catalytic peptide | Isoform junction RNA + proteotypic peptide + activity | Relevant organ | All stop codons always cause NMD |
| Normal CYP3A5 genotype but lower CYP3A activity | Enzyme function altered through inhibitors, inflammation, expression, disease, competing pathways | Targeted functional metabolites / exogenous probe under controlled conditions | Baseline plasma/urine with renal/liver covariates; source organ for mechanistic confirmation | CYP3A4+3A5 systemic biomarkers are isoform-specific without experimental evidence |

**Context for tacrolimus:** liver + enterocyte CYP3A4/5, product 13-O-demethyltacrolimus; cofactor CYP3A4, ABCB1, inflammation, azoles, hepatic function and hematocrit are biologically credible covariates. Liver microsomal studies quantify variable contributions of CYP3A5 to 13-O-demethylation, so a CYP3A5 SNP only partly explains net CYP3A function: https://pubmed.ncbi.nlm.nih.gov/15951320/. In a **kidney transplant recipient**, recipient liver/intestine matter for systemic oral PK and donor kidney for local graft drug biology; in a **liver transplant**, donor liver and recipient intestine genotype must be separated. **Causal markers upstream of drug dose** and endogenous substrate/product ratios should be distinguished from *postdose drug metabolite concentrations*, which can leak contemporary PK target.

### Mechanism-admissible candidate assays in real fluids

- **Urine:** 6β-hydroxycortisol/cortisol or 6β-hydroxycortisone/cortisone: CYP3A-associated endogenous flux, but **principally hepatic CYP3A4 not direct CYP3A5**, circadian variation, renal clearance, collection interval, dilution and LC-MS/MS identity matter. Pre-dose t0 urine MS1 is only a candidate search; two MS1 ions are not automatically a verified substrate-product ratio.
- **Plasma:** 4β-hydroxycholesterol/cholesterol (slow-response marker, weeks-long timescale, dyslipidemia/statin confounding); CYP3A-associated bile-acid ratios potentially affected by transporter and microbiota. None alone measures CYP3A5 activity or intestinal first-pass. For early tacrolimus dosing, previous kidney transplantation plasma 4β-OHC studies reported no useful extra early predictive gain, so explicitly test the failure hypothesis.
- **Tissue:** splice-specific mRNA + CYP3A4/5 protein + CYP substrate conversion (direct biological linkage). Clinical availability much lower. Whole blood / PBMC RNA should **not** be treated as surrogate hepatic CYP3A by default.
- **Urine versus plasma comparison** is not simply choosing whichever gives smaller MAE: differences in matrix, renal handling, sampling, stability, metabolite detectability, and time must be documented. Same patient × visit PK target and nested patient-disjoint comparisons are required.

Background functional marker review: https://pmc.ncbi.nlm.nih.gov/articles/PMC13688978/; urinary ratio circadian variability: https://pubmed.ncbi.nlm.nih.gov/10805065/.

## 3. What has actually been evaluated on ORIGINAL PUBLIC DATA?

All numerical endpoints are DIFFERENT and cannot be pooled as an across-cohort MAE.
The gene/protein location and sample time must accompany each result.

| Source, n biological donors | Original specimen & collected omics | Original endpoint | Genetic comparison? | Observed original-source result | Selection rationale + limits |
|---|---|---|---|---|---|
| **HLC / Yang**, 171 matched | **Human liver** DNA SNP rs776746, CYP3A4/5 RNA and separately measured liver microsomal CYP3A probe activity | In-vitro CYP3A midazolam probe activity | **YES: one genuine CYP3A5 SNP only** | Fixed repeated CV age/sex+SNP MAE **697.803**; age/sex+RNA **641.842**; age/sex+SNP+RNA **641.836** | RNA CYP3A4/5 biologically preselected. **No full CYP3A PGx and no original nested replay yet** (authorized input ZIPs absent). CYP3A5*3 source strand/function needs confirming. |
| **EPA original human liver**, 93–95 for 8 enzymes | **Human liver** CYP RNA qPCR; microsomal CYP measured activity | 8 in-vitro enzyme activities | **NO** | CYP3A4 clinical age/sex MAE **2080.332** ->+target RNA **1766.969** ->+8 CYP RNAs **1755.430** (Ridge alpha1). Broader 8-CYP worsened 5/8 vs matching CYP RNA at alpha1; strong penalty altered directions. | One target CYP selected a priori vs all 8; cohort contains no germline genotype or patient PK. [Original CI](https://github.com/rberrah/Omics2PK/actions/runs/42101612511). |
| **Wegler liver donor**, 15 | **Human liver hepatocytes and microsomes**, 5 CYP protein measurements, *in-vitro* measured drug CLint | 5 real drug intrinsic metabolic clearances (midazolam, omeprazole, diclofenac, bupropion, bufuralol), each HH/HLM | **NO** | In **hepatocytes** corresponding-CYP Ridge alpha1 wins 5/5 over 5-protein Ridge alpha1; midazolam logCLint MAE **0.162640 vs 0.192036**. HLM results heterogeneous; XGB sometimes wins. | Corresponding substrate-enzyme probe selected a priori, no clinical PK/PGx. Only n15; CYP4/5 shared reaction specificity not proven. [Original CI](https://github.com/rberrah/Omics2PK/actions/runs/42105373917). |
| **Jarmusch**, n13 day1, n12 day9 | **Peripheral blood plasma**, time0 **day1 before oral probe and before antibiotic treatment**; original clinical weight/sex, available CYP3A4/5 genotype; GNPS untargeted MS1 | True in-vivo oral midazolam CL/F on days 1 and 9 | **YES: limited source PGx** | Day1 clinical+PGx MAE **6.195090**; +18 source-putative steroid/bile ions **6.195090**; +207 library ions **6.195090**; +6349 untargeted **6.612286**. Baseline J1 omics -> day9 CL/F n12 baseline **3.310359**, +18 **3.252587** (−1.75%), +207 **3.780407**, +6349 **3.675309**. | k=0/1/3 chosen nested with always-available zero, 4-fold inner Ridge alpha1/20/100, outer person-LOO. 18 annotations only tentative, no strict molecule IDs. Day9 outcome with J1 omics is **baseline prognosis across cefprozil**, not contemporaneous postantibiotic phenotype. [Original CI](https://github.com/rberrah/Omics2PK/actions/runs/42129571169). |
| **Jarmusch same 10 paired**, n10 (of 11 with both J1/J9 plasma, only 10 with true day9 PK) | **Plasma time0 day1 vs plasma time0 day9**, same clinical+PGx, source-verified subject ID join | True midazolam day9 oral CL/F | **YES: limited source PGx** | Baseline **3.356567**; untargeted J1 MS1 **3.081850** (post-hoc gain 8.19%, descriptive CI across zero), untargeted **J9** MS1 **3.705224** (worse). 18 targeted annotations J1 **3.535090**; J9 **3.611280** (worse). | Apparent positive J1 untargeted n10 not robust and exploratory; timepoint critically changes result, likely timing/measurement/finite sample effects. [Original CI](https://github.com/rberrah/Omics2PK/actions/runs/42132075282). |
| **Rotroff**, n32 | Predose **human blood/plasma** GC-TOF metabolomics, 384 spectral features, metformin experiment | Observed human metformin AUC and renal clearance | **NO** | Source reproducible prior study: AUC MAE **2826.047 ->2825.024**, essentially unchanged; renal CL MAE **115.111->118.145** worsened. | Small n, many unidentified spectral peaks, no genotype, renal-dominant different biology; do not call it CYP3A evidence. |
| **EATRIS-Plus** n88 (MultiOmics, supporting negative control) | Same-donor WBC RNA and circulating protein/metabolite omics | Measured cross-omic reconstruction, **not PK** | **NO** | Underregularized addition of views can degrade heldout interomic reconstruction; gene/protein tissue mismatch is mechanistically important. | Not evidence of human drug clearance or gene→PK improvement. |

**Access gates (not negative studies):** Wegler n54 human rosuvastatin paired liver transporters + clinically measured PK was published but patient-level AUC table not yet source-qualified in the public numeric supplement. Silveira n40 predose rosuvastatin raw spectra public but verified complete individual source PK crosswalk unavailable. Shin n24 genotype and urinary endogenomics + IV midazolam clearance published but original 24×period source measurements not yet verified. Never manufacture same-patient joins across these sources.

## 4. Detailed original Jarmusch model selection: WHAT was and was not chosen

**Important:** source library hits describe *putative ions*; without actual validated retention/mass standards these are not confirmed distinct steroids, bile acids, CYP3A reaction products or measured intracellular functions. The algorithm's k counts record number of ion features selected **for each training fold**, not a prescription for a particular heldout patient.

**Day1 CL/F / day1 plasma n13** (outer person LOO, inner 4-fold, penalty fixed candidate grid):

| Eligibility gate for original molecular signals | k0 folds | k1 | k3 | Median difference best k>0 minus k0 inner log error | Heldout CL/F MAE | Why |
|---|---:|---:|---:|---:|---:|---|
| Source clinical weight+sex+PGx | 13 | 0 | 0 | — | **6.195090** | No molecular inputs |
| 18 putative steroid/bile GNPS annotations | **13** | 0 | 0 | **+0.016277** | **6.195090** | k>0 has higher inner heldout error in every training set, **despite all 18 measured** in ≥50% of predose plasma |
| 207 source MS2-annotated ions | **13** | 0 | 0 | **+0.016325** | **6.195090** | Same inner crossvalidated negative margin for extra molecular features |
| All 6349 raw MS1 ions | 5 | 3 | 5 | −0.008383 | 6.612286 | Nonzero k chosen 8/13 based on internal error, but **outer test degrades**: overfit risk |

**Day9 CL/F / predose day1 plasma n12**: 18 putative-target ions k0=9, k1=3, k3=0; no stable positive clinical increment, despite numerical ΔMAE −0.057772. The person n10 same-cohort two-timepoint comparison produces a larger *exploratory* J1 full-MS1 improvement, but **not for J9 MS1 or source annotations**.

**What is NOT YET known:** specific ion feature IDs *frequently* selected across training folds and a confirmed original chemical identity/functional pathway. The existing benchmark reports feature counts and source-annotated candidate categories, not a verified compound-to-enzyme list. Do **not** name “cortisol selected” or “CYP3A4 functional metabolite retained” from generic GNPS steroid keyword matching.

## 5. A formal selection/rejection log is needed — not just model MAE

Every candidate (gene/protein/transcript/ion/metabolite) gets a row:

| Required provenance and evidence | Allowed field |
|---|---|
| \`gene, allele_haplotype, molecular_mechanism\` | Evidence-qualified variant and its splicing/missense/PTC/regulatory mechanism |
| \`source_compartment, organ, tissue, matrix, donor_recipient, collection_time\` | E.g. recipient peripheral plasma pre-dose day1; donor liver protein; spot urine vs timed excretion |
| \`molecule_identity, confidence_level, assay_quality\` | Chemically standard-validated vs putative MS2 vs unknown MS1 ion |
| \`mechanistic_pathway, expected_direction, evidence_reference\` | CYP3A5 splice -> full-length enzyme; CYP3A4 hydroxylation of endogenous substrate; competing CYP/transport routes |
| \`eligible_for_hypothesis, excluded_reason\` | Wrong tissue, wrong time (postdose for current dose prediction), known non-specificity, non-detected, excessive missingness, genomic-effect contradiction, uncertain chemical identity |
| \`n_training_folds_k_selected, selection_stability\` | **Training-only nested** frequency; no heldout leakage |
| \`delta_inner_conditional_PGx_error, delta_outer_PK_error, uncertainty\` | Separate tuning from validation; compare baseline measured clinical+*complete* PGx |
| \`mechanistic_test, independent_replication, confidence\` | Splice-junction assay, targeted protein, validated substrate/product ratio, protein activity, external cohort |
| \`decision\` | Retain as hypothesis, measure targeted, exclude from assay, keep exploratory, or no supported extra biomarker |

**Do not exclude a transcriptome by gene variant class alone:** e.g. CYP3A5*3 strongly motivates *isoform-specific RNA* even when bulk transcript abundance cannot interpret CYP function. A late nonsense variant may escape NMD, so require a transcript-specific test, not a categorical skip.

### Proposed biological candidate gating

1. Check actual measured original **variant consequences**, homozygous/compound heterozygous genotype and strand/haplotype. Detect missing CYP3A4/5 alleles and subject ancestry implications before mechanistic assignment.
2. Draw the pathway DAG **gene -> transcript isoform -> active protein -> enzyme/transporter activity -> endogenous substrate-product -> drug PK**, with context nodes (inhibitors, inflammation, renal excretion, HCT, food, gut microbiome, hepatic function).
3. A priori select **assay layer and compartment** for each testable edge, only when original source has that specimen/time. Keep plasma and urine evidence labels separate, not interchangeable.
4. Run orthogonal **molecular function consistency** tests, not only PK regression: genotype->isoform/protein where measured, direct targeted enzymatic activity, analyte identity/quantification and biologically expected direction. Infer/qualify confounding.
5. For **separate incremental PK question**, require nested selection conditioned on actual clinical+PGx. Prefer residualized predictor screening inside each training fold (not yet implemented in earlier Jarmusch benchmark), include k0, compare n features vs information gain, select with complexity/1-SE parsimony, and report paired **heldout** CL/AUC error with CIs/calibration.
6. Independently freeze and replicate a retained biomarker panel. Separate new hypothesis generation and scientific discovery from confirmatory use of the **same** original people/outcomes.

**Desired scientific message:** one may have highly informative omics about CYP3A biology but still have *no* incremental improvement in finite-sample clinical PK. Conversely, a moderate reproducible PK improvement does not prove the proposed CYP mechanism without orthogonal molecular/function evidence. Both axes are needed.

## 6. Urinary Jarmusch follow-up

Original Jarmusch GNPS metadata already source-qualified true day1/day9 **time=0 urine samples** in 14 original subjects on both days. Study collected multiple post-administration urine samples over 8 h, but **these postdose samples must NOT be used to predict current-dose CL/F** or they would contain downstream drug/metabolite consequences of that CL/F.

New pipeline \`scripts/benchmark_original_jarmusch_urine_t0_nested.py\` and dedicated real original source workflow \`.github/workflows/jarmusch-original-urine-vs-plasma-day9-pk.yml\` compare **urine J1 t0**, **urine J9 t0**, **plasma J1 t0**, **plasma J9 t0** on the same source-qualified original n10 people with both plasma times and actual day9 CL/F; strict matching by original source biological person and analytical MS1 sample ID. Output is aggregate only. No stable quantitatively calibrated cortisol ratio should be claimed from generic MS1 peaks. If the original BIOM does not include all matching urine assay rows, **ABSTAIN**, do not impute or synthetically join.
