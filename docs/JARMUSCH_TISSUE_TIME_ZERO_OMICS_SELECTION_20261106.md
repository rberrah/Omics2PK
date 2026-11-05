# Source tissue, assay location, and defensible selection of 0/1/3 variables
_2026-11-06 — Omics2PK research-only, no proven tacrolimus PK improvement_

## Key principle: "omics" is not a biological compartment

The **biological source organ and the sampling specimen are different fields**.
Plasma ions have no inherent identity as hepatic enzyme activity;
multiple organs can produce, metabolize, transport and clear them.
Same-person AND time-resolved AND specimen-qualified associations are necessary.

### Jarmusch 2020 actual original study ([DOI](https://doi.org/10.1111/cts.12785))

- **Original trial:** 14 recruited healthy adults. Oral Cooperstown probes including
  midazolam at days **1 and 9**, with cefprozil 500 mg BID on days **2–8**.
- **Original source specimen:** venous blood sampling -> **blood plasma**,
  collected before probe dosing (t0) and repeatedly until hour 8 on
  days 1 and 9. Predose plasma LC-MS1 is source true
  \`ATTRIBUTE_Sample_Type=blood_plasma\`,
  \`ATTRIBUTE_Time_Point_Mins=0\`, \`MS_SampleType=sample\`;
  this is not liver protein, liver RNA or a hepatocyte enzyme assay.
- **Other collection sites:** urine (urinary metabolites), feces (gut
  metabolites, stool 16S microbiota); source GNPS job metadata also contains
  \`skin\` analytical samples, which are **excluded** from every plasma benchmark.
- **Genomic locus source:** trial genotyped CYP3A4/5 by PCR/TaqMan,
  but original article methods consulted do **not** establish the DNA
  specimen's anatomical origin; do not invent it as blood DNA.
- **Original molecular matrix:** GNPS 6,349 MS1 peak features from
  the whole study, annotation by approximate spectral MS2 library
  matches; many unannotated, isobaric or nonspecific, no analytic-standard
  confirmation of the source 18 tentative steroid/bile candidate ions.
- **Original source-qualified sampling audit:** [GH Actions qualified run
  42131046154](https://github.com/rberrah/Omics2PK/actions/runs/42131046154).
  Exactly **13 unique** predose day1 plasma people; **12 unique** predose
  day9 plasma people; **11 same people with BOTH** predose plasma assays.
  The 18 putatively annotated steroid/bile ions (selected across original
  mixed specimens) are all detected in the original predose day1 plasma
  and all are present in >=50% of day1 donors, but their source biochemical
  identities and CYP3A specificity remain unvalidated.
- **PK endpoint:** the article measures oral midazolam CL/F by plasma
  concentration-time exposure, reflecting combined gut and hepatic
  first-pass effects and systemic elimination. A genuine clinical drug
  PK phenotype, **not** a direct biochemical measurement of
  isolated hepatic CYP3A4 or CYP3A5 activity.

**Critical temporal distinction in prior source replay:**
Earlier day9 CL/F was predicted using predose **day1 plasma**, taken
before the week of antibiotics. This is a real **baseline prognostic**
prediction over an 8-day treatment perturbation, **NOT** contemporaneous
post-antibiotic metabolic phenotyping and cannot be interpreted as
representing day9 enzyme function. A separate same-11-patient
day1-versus-day9 source-clock experiment is now implemented in
\`scripts/benchmark_jarmusch_temporal_plasma_nested.py\`, with day9
pre-dose plasma as an *independent molecular predictor design* and
the same measured original day9 CL/F outcome.

### Human Liver Cohort (Yang et al 2010, DOI 10.1101/gr.103341.109)

Human liver tissue 1–2 g from Vanderbilt, Pittsburgh and Merck
collections; original study extracts **DNA** (genotyping),
**RNA** (liver transcript microarray), and **human liver microsomes**
for original CYP probe activities **from liver samples**.
Thus liver tissue matches the enzymatic source organ. Donors include
postmortem samples and surgical resections, limiting clinical
transplant transferability. Site/batch differs between collections.
For our genuine original 171–172 donor genotype/RNA/CYP overlap,
the only source-verified variants in the two original SNP contrasts
are CYP3A5 rs776746 and CYP2C9 rs1799853. **Not full CYP3A genotyping**.
Outcome is hepatic *in-vitro* enzyme probe activity, **not**
patient midazolam/tacrolimus AUC or CL/F.

### Wegler donor-matched 15-liver E2

Matched isolated hepatocytes and liver microsomes with targeted
CYP protein abundances, plus *in vitro* probe drug intrinsic
clearance; same liver donor but different intracellular assay
system. This directly interrogates CYP amount versus in-vitro
drug metabolism, but samples are not predose plasma and endpoints
are not systemic patient PK.

## What "selecting zero omics" **actually means**

In \`scripts/benchmark_original_jarmusch_nested_pk_germline.py\`:
1. For each biological person independently held out (LOPO outer
   validation), fix training patients only. Never use held person's
   plasma peak intensities to select feature identity or to tune the
   predictive model.
2. Within that outer training set, perform **4-fold inner CV**.
   On each inner training fold, use source plasma peak prevalence
   >=50%, train-only peak transforms, and candidate sets constrained
   either to 18 putative steroid/bile GNPS annotations, 207 generic
   MS2 library matches, or all 6,349 MS1 peak signals.
3. For each arm, compare **k=0, 1 or 3** variables and fixed Ridge
   penalty \`alpha ∈ {1,20,100}\`. For k>0, features are ranked on
   **inner training labels only** by centered spectral covariance
   with log midazolam CL/F; categorical CYP3A4/5 genotypes are
   encoded using the inner training people.
4. Choose \`(alpha,k)\` minimizing **inner validation mean absolute
   log(CLF) error**. Exact score ties prefer fewer features, and
   stronger penalty. The **outer primary MAE** is on the original
   numeric CL/F scale: optimization and evaluation scales DIFFER
   intentionally. A no-omic arm optimizes alpha with k fixed at 0.
5. k=0 is a valid result indicating that **no measured molecular
   feature subset won this particular small-sample inner
   generalization contest**. It does NOT prove biological
   irrelevance, absence of PK association, or absence of CYP
   activity signal. k=1/3 means omics improved this inner
   training/CV risk, not that the molecule causes CL/F.
6. Each outer fold independently chooses k and alpha. The counts
   (e.g. day1 k0=13/13 for annotated molecules) describe
   **13 independently trained models**, NOT patients individually
   assigned different assay strategies. A true deployable
   frozen model requires a training-only final model and an
   external cohort.

**Quantitative audit added:** for each source-verified original
outer experiment, report both best inner k0 score and best inner
k>0 candidate score **as an aggregate margin**:
\`Δ_inner = MAE_log(best k>0) - MAE_log(best k=0)\`.
When Δ_inner >= 0, k0 is preferred; when < 0, the
nested CV prefers at least one molecular feature. The
report includes counts, mean, median, minimum and maximum
of Δ_inner across heldout donors **without exposing subject rows**.
Outer paired ΔMAE and patient-resampled uncertainty then determine
whether that inner training choice improved heldout CL/F.

### Scientific methodological limitation / next protocol amendment

Current supervised feature ranking is **marginal covariance with
measured CL/F**, not a true residualized partial association
conditional on clinical+PGx. It is therefore a statistical
feature screen, **not a mechanism-proven CYP3A functional
marker** or established molecular contribution beyond genotype.
The head-to-head nested heldout arm *does* test incremental
predictive performance against PGx, but individual ion
selection cannot be interpreted as genotype-independent
causation. Proposed follow-up: train-only residual screening,
repeat-measure paired day9 plasma models, external source and
orthogonal targeted metabolite assay validation. All exploratory
changes after seeing old fold results must be labelled such,
not retroactively pre-registered.

### Why blood does not automatically report hepatic CYP3A5

Plasma metabolites represent production by multiple tissues,
microbiota, protein binding, renal/biliary elimination,
inflammation and substrate availability. CYP3A4/5 are
principally relevant in **liver and intestinal mucosa**,
whereas the measured plasma MS1 peaks could reflect numerous
other pathways. Genotype is stable germline context; plasma
profile can fluctuate, and the plasma sample time is crucial.
A functional PK bridge needs an independently validated
circulating substrate/product measure tied to hepatic/intestinal
enzyme expression/activity (or a direct enzyme probe when
ethically justified), not arbitrary MS1 peaks.

**Primary cited literature:**
- Jarmusch et al. 2020, PMC7485946, DOI 10.1111/cts.12785.
- Yang et al. 2010, PMC2909567, DOI 10.1101/gr.103341.109.
- Original GNPS qualified tasks and code under PR #141.
