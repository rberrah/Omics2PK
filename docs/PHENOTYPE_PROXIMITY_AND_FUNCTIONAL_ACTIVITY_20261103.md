# Phenotype-proximity and CYP3A5 activity ablation — research v0.2 (2026-11-03)

## Revised central question

**On the same people with authentic CYP3A4/5 variants, matched omics, and a genuine independently assayed phenotype, does external biological knowledge improve out-of-sample enzyme activity prediction over genotype alone and direct statistical omic fusion? Does that improve actual measured drug PK on another, correctly linked benchmark?**

The two questions (function and PK) are **separate endpoints**. A positive enzyme-activity result does not establish a positive in-vivo tacrolimus CL/F result. Clinical+PGx remains the main clinical comparison; any intermediate-function model is exploratory.

## Biological orientation — proximal is not specific

Genotype defines a *molecular predisposition*, not patient-specific percentage enzyme activity. RNA measures transcripts and isoforms; targeted protein measures abundance/identity, not automatically catalysis. An independently measured enzyme-specific reaction is the closest experimental CYP3A5 activity reference. Metabolomic abundances integrate formation, utilization, transport, distribution and excretion: they can be **close to the realized global metabolic phenotype**, but an individual metabolite is NOT a CYP3A5-specific readout without independent substrate/route qualification.

4β-hydroxycholesterol or its cholesterol-normalized ratio is a candidate CYP3A biomarker, particularly for sustained enzyme induction, **not an exact CYP3A5 activity measurement**. Long kinetics, CYP3A4/CYP3A5 overlap and liver/intestine differences preclude immediate or isoform-specific interpretation. Similarly, drug metabolites measured **after administration or after the target PK landmark** cannot enter a *predose* predictive model as covariates; they may be a separately defined postdose outcome.

Interpretation hierarchy for **CYP3A5-specific** prediction:
1. Complete original diplotype *1/*3/*6/*7 (or a properly curated expanded panel): categorical expected function; no arbitrary 0/75/95% constants.
2. Full-length CYP3A5 proteomics or validated canonical splice-junction RNA: evidence for an upstream functional mechanism, with independent assay specificity.
3. Independently validated CYP3A5-specific metabolite/tracer: potentially proximal functional proxy; shared CYP3A4/CYP3A5 metabolite remains *exploratory*, not an enzyme-specific curated predictor.
4. Real assayed CYP3A5-specific catalytic rate: the supervised E2 target. Shared CYP3A probe activity and CYP3A5-specific activity are NOT interchangeable.
5. Original tacrolimus dose-concentration/AUC/CL/F: distinct E3 endpoint where dose/formulation, DDI, hematocrit, tissue and timing must be handled.

No omic layer is necessarily sufficient by itself, or guaranteed to provide incremental gain. The *measured function* is the target; omics can act as explanatory intermediates and predictive state proxies.

## Implemented software: strict held-donor activity ablation

`src/omics2pk/functional_activity_benchmark.py` adds `benchmark_cyp3a5_activity` with a **single pre-specified fixed Ridge(alpha=1)** model on log measured specific CYP3A5 activity. Uses 3–10 outer donor-heldout folds; **all model scaling and response fitting happen inside training folds**. Original donor IDs, exact per-assay alignment, positive measured activity in source units, explicit genotype completeness, validated outcome declaration, timestamp and hepatic tissue are required. Omic features must use individual-preserving preparation or externally frozen scaling; arbitrary whole-cohort correction is disallowed. The code never transfers individual donor rows or predicts clinical doses.

Comparison arms on the *same* eligible participants:
- authentic CYP3A5 genotype **alone**;
- genotype + each available RNA, protein, or metabolite layer separately;
- genotype + all measured omics;
- genotype + subset tagged by **externally prespecified** mechanistic annotation.

The curated arm accepts only declared validated canonical splice RNA, full-length CYP3A5 protein, or independently validated CYP3A5-specific metabolic proxy, and requires an external evidence string for included features. Total-gene RNA and shared CYP3A4/5 metabolites remain available in the **uncurated comparator** but cannot silently become qualified mechanistic evidence. A string and a validation flag are **declarations**: the implementation does not independently authenticate DOIs, clinical collection clocks or original source rows, and its reported research source SHA is not cross-checked to raw files.

This is an *association-based, regularized exploratory comparator*, **not a calibrated Bayesian model**. Calibrated activity intervals, quantitative allele penetrance and mechanistic pathway inference have **not** been implemented.

### Executed synthetic qualification (60 artificial donors)

All biological states, functional signals, assay times and outcome values were explicitly generated. There are **no real human omics or PK observations in these tests**.

| Crossvalidated arm | Synthetic MAE, arbitrary simulated activity units |
|---|---:|
| Genotype only | 0.45188548 |
| Genotype + transcriptomics | 0.11955774 |
| Genotype + proteomics | 0.14072428 |
| Genotype + metabolomics | 0.09450018 |
| Genotype + all omics | 0.09883206 |
| Genotype + curated omics | 0.09327349 |

These numbers are **deliberate software test fixtures** constructed with an informative metabolic signal. They **must not be interpreted as evidence that metabolomics is better in real humans, or that knowledge improves PK**. In particular, the synthetic response was designed to correlate with the phenotype-like metabolic proxy. Omics may worsen real heldout performance, as already observed on genuine CYP3A-probe clinical cohorts.

Tested software guards: reversed row-order invariance with *all* original donor arrays aligned, donor-ID mismatch rejection, postoutcome sampling rejection, incomplete genotype abstention, unvalidated activity rejection, missing external evidence rejection, exclusion of CYP3A4/CYP3A5-shared metabolite from CYP3A5-specific curated arm, non-frozen upstream preprocessing rejection, finite-value checks and aggregate-only outputs.

## Data readiness and next scientific experiment

Current genuine Omics2PK examples **cannot run the central four-way causal/functional comparison**:
- EPA human liver RNA + shared CYP activities: original RNA and measured activity exist, but **no same-donor complete CYP3A5 PGx + specific metabolome** is established.
- Latham 21-donor HLM studies: targeted CYP proteins + functional CYP probe activities + drug-specific CLint are present; there is **no matched genotype, RNA or metabolomics** in the qualified tables.
- Jarmusch n13 CYP3A probe with predose MS1 and authentic PGx: clinical midazolam PK is present, but no independently measured CYP3A5-specific activity; raw MS1 peaks are not verified metabolites or enzyme proxies.
- BIOMARGIN–ABIS/tacrolimus: original valid same-patient crosswalk remains missing.

Thus, run E2 only on independently assayed enzyme activity with matched people, then E3 on observed PK without target leakage. Freeze feature annotations using CPIC, PharmVar, verified metabolite identity, published reaction/isoenzyme evidence **before** testing heldout outcomes. Explicitly test knowledge-ablation, shared-metabolite negative controls, RNA isoform ambiguity and partial panel coverage. Require paired heldout intervals and external-site validation before making an improvement claim.

**Scope:** Omics2PK only. MultiOmics frozen research core and `DOWNSTREAM_EXCHANGE_V1` unchanged; optional annotated exchange sidecar is future work.

### Supporting scientific references

- CPIC CYP3A5–tacrolimus, https://pmc.ncbi.nlm.nih.gov/articles/PMC4481158/
- 4β-hydroxycholesterol critical review (2019), https://doi.org/10.1002/jcph.1391
- Recent review (April 2026), https://doi.org/10.1016/j.dmpk.2026.101524
- PharmVar CYP3A5, https://www.pharmvar.org/gene/cyp3a5
