# Public-only Omics2PK: actual omic collection site ≠ target organ

**2026-11-05. Scope frozen to original public datasets or already
previously replayed public-source evidence.** No imaginary future
BIOMARGIN/ABIS or tacrolimus clinical multiomic cohort; no new data
request as a prerequisite. Owner-authorized Human Liver Cohort
historical results are identified as *not a new open-data run*.

## Core rule

```text
sampled biological material (what/where/from whom)
  + physical fraction (plasma vs white cells vs liver tissue vs microsomes)
  + method (gene RNA / exon RNA / targeted peptide / MS metabolite / probe activity)
  + time relative to dose (predose vs postdose; original culture time)
  + genotype material source (original recipient/organ donor; allele coverage)
  -> physiology/PK process (hepatic metabolism, gut first-pass,
       renal tubular secretion, systemic erythrocyte partition, immune PD)
  -> *source qualification* and prospective omic hypothesis
  -> original heldout drug outcome + PGx reference, if actually available
```

**DNA source and expressed biology are different.** An authentic
germline variant from peripheral nucleated blood cells may reflect a
person's inherited genome (subject to transplant origin/chimerism),
but **CYP3A4 mRNA measured in blood is not CYP3A4 mRNA in the liver**.
For liver transplantation the liver donor's CYP3A5 genotype can
influence transplanted hepatic enzyme expression, whereas recipient
intestinal CYP3A can affect oral first-pass. Hematopoietic stem-cell
transplantation can make blood-cell DNA donor-derived, so apparent
blood genotyping may not represent the recipient's germline tissue
without chimerism/source qualification. No claim these effects were
measured in the current original public sets.

## What is actually sampled in existing evidence

| Original dataset | Actually sampled tissue/compartment, preparation, technique | Available biological endpoint | Current source and endpoint limit |
|---|---|---|---|
| [HLC / Sage syn4499](https://www.synapse.org/Synapse:syn4499) | **Human liver tissue** for gene microarray; same-donor **liver microsomal preparations** for CYP substrate activities; genuine limited array SNP calls in donor tissue. | 10 observed microsomal CYP substrate turnover readouts, including shared CYP3A. | Previously owner-authorized matching results n387–397 and n171–172 RNA+SNP, **not newly open source patients**; CYP3A substrate ≠ CYP3A5-specific; no patient plasma PK. |
| [Latham 2025 PHH](https://pmc.ncbi.nlm.nih.gov/articles/PMC12990092/) | Isolated **primary human hepatocytes**, original donor **cell pellets** for genotyping, **targeted CYP protein QTAP**, vendor CYP2C8 and shared CYP3A enzyme probes. Drug depletion in PHH incubation medium is outcome, **not endogenous metabolomics**. | Original 14 PHH donors; n13 complete targeted protein, independent imatinib *in-vitro* apparent clearance. | All source replayed; original hepatocyte assays ≠ liver biopsy in patient or systemic CL/F; only n9/3/1 CYP3A5 phenotype groups in protein-complete subset. |
| [Rotroff metformin 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4906013/) | **Heparin plasma**, GC-TOF untargeted metabolomics: **A** 8 a.m. day 1 after overnight fast *before metformin*; **B** 7:30 a.m. day 2, 12.5h after first dose; **C** 2h after second dose. | Original plasma metformin concentrations, n32 measured AUClast0–24 and CLr among n33 sampled persons, 384 signals. | Only A is baseline/pre-treatment. B/C cannot predict baseline PK as if predose; they could be a separate prospective on-treatment landmark problem with a later endpoint and original time validation. Original full donor PGx **not available**; plasma features are NOT directly measured renal OCT2/MATE activity. |
| [EATRIS-Plus v4 125](https://zenodo.org/records/19136013) | **White blood cells**: WGS/full-germline access restricted; mRNA-seq and miRNA. **Plasma**: proteomics, targeted metabolomics, lipidomics. Some cell composition metadata public. | Healthy-person multiomic reference and immune/plasma network integration. | Open preprocessed multiomics only, full genomic genotype/phenotype controlled; **no independently measured CYP3A drug PK**. Plasma protein is not liver intracellular CYP. |
| [GTEx V8 original medians](https://gtexportal.org/home/downloads/adult-gtex/) | Bulk **liver**, bulk **whole blood**, terminal **ileum** RNA gene-TPM *tissue aggregate medians*. | Relative gene-tissue expression context only. | Original public DCC source benchmark reproduced; not PBMC, purified hepatocytes, enzyme function, source same-donor RNA correlation or individual drug PK. |
| [Zhong 2020 GEO RNA](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124074) + [GSE147628 genotype](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147628) | Original **human primary hepatocyte RNA-seq** and genotyping-array sample metadata; 60 RNA / 96 genotype titles with 60 title overlaps. | Published CYP3A5 rs776746 exon-/splice-related association, not our own reconstructed donor-level functional model. | RNA donor sample-title matching alone does **not** qualify actual individual rs776746 calls, exon junctions, enzyme activity or drug PK. |
| [NCI-60 original 3-omics](https://pmc.ncbi.nlm.nih.gov/articles/PMC6889472/) | Cultured **cancer cell lines**, original transcriptome/metabolome from IntLIM and SWATH-MS proteins from Guo *different laboratories, batches and acquisition periods*. CellMiner somatic mutation flags only. | 47 common cell-line identities and measured in-vitro doxorubicin GI50 response. | Exact cell-line label match is not one simultaneously aliquoted specimen; not human patient PK or germline CYP3A5 response. |
| [Amorosi CYP2C9](https://doi.org/10.1016/j.ajhg.2021.07.001) | **Engineered yeast cells** for Click-seq variant activity; **engineered cultured human cells** for VAMP-seq variant protein abundance. | 4,421 original matched *protein constructs*, activity versus abundance. | Same engineered amino-acid substitution ≠ same biological specimen/patient and not measured drug plasma CL. |
| [Suiter NUDT15](https://pmc.ncbi.nlm.nih.gov/articles/PMC7071893/) | Independent engineered **HEK293T cell** protein stability VAMP-seq and post-thioguanine cellular tolerance assays. | 2,844 measured matched *construct variants*. | Cell survival PD ≠ substrate-specific hydrolase assay, measured patient toxicity or PK. |

### Direct public tissue check: GTEx V8 (measured, not predicted)

Original bulk-tissue **median** TPM, from the qualified DCC source:
- CYP3A4 liver **335.297**, whole blood **0.0501182**, terminal
  ileum **48.2465**;
- CYP3A5 liver **155.853**, whole blood **0.0870020**, terminal
  ileum **61.4648**.

See [original live GTEx replay evidence](GTEX_V8_CYP3A_BLOOD_LIVER_ORIGINAL_TISSUE_SCREEN_20261104.md).
This supports the *negative-control design*: circulating leukocyte
CYP3A mRNA must **not** be relabelled liver CYP3A expression. It does
**not** show whether an indirect blood inflammatory gene signature
could predict hepatic CYP function or drug PK in a paired patient study.

## Pharmacological interpretation: physical measurement versus model target

| Target drug and process | Preferred explanatory tissue and assay | Actually feasible accessible patient measurement | What must be marked indirect |
|---|---|---|---|
| Tacrolimus hepatic CYP3A4/5 CL/F | Liver hepatocyte RNA **splice-resolved** (for CYP3A5*3), isoform-specific active CYP protein, independently verified substrate activity; genotype origin matters | Whole-blood C0, genotype from qualified recipient/donor source, plasma/blood signatures at measured time | Whole-blood leukocyte RNA; plasma proteins or metabolite peaks; intestinal CYP expression is first-pass, not the liver. |
| Tacrolimus intestinal first-pass bioavailability | Proximal intestinal enterocyte CYP3A/ABCB1 assays, if authentically measured | Oral exposure and genotype with appropriate donor origin, diet/drug interactions | Ileum GTEx whole-tissue average is **not** proximal enterocyte protein or individualized F. |
| Tacrolimus whole-blood distribution | Hematocrit, erythrocyte binding, measured plasma and whole-blood tacrolimus with documented matrix | Blood C0 and Hct | Plasma and whole-blood C0 are not numerically interchangeable; a whole-blood drug concentration is an outcome, not a new omic. |
| Metformin renal active secretion | Renal tubule OCT2/MATE1/MATE2-K transport, qualified kidney-derived assay | Predose plasma metabolome, measured metformin AUC/CLr and renal function | A plasma signal represents mixed production, distribution, diet, microbiome and excretion. It is not renal transporter activity without validation. |
| Drug response in NCI-60 | Cell-type-specific cancer biology | Cultured cell-line RNA/protein/metabolome and GI50 | Predictive signal in culture cannot be renamed patient enzyme CL or clinical efficacy. |
| Tacrolimus immune pharmacodynamics (distinct from PK) | Activated recipient/donor immune-cell subsets, e.g. T cells; appropriately phenotyped intracellular signaling | Qualified blood immune-cell RNA/protein and clinical immune endpoints **if** truly available | White-blood-cell aggregate mRNA is a mixed-population proxy, not purified T-cell calcineurin activity, rejection diagnosis or measured drug clearance. Public EATRIS has no tacrolimus response endpoint. |

Other measurements: **urine metabolites** can reflect excretion,
renal physiology and normalized analyte recovery but do not localize
which renal transporter changed; **serum proteins** do not equate
intracellular abundance; **PBMC** are nucleated immune cells while
**whole blood** mixes cell populations; hemoglobin/hematocrit can
influence tacrolimus distribution without encoding hepatocyte
function. Isolated intestinal mucosa, tissue biopsies and liver
microsomes are biologically informative but not routine predose
patient blood samples; do not assume they exist in a public study.

## Two independent axes: mechanistic value and clinical feasibility

A direct hepatic CYP3A5 exon-junction or specific active-protein assay
may be biologically **more specific** than circulating plasma, but a
liver biopsy is not a clinically trivial routine dosing measurement.
Conversely **blood/plasma/urine** can be obtained noninvasively or
with comparatively low invasiveness but may only carry **indirect
host-wide physiology**. Accordingly `qualify_public_specimen_for_PK()`
reports **mechanistic source-to-organ proximity** separately from
`potential_clinical_collection_feasibility`. Neither criterion
automatically chooses the omic, estimates sample costs, measures
assay accuracy, or proves clinical dosing utility.

| Biological source | Physiopharmacological proximity | Clinical practicality |
|---|---|---|
| Hepatocyte CYP3A5 splicing / independently selective CYP enzyme activity | High for **hepatic CYP3A-specific** mechanism if assay original and specific | Direct original liver specimen or ex-vivo hepatocytes; **not** assumed routinely available in PK care |
| Whole-blood / PBMC gene expression | Directly reflects assayed blood-cell biology, **not** hepatocyte CYP3A function | Potentially collectable peripheral blood; cell mix and assay analytical validity must be established |
| Predose plasma metabolite | Systemic functional proxy, can reflect multiple tissues and exposure contexts | Clinically accessible in principle; metabolite structure, source timing and metabolic specificity still needed |
| Urine metabolite | Excretion/renal/host-wide proxy, not OCT2/MATE2 activity itself | Noninvasive collection in principle; collection interval and urinary dilution/creatinine interpretation essential |
| GTEx postmortem/reference-tissue median | A tissue atlas context, **not a patient measured biomarker** | Not an individual predose clinical assay |
| NCI-60 / engineered DMS constructs | Experimental cell-culture mechanism | Not the original patient's biospecimen |

## Implementation: fail-closed public-only source gate

`src/omics2pk/specimen_physiopharmaco_gate.py` contains
`SOURCES` fixed source-specific original tissue/fraction/method/time
annotations and `qualify_public_specimen_for_PK(...)`. It returns
whether a measured source:
- is anatomically in the target organ, an **ex-vivo human donor preparation**,
  a circulating systemic proxy, a tissue mismatch, an engineered cell line
  or a published aggregate without patient samples;
- was collected before treatment or **after** dosing;
- actually carries an independently measured **human plasma PK** outcome;
- has authentic **full pharmacogenetic coverage** (none in the listed
  fully open original cohorts).

`qualify_public_omic_model_arms(...)` additionally refuses to
merge liver RNA from an HLC donor and PHH protein measured in a
different Latham donor, to treat GTEx medians as patient covariates,
or to pass Rotroff postdose B/C as pretreatment predictors.
The original Latham n13 in-vitro heldout arms now attach the
actual **primary hepatocyte** sampling evidence next to the
MAE/log-MAE/P20 comparison and explicit insufficient genotype
strata. Source-level model-arm qualification is **not itself**
a verified sample-level join or a validation of clinical benefit.

**Read-only cross-project handoff now implemented:**
`qualify_multiomics_specimen_handoff(...)` accepts the optional
MultiOmics `omics_specimen_context_v1` (or the integrated
`reconcile_variant_with_original_specimens` report). It independently
checks allowed original organ/fraction combinations, source assay
methods, pre/postdose timing, and whether material is a
human specimen, cancer cell line or genetically engineered construct.
Output labels direct local tissue **context**, indirect systemic
proxy, or ineligible experimental/post-treatment measurements; it
never creates numerical human omic values, clinical genotype,
measured tacrolimus PK or patient-level model validation.

The matching MultiOmics component in
[MultiOmics PR #121](https://github.com/rberrah/MultiOmics/pull/121)
requires original assay fraction (e.g. white cells vs plasma,
liver tissue vs microsomes), method, dose-relative time,
biology provenance (actual patient vs engineered construct)
and sample pairing proof. It must *not* change frozen
`DOWNSTREAM_EXCHANGE_V1`.

## EATRIS-Plus: public original cross-fraction analysis now quantitative

[MultiOmics primary original-source evidence](https://github.com/rberrah/MultiOmics/blob/feat/variant-omic-functional-evidence-sidecar-20261105/docs/multiomics/EATRIS_PLUS_V4_PUBLIC_WBC_PLASMA_EVIDENCE_20261105.md)
now includes a source-verified original MuData benchmark. The
published Zenodo v4 141.7 MB `EATRIS_PLUS_mudata.h5mu`
(MD5 `bbf0342797adcbcd6a3c5b242ca64573`)
contains **125 overall subjects**, but **90** exact original
persons match all four study assays used for prediction:

- **white blood cell** `mRNA`: 16,386 original quantified
  gene features;
- **plasma** acylcarnitines (33) and amino acids (46):
  79 targeted circulating metabolite features, **not**
  hepatocyte enzyme flux and not direct hepatic metabolite
  levels;
- independently assayed **plasma** `proteins`:
  1,393 original source feature columns, of which **80**
  measured protein targets are preselected from training
  samples alone in each outer fold;
- all physical specimens belong to the same source person
  by original index, **not proven same aliquot or
  simultaneous collection**.

The phenotype is the original **standardized plasma protein
profile** (cross-omic reconstruction), **not** C0, AUC,
individual CYP3A function, tacrolimus clinical PK or PGx.
The earlier n90 repeated subject-heldout benchmark showed
a mean-only baseline 0.760615 versus a poorly regularized
RNA+metabolite model at 0.819299 (different original
PC-per-layer design).

**New independently executed original-source 2×2 model
sensitivity** [CI SUCCESS #42068596760](https://github.com/rberrah/MultiOmics/actions/runs/42068596760)
holds the *same total* 3 versus 9 training-only
omics PCs for each predictive arm, and varies fixed
Ridge alpha 10 versus 1,000,000:

| WBC and/or plasma molecular predictor | 3 PCs, alpha=10 | 9 PCs, alpha=10 | 3 PCs, alpha=1e6 | 9 PCs, alpha=1e6 |
|---|---:|---:|---:|---:|
| None, source train mean | **0.760615** | **0.760615** | **0.760615** | **0.760615** |
| WBC mRNA alone | 0.783413 | 0.807777 | 0.760637 | 0.760636 |
| Plasma targeted metabolites alone | 0.773407 | 0.794645 | 0.760611 | 0.760611 |
| WBC + plasma | 0.777599 | 0.808840 | 0.760627 | 0.760633 |

**Within-model factor isolation:** at alpha 10, adding
six PCs to the combined WBC/plasma predictor worsened
the same-person standardized protein MAE by
**+0.031241** (original-person paired descriptive
95% bootstrap [+0.022673,+0.039824]). At 9 PCs,
increasing fixed alpha from 10 to 1e6 reduced
joint error by **−0.048207**
[−0.060455,−0.036572], but reduced the molecular
prediction almost to the constant **no-omics** mean.

**Interpretation:** poor out-of-sample performance in
n90 is consistent with finite-sample feature/model
estimation and weak or redundant **cross-compartment
signal**, not intrinsic information toxicity. Excessive
shrinkage suppresses usable molecular variation and
must **not** be interpreted as successful target-organ
functional inference. The choice of alpha 1e6 was
motivated by the earlier n90 model behavior, making
this test **exploratory/post-hoc, not an independent
confirmatory cohort**. Original-person bootstraps do
not adjust model-choice multiplicity. The comparison
cannot demonstrate that WBC RNA, plasma proteins or
plasma metabolites help human tacrolimus PK **beyond
true clinical+complete PGx**, because the public
source has **neither PK exposure nor full public
germline genotype**.

A critical Omics2PK implication is that **biological
proximity and assay-method practicality** should be
recorded as candidate priors for feature selection,
not a substitute for heldout drug-specific
incremental predictive gain.

## No additional human cohorts assumed: what to do using public data

1. **E1: tissue/site mechanism and assay consistency**, already
   testable on original GTEx V8 medians, public original
   CYP2C9/NUDT15 variant/function screens, and GEO CYP3A5 hepatocyte
   source metadata. Publish both positive and negative examples
   *without claiming individual clinical PK*.
2. **E2: within-same-PHH origin sensitivity and endpoint metric
   comparison**, original Latham n13; compare source PGx vs
   proteins vs measured shared probe activities using log-error,
   RMSE, P20 and paired interval with explicit n9/3/1 genotype
   class abstention. No full RNA or metabolomic matrix.
3. **E3: genuine human systemic PK with blood-derived metabolome**,
   original Rotroff **predose plasma A, n32** complete PK. The
   original source has no qualified full PGx comparator: report
   **clinical-only vs predose plasma** as an *initial E3*, not the
   project's final genotype-over-omic claim. Use B/C strictly
   for separate longitudinal effects.
4. **E1 cancer PD:** original NCI-60 47 source-matched cell line
   RNA+protein+metabolome+somatic mutation; compare models and
   perturbation nulls, but do not import those cell-line features
   into real patient PK.
5. **Public EATRIS-Plus:** evaluate multiomic integration
   conditional on **white-cell vs plasma origin**, original
   missingness and cell composition; this is not a tacrolimus
   PGx/PK cohort, and restricted full genome must not be
   imputed or treated as public.

**No new BIOMARGIN/ABIS/tacrolimus patient materials assumed.
No follow-up human cohort acquisition is a prerequisite to finishing
the research-only methodology and source-replay qualification.**
