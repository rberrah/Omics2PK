# Omics2PK — public-source evidence usable to justify cohort access and a clinical biomarker pilot
**Research brief · 2026-11-05 · validation stage, not clinical implementation**

## The project hypothesis

Among the same patients, with genuine time-stamped clinical PK,
pharmacogenetics, and baseline specimens, **does selecting functionally
relevant molecular information improve independently held-out tacrolimus
AUC or apparent clearance prediction** beyond a credible
clinical+PGx PopPK baseline? Is that gain stable across assays and sites,
or does indiscriminate omic expansion degrade finite-sample estimation?

The incremental target is *predictive* performance, not
the number of omic predictors nor global TreeSHAP attribution.

## Three independently different public evidence bridges

| Stage | Original, verified human source | Directly observed endpoint | Main bounded result | Does it show better tacrolimus PK? |
|---|---|---|---|---|
| E1 assay sensitivity | EATRIS-Plus public v4 n88, MD5 confirmed, same-donor WBC mRNA/plasma metabolites and plasma proteome | Measured independent plasma protein assay profiles, two non-identical matrices | Under-regularized extra omics worsen interomic heldout reconstruction; strong shrinkage approaches mean-only. Relative error differs by protein target assay. | **No** |
| E2 enzyme function | US EPA n93–95 / each of 8 CYP activities, source SHA confirmed | Original human hepatic probe enzyme activity, not PK | CYP3A4 target qPCR adds ~15.1% reduction in MAE vs age/sex for the 8-RNA-complete donor subset; 7 extra CYP RNAs yield only ~0.65% additional relative reduction at alpha1, while harming 5/8 enzymes | **No** |
| E2 drug biotransformation | Wegler n15 donor-matched hepatocytes and microsomes, source SHA confirmed | Original measured in-vitro drug-specific CLint, five probe drugs | In intact hepatocytes: target CYP protein Ridge1 outperforms five-CYP Ridge1 in all 5 probe-drug endpoints; limited n, model- and assay-system-dependent | **No** |
| E3 real human PK with genotype | Jarmusch original n13 day1 / n12 day9, real predose plasma untargeted MS and original CYP3A4/5 genotypes | Independently measured in-vivo midazolam CL/F | Untargeted MS1 2-PC fusion **worsens** actual heldout CL/F MAE vs clinical+limited PGx in this exploratory small study | **No** |
| E3 real human PK without PGx | Rotroff original n32 matched predose plasma GC-TOF, verified clinical measures | Measured in-vivo metformin AUClast and CLr | No confirmed robust gain from broad/nested omic features over clinical; 384 GC-TOF features with 229 unidentified | **No** |

**Not comparable numerically:** each row uses a DIFFERENT drug,
biological specimen, analytical modality, endpoint and patient
population. Do not meta-analyze MAEs as a single efficacy metric.

Original source execution and detail:
- MultiOmics original EATRIS [#42078275344](https://github.com/rberrah/MultiOmics/actions/runs/42078275344);
- EPA [#42101612511](https://github.com/rberrah/Omics2PK/actions/runs/42101612511);
- Wegler [#42105373917](https://github.com/rberrah/Omics2PK/actions/runs/42105373917);
- [Jarmusch original GNPS report](JARMUSCH14_ORIGINAL_GNPS_BIOM_TO_REAL_MIDAZOLAM_PK_20261019.md);
- [Rotroff n32 original source and information-budget report](SHAP_INFORMATION_GAIN_PARETO_NETPK_CASCADE_20261104.md).

## Scientific inference that is justified today

- **Positive mechanistic feasibility:** matching CYP expression or protein
  abundance contains genuine predictive information for some directly measured
  liver enzyme functions and some *in-vitro* drug-metabolism rates.
- **Observed negative controls:** increasing the number of measured
  molecular covariates can increase heldout prediction error and the result
  depends on regularization, source tissue, target and substrate.
- **Clinical incremental utility:** not established. There is no qualifying
  public tacrolimus cohort with independent human PK + full validated PGx
  + measured predose molecular profiles on the same people. Do not
  retrofit E1/E2 in-vitro into tacrolimus clinical efficacy.

## Concrete pilot to request (avoid full multiomics initially)

**Patient-level minimum evidence package:** one stable deidentified
patient/patient-visit key, all tacrolimus doses/routes/formulations with
administration time, real concentration/timestamp profiles and
an independently evaluated AUC/CL/F target, true CYP3A5*3/*6/*7
and relevant CYP3A4 genotyping quality, liver and renal
biochemistry, HCT, CRP, azole exposure, transplant timing,
and a **stored pre-exposure/planned-prediction blood/plasma
aliquot** (with tube type, handling/storage history and
sample timestamp). CYP3A5 *3-only is not equivalent to full
PGx. Within-run hepatic CYP3A4/5 gene expression is not a
consequence of the patient's PBMC values.

**Measurement priority:** start with source-qualified plasma
endogenous targeted metabolomics/validated CYP3A-related
metabolite ratios sampled **before the PK prediction window**;
consider untargeted spectral features as a secondary pathway;
quantify protein abundance/proteoforms only if assay method,
specimen, independent biological functionality and cost justify it.

**Primary statistical comparison:** test, on exactly the same
heldout patients: (1) clinical existing PopPK, (2) clinical+PGx,
(3) +predeclared physiological metabolite panel,
(4) +nested selected broad metabolome, (5) fixed black-box comparator.
Within each external fold, all feature selection, normalization,
imputation, dimension choice and hyperparameter tuning remain
training-only. Do not estimate the test AUC with the same augmented
model then score the model against its own estimate; use measured
reference exposures or independent external PK evaluation.
Patient-grouped CV and site/time external replication if possible.
Compare paired MAE, bias, calibration, interval coverage,
fraction of patients meeting meaningful prediction
accuracy, and laboratory cost per incremental information gain.
A zero-feature solution is allowed.

**Funding argument:** existing public E1/E2/E3 research supports
feasibility and explains why it is scientifically necessary to
measure the **correct** molecular phenotype at the **correct
time/site** and test actual incremental clinical value. Funding is
needed primarily for qualified same-patient omics+PGx+PK specimens,
analytical measurements, independent reference PK, external
validation—not for a platform that already claims a proven
improvement in tacrolimus dosing.

## Specific source access blockers (not negative biological results)

- Wegler rosuvastatin n54: published original 54-subject
  hepatic transporter protein + individual in-vivo rosuvastatin
  PK study; source supplementary spreadsheet provides proteins
  and covariates but **not** original 54 individual AUC/C(t)
  columns. The paper's published 78% within-2-fold result
  is the *authors'* result, not independent Omics2PK validation.
- Silveira rosuvastatin n40: original raw predose
  metabolomic spectra exist at MTBLS3128; source processed
  five CSVs and individual PK crosswalk are not publicly
  retrieved. Original parent-linked PK table required.
- Shin CYP3A5 n24: genuine IV midazolam genotype,
  urinary endogenous markers and PK in the same study
  (DOI 10.1038/clpt.2013.128); per-person three-period
  numeric phenotype/source crosswalk not acquired, no
  independent comparative MAE.
- No known valid BIOMARGIN multiomic-to-ABIS full
  crosswalk. Require authorized source coordinator confirmation
  before making any linkage claim.

**Decision:** This evidence merits a limited authorized pilot,
but does **not** merit a claim that Omics2PK already improves
patient tacrolimus clearance/AUC, nor justify routine dosing.
