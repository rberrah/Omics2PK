# Tacrolimus-first CYP3A4/CYP3A5 Omics2PK research qualification (2026-10-19)

## Frozen scientific priority

**PRIMARY DRUG: tacrolimus; enzyme system CYP3A4 + CYP3A5 (with intestinal/hepatic first pass, blood partitioning and efflux transport context).** Midazolam remains an **auxiliary CYP3A functional probe** with measured clinical oral PK, and is NOT itself a transplanted patient's omic phenotype or tacrolimus PK result. Other CYP3A substrates (e.g., ciclosporin, everolimus) are lower priority pending exact patient-matched source cohorts. No pharmacometabolomic benefit is assumed in advance.

The aim is to test `clinical + real CYP3A5/CYP3A4 PGx` vs `clinical + PGx + prospectively measured omic phenotype` on **identical true tacrolimus-treated people with observed PK**. Full dosing/formulation, hematocrit, time-after-transplant, drug interactions, age, weight, kidney/hepatic function and patient-level nonduplicated PK observations must be accounted for. For oral tacrolimus usually `CL/F` is identified, not `CL` and `F` separately. CYP3A5 *3 alone is an incomplete multiancestry diplotype panel; *6/*7 and CYP3A4 *22 should be explicitly measured or flagged missing when relevant. A CYP3A5 genetic allele is not a measured enzyme activity.

## Native MultiOmics V1 ZIP already generated locally (private)

Owner-provided original Jarmusch GNPS QIIME2 study n13 pre-dose plasma MS1 samples have a same-person original CYP3A4/CYP3A5 PGx and **MIDAZOLAM**, not tacrolimus, clinical PK outcome. A separate actual **native source-equivalent no-signature exporter** was executed in Node using the pinned `rberrah/MultiOmics` `src/lib/multiomics/downstream-exchange.js` no-score branch at commit `cd101839e8a9fae29c88b3bfd86ea350088af282`, source Git blob `8aad6988c6c6ec9788bd96c75aa6ea0a2887f52e`.

The private user-facing archive `CYP3A_MultiOmics_to_Omics2PK_research_bundle_20261019.zip` contains genuine `exchange/samples.csv`, `features.csv`, `matrices/metabolomics.csv`, `manifest.json`, plus a separate private clinical-midazolam sidecar and a tacrolimus-specific **study design spec**, NOT invented tacrolimus outcomes. Local **numerical source-to-native-export equality: 82,537 / 82,537 areas identical**, 6,349 MS1 features, 13 unique true biological people; 6,089 features positive in at least one baseline plasma assay. SHA256 and referential schema checks passed. Original per-person data, ZIP and biosample keys are deliberately NOT part of public GitHub, Actions logs or artifacts.

**Producer provenance limitation:** Vite-specific `import.meta.env` is absent outside a browser build; the local Node no-signature execution recorded `software.commit='unknown'` instead of forging a Vite commit. The actual upstream Git source commit and blob IDs are separately source-qualified in the local bundle's technical QC. This is a **real source-compatible V1 ZIP built outside the browser**, not evidence that the deployed MultiOmics browser exported exactly the same original bytes; it is neither a true user-interface E2E test nor a PopPK model.

At baseline only (`T0` relative sample label), with no verified original absolute collection clock; no dates are fabricated. **GNPS MS1 peak IDs are unidentified spectral features and cannot be fed to RNA-only SWEET or promoted to confirmed CYP3A enzyme abundance.**

## Existing real-PK CYP3A pilot and failed hypotheses

The exact original midazolam patient-level PK LOPO comparator (`scripts/validate_local_jarmusch14_gnps_biom_pk.py`) and its independent input package replay (`scripts/validate_private_jarmusch13_multiomics_pipeline.py`) show measured midazolam J1 (n13) and J9 (n12) numerical source scale:

| Study endpoint | Clinical sex+weight | Clinical + authentic CYP3A4/5 PGx | + two fold-local unsupervised MS1 PCs | Interpretation |
|---|---:|---:|---:|---|
| Midazolam oral CL/F day1 (n13) | 6.8753 | 7.3434 | 9.8880 | MS1 PCA made prediction worse |
| Midazolam oral CL/F day9 (n12) | 3.9843 | 4.6012 | 4.7488 | MS1 PCA made prediction worse |

The additional **exploratory 9-peak bile-acid-class MS2 spectral match proxy** was also worse (J1 +2.6318, J9 +0.5182 MAE vs PGx). That panel represents putative co-annotated MS1 signals, not nine proven distinct biological molecules, not measured CYP3A flux, and was explored after this dataset was examined. It cannot furnish a confirmatory result. All patient-level data stay private.

**Conclusion:** we have a working pipeline and a *negative* true-human CYP3A-probe PK benchmark, not proof that biological phenotype improves tacrolimus PK.

## Tacrolimus data acquisition and modelling priorities

| Priority | Evidence source | True existing measurement | Missing to prove the primary paired claim |
|---|---|---|---|
| P0 — Tacrolimus direct PK+pre-dose metabolome | Phapale 2010, DOI `10.1038/clpt.2009.296`, n29 healthy volunteers | Pre-dose urine LC-MS 1,256 metabolite ions; tacrolimus plasma C(t) up to 72h, real observed individual PK in publication | **Original numerical 29 × 1,256 matrix; dose/time/C(t) and/or individually observed AUC, shared original subject IDs; authentic CYP3A4/CYP3A5 genotype only if actually available**. All 13 historical Wiley DOC supplements obtained previously did **not** include those complete individual arrays. |
| P0 — Tacrolimus PK and full CYP genotype comparator | Brazeau 2025 `10.1111/cts.70370`, n65 transplant; Tsironi 2026 `10.1038/s41397-025-00395-w`, n94 with 5 visits | True dose-adjusted AUC or serial normalized troughs, source CYP3A genotype information; **no full patient-level omic abundance crosswalk independently secured** | Same-person raw original PK+PGx table and separately measured pre-dose molecular state; avoid treating article group summary as n65/n94 patient matrix. |
| P0 — Public tacrolimus clinical+PGx C0/D reference | Yang 2022 `10.1111/cts.13388`, **GSE196772**, n251 kidney transplant | Public SNP-array genotypes; published stable C0/dose/weight phenotype and clinical covariates. Published in-sample R2: clinical age+sex+HCT **17.89%**, +CYP3A5*3 **35.12%**, +rs75125371 **37.34%** | GEO does **not** establish the individual 251-row clinical C0/dose/weight/HCT vector. Use as comparator definition/source gate only until same-subject clinical rows are legitimately obtained; do not call published R2 held-out prediction or PopPK. |
| P1 — Tacrolimus CYP3A5 genotype measured drug/metabolite PK | Thummel et al., healthy n24, `10.1038/clpt.2012.195` or source original PMID `23073208` | Original observed blood tacrolimus PK and 31-DMT/12-HT/13-DMT metabolite:parent AUC ratios by CYP3A5 genotype | Pre-dose endogenous molecular abundance layer and authentic subject-level PK arrays; **postdose tacrolimus metabolites must not be called baseline omics predictors**. |
| P1 — Liver enzyme quantity/activity to mechanistic score | Public Wegler, HLC/EPA liver datasets | Paired CYP proteins, enzymes/activity / CLint in human liver donors (each source independently qualified) | Original same-person tacrolimus in-vivo AUC/CL/F **not present**, cannot merge unrelated people or call E2 evidence E3. |
| P2 — Prospective transplant longitudinal clinical study | He 2022 n109 tacrolimus and longitudinal metabolomics D7/M1/M3, DOI `10.3390/molecules27113517` | Repeated molecular and clinical outcomes reported | Individual timepoint matched actual tac dose/C0/AUC, specific treatment timing, access and governance; repeated measures must split by whole patient. |

### Predeclared representative model arms (one same-person cohort)

- `C0`: a fixed reference tacrolimus clinical+PGx PopPK (e.g., Woillard-style structure) fitted to actual **dose–concentration–time** data, with dosage, formulation and relevant clinical factors.
- `M1`: clinical+PGx + expression/protein/metabolite **direct train-only representation**.
- `M3-SCORE`: clinical+PGx + **organ/drug-specific, frozen out-of-sample individual molecular CYP3A4/5 functional score**; provenance and compartment must be explicit.
- `M3-GRAPH`: same PopPK with **RNA SWEET** feature calibrated only from training donors when measured transcriptomics exists, with/without Pillar A biological prior. Negative prior/graph controls required.
- `M4`: more complex graph/PK hybrid **only after M3 evidence is robust**, to avoid p≫n overfitting, structural identifiability issues and genotype double-counting.

For an appropriately measured molecular state score `S_i`, a hypothesis is
`log(CL/F_i)=log(CL/F_pop)+beta_PGx*G_i+beta_omic*S_i+eta_i`.
This formula is a **research hypothesis**, not a validated equation or clinical dosing engine; for an oral drug, exposure can also depend on `F`, transporter actions and whole-blood binding.

**Required primary decision:** held-out per-person performance and uncertainty **increment of omics over the identical clinical+full PGx arm**, using strictly training-fold feature normalization/feature choice/score calibration, not comparison of unrelated studies or improvement of genotype prediction on CYP3A activity only. Fully independent external validation before any transportable superiority claim.

No study-source individual values, patient material or copyrighted published supplementary binaries should be pushed to this public GitHub repository.
