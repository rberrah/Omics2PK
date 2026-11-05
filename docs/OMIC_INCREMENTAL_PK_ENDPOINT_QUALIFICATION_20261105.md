# Genotype baseline + omic layers: outcome-aware, fold-safe incremental evidence

**2026-11-05 · Omics2PK method · research-only software qualification.**

The scientific question is **not** whether a mutation alters a transcript
or a protein; that belongs to **MultiOmics and orthogonal functional evidence**.
The Omics2PK question is whether measuring a molecular layer improves
**an independent drug- or enzyme-specific heldout endpoint above the
original patient's genuine measured clinical+pharmacogenetic baseline**.
Hypotheses from the MultiOmics variant sidecar are *testable predictors*
only when the underlying quantitative original sample-level
molecular values, timing and provenance are independently validated.

## 1. Input contract

Implemented in
`src/omics2pk/omic_incremental_endpoint_audit.py` as
`audit_omic_incremental_PK_prediction(...)`.

A qualifying input must provide **one original biological subject per
record**, the original heldout fold, one independently observed strictly
positive PK/function endpoint and exactly the same corresponding
**heldout** predictions for clinical+actual PGx and every measured
additional omic arm. This source-attestation API does not independently
reconstruct/fingerprint the model training splits. It refuses:

- **No measured genotype** mislabelled as a genotype baseline; one source
  SNP is marked explicitly as **incomplete PGx**, not silently upgraded
  to a star-allele diplotype.
- An omic measured **after** the intended prediction landmark,
  predicted on another subject or fitted on that heldout subject;
  these are *caller-asserted original qualification flags* whose
  scientific truth still requires original sample and split audit.
- Different model arms evaluated on differing subjects, negative/zero
  targets in a positive-clearance log metric, duplicated biological
  subjects, missing fold identities or incompatible units/endpoint
  names.
- Arbitrary inference of a **mutation's function** from raw RNA/protein;
  known-function variant strata must be externally established
  **before** analyzing the response.

## 2. Predetermined pharmacological outcome and evaluation

### Clinical PK concentrations/exposure

For positive endpoints like true plasma AUC, CL, C0 or Cmax:
- **Primary candidate: mean absolute log(predicted/observed)** when
  relative exposure error is clinically meaningful, robust to the
  choice of concentration unit. A 2-fold error gives the same absolute
  log loss whether it is 2× or ½×.
- Linear **MAE** in the original outcome units and **RMSE**
  as secondary outcomes; these weight large absolute errors and
  extreme high-concentration/clearance subjects differently from
  the log metric.
- **P20** = fraction with predicted/observed ratio in [0.8,1.2],
  and **within factor two** as intuitive accuracy measures. **These
  are not pharmacological target-attainment probabilities**.
- Signed log-bias and log RMSE to distinguish systematic direction
  and outliers.
- If the scientific objective is dosing harm, target attainment,
  exposure-toxicity or clinical utility, specify these outcomes and
  eligibility explicitly and add their validated decision metrics:
  none are inferred from generic positive exposure alone.

**MAE is never automatically “the” best metric.** Predeclare the
primary criterion for each endpoint based on drug dosing context and
what an error means clinically, *before* outcome-based omic selection.

### Biological function and cell PD

The same numeric error module can qualify **original positive
enzyme activity** or **PHH drug depletion clearance**, but clearly
labels them `BIOLOGICAL_FUNCTION_NOT_CLINICAL_PK` and
`IN_VITRO_DRUG_METABOLISM_NOT_CLINICAL_PK`. It **cannot** evaluate
binary toxicity or survival as if they were positive PK; such clinical
PD endpoints need a separate classification/survival module
(AUROC/AUPRC, Brier/calibration, net clinical benefit and time
bias controls). A cell-line thiopurine tolerance score is neither
patient toxicity nor drug AUC.

## 3. Original heldout, paired scientific comparisons

All prespecified models get the same biological subjects and
same endpoint, with training-only imputation, feature selection
and any genotype×omic interactions; the outcome must never be
used to choose the omic ahead of the outer validation fold.
The generic evaluator returns:

- MAE, RMSE, mean absolute log ratio, log-RMSE,
  signed log bias, P20 and factor-two accuracy *for every arm*;
- delta of **the predeclared primary loss** for each addition
  versus clinical+PGx (negative is better);
- exploratory subject-level *paired* bootstrap
  (499 resamples), **not** model refits, multiplicity-adjusted
  inference or proof of generalization;
- descriptive variant-function-stratum comparisons only if the
  source has **at least 20 independent subjects in that stratum**;
  otherwise explicit **ABSTAIN**;
- **no automatic best-omic / clinical-benefit declaration**, even
  when one numeric point estimate is smaller.

The module distinguishes a true independent external cohort from
cross-validation on one original source. It cannot itself certify
which participant rows were used to train externally supplied
predictions; source and model logs remain necessary.

## 4. Real original example: Latham primary hepatocytes

The existing original Table 4
[official source replay #42015617333](https://github.com/rberrah/Omics2PK/actions/runs/42015617333)
has 13 same-donor CYP2C8+CYP3A5 source PGx + targeted proteins
+ independently measured shared CYP3A/CYP2C8 probe activities
+ *in-vitro* apparent imatinib clearance. Models used
the same leave-one-donor-out splits, predictions were from
train-only Ridge on log clearance; outcomes and source individual
rows were not committed.

| Original predeclared arms | MAE original units | Log MAE | P20 |
|---|---:|---:|---:|
| Source PGx-only | 2.906576 | 0.432147 | 5/13 |
| + proteins | 3.003704 | **0.408474** | 5/13 |
| + probe activity | **2.762944** | 0.415360 | 5/13 |
| + both | 3.262812 | 0.445258 | 5/13 |

All previously qualified original 95% donor-bootstrap differences
versus source PGx include zero. **Linear MAE chooses activity;
log-MAE chooses protein; P20 ties. The output does not support
a superior modality**. The reusable evaluator was **successfully rerun against official
original Latham Table 4** in
[GitHub Actions #42017839096 (retry SUCCESS)](https://github.com/rberrah/Omics2PK/actions/runs/42017839096),
after one transient external PMC source-unavailable abstention.
It reproduces the same 13 source donor-level heldouts and
refuses all CYP3A5 genotype-stratum winners:
**n=9** known no-function, **n=3** reduced-function and
**n=1** function-allele carrier among the 13 protein-complete
donors; each is below the minimum n20 per stratum.
The primary log-error paired 499-donor bootstrap for
**probe activity** minus PGx baseline is **−0.016787**
[−0.131239,+0.106452] and for **proteomics** is
**−0.023673** [−0.144389,+0.145689], both including zero.
These are descriptive heldout error contrasts and do not
justify an assay-choice decision or clinical dose effect.
The original PGx panel is only partial and the CYP3A
midazolam probe is CYP3A4/5-shared, not CYP3A5-specific.
Neither the outcome nor this 13-donor comparison is clinical PK.

## 5. Larger human E2: original HLC already analyzed, source locally controlled

An earlier owner-authorized
[Human Liver Cohort source-matched RNA+SNP+activity analysis](HLC_ORIGINAL_400_LIVER_CYP_RNA_GENOTYPE_RESULTS_20261018.md)
found n171–172 donors for genuine single-SNP+hepatocyte
RNA+microsomal CYP function comparisons. Incremental RNA
reduces heldout MAE over age/sex+SNP for the three measured
substrate-activity endpoints, **but** a single SNP is not
complete CYP3A5/CYP2C9 PGx, original function assays are
not clinical drug PK, and the source authorized data cannot
be committed as original benchmark assets. This previous
analysis was not newly recomputed from original bytes in
this PR.

## 6. Two-way scientific handoff to MultiOmics

```text
Genome allele (+ transcript, strand, splice/codon consequence)
  ↓
MultiOmics: RNA/protein/metabolite/functional readouts
  ├─ observed original feature/assay signal, tissue and time
  ├─ discordance and biochemical *hypothesis*
  └─ independently graded literature evidence / mechanism uncertainties
  ↓ [Research-only variant evidence sidecar AND source original matrices]
Omics2PK:
  clinical + genuine PGx/PopPK baseline (zero omics)
  + variant-conditioned measured molecular representation
  + no/one/all modality black-box alternatives
  ↓ independent patient-heldout genuine PK/PD
  paired endpoint-appropriate errors, uncertainty and generalization
```

A plausible mechanism can be **predictively useless** after the
genotype is already modeled. Conversely a metabolite of
unknown identity can predict PK but cannot acquire a fictitious
enzyme identity. A predeclared omission/no-omics arm is
required for every comparison.

**Status:** the generic evaluator and synthetic anti-leakage
tests are implemented, with a source-real PHH adapter.
No full same-person genotype+RNA+proteome+metabolome+real
patient tacrolimus AUC/CL/C0 external study has yet been
validated, and the evaluation module does not change that.
