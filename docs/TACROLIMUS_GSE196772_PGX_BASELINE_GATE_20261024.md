# Tacrolimus clinical + PGx C0/D baseline source gate — GSE196772 / Yang 2022

**Updated:** 24 October 2026  
**Study:** Yang et al., *Clinical and Translational Science* 2022  
**DOI:** `10.1111/cts.13388`  
**GEO:** `GSE196772`

## Why this matters for Omics2PK

The primary Omics2PK tacrolimus question is not whether omics alone correlates
with tacrolimus exposure. It is whether a molecular phenotype adds predictive
information **beyond an already strong clinical + pharmacogenetic reference**.

Yang et al. provide a useful published reference in **251 Chinese kidney
transplant recipients**.

The article reports stable tacrolimus trough exposure using:

`log(C0 / (daily dose / body weight))`

and records age, sex, body weight, hematocrit and other clinical variables.
The public GEO accession contains the genotyping source.

Published multivariable explanatory R2:

| Published arm | R2 |
|---|---:|
| Age | 11.03% |
| Age + sex | 15.76% |
| Age + sex + HCT | 17.89% |
| + CYP3A5*3 / rs776746 | 35.12% |
| + rs75125371 | 37.34% |

The approximately **17.2 percentage-point increment from CYP3A5*3** is exactly
why an Omics2PK study must not compare an omic score only against a clinical
model that omits PGx.

## Critical evidence boundary

These R2 values are **published in-sample explanatory regression results**.

They are not:
- held-out prediction performance;
- external validation;
- dense concentration-time PopPK;
- AUC or CL/F;
- evidence that an omic phenotype improves over genotype.

GEO `GSE196772` publicly exposes SNP-array data for 251 recipients. The paper
states that the individual tacrolimus dose, C0 and clinical information were
retrieved from medical records, but the GEO release and the inspected public
supplements do not establish a same-subject patient-level clinical PK table.

The supporting information is described as:
- Figures S1-S4;
- Table S1;
- Table S2: SNPs/association results below the paper threshold;
- Table S3: candidate-SNP association analyses.

None is a verified 251-row patient × C0 × dose × weight × HCT table.

Therefore the current source status is:

`PGX_GENOTYPES_PUBLIC__INDIVIDUAL_CLINICAL_C0D_VECTOR_NOT_ACQUIRED`

No individual C0/D values may be reconstructed from means, SDs, figures,
regression coefficients or genotype-group plots.

## Minimal future source table

The source gate
`src/omics2pk/tacrolimus_pgx_c0d_source_gate.py`
accepts one row per biological recipient containing:

```text
subject_id
age_years
sex
weight_kg
hematocrit_fraction
c0_ng_ml
daily_dose_mg
stable_maintenance
rs776746_effect_allele
rs776746_effect_allele_dosage
[rs75125371_effect_allele]
[rs75125371_effect_allele_dosage]
```

Important:
- hematocrit must be converted explicitly to a fraction; percent/per-mille is
  never guessed;
- genotype effect-allele orientation must be constant across the cohort;
- stable-maintenance status must be source documented;
- biological subject identity must exactly bridge the genotype and clinical
  record;
- no row-order merge is allowed.

The gate derives:

`ln(C0_ng_ml / (daily_dose_mg / weight_kg))`

only from observed same-subject source values.

## Role in the Omics2PK hierarchy

This study is a **PGx baseline source**, not the final Omics2PK cohort.

When a genuinely paired tacrolimus omics cohort becomes available, the minimum
comparative structure should be:

1. clinical;
2. clinical + CYP3A5 PGx;
3. clinical + full available CYP3A PGx;
4. clinical + PGx + direct measured molecular phenotype;
5. clinical + PGx + SWEET / mechanistic molecular score, if RNA exists.

All arms must use identical patients and held-out splits.

If the molecular arm does not beat the PGx arm, that is a valid negative result.
The objective is incremental information, not forcing an omics benefit.

## Population-specific caution

rs75125371 was reported as a population-specific variant in this Chinese cohort.
It should not automatically be required in a European tacrolimus reference
model. CYP3A5*3 is the more portable minimum comparator, while ancestry-relevant
CYP3A5 *6/*7, CYP3A4*22 or other variants should be incorporated only when
genuinely measured and justified for the target population.

## Current conclusion

GSE196772 materially strengthens the **definition of the comparator** for
Omics2PK, but does not remove the paired omics+PK data bottleneck.

The claim currently supported is:

**clinical + CYP3A5 PGx is a strong tacrolimus C0/D reference that any molecular
phenotype should be tested against.**

The claim not yet supported is:

**molecular phenotype outperforms genotype for tacrolimus PK.**
