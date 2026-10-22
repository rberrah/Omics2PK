# Wegler 2021 — original donor CLint/protein crosswalk and held-donor benchmark

Updated: **22 October 2026**.

## Source provenance

Wegler C et al., *Molecular Pharmaceutics* 2021, DOI
[10.1021/acs.molpharmaceut.1c00053](https://doi.org/10.1021/acs.molpharmaceut.1c00053),
PMCID PMC8041379.

The study contains donor-matched human hepatocytes (HH) and human liver
microsomes (HLM) from **15 human liver donors**, five CYP probe substrates and
quantitative CYP/proteomic measurements. This is an **in-vitro drug-metabolism
E2 benchmark**, not clinical plasma PK, PopPK or a dosing model.

Official PMC AWS source files are SHA-locked:

- supporting PDF: `5beaab1117a32bb281734504dc19237bd264b5db0857837dfd794a9319892bab`;
- Data S1 proteomics: `8ce5df4f3e167915bbb963c8706b747360cf8ea9dd9538e7184d2a454d9c6ea6`;
- Data S2 bioinformatics: `623b6f4813193c348aa4a213ea9c4f5d380230a8932145cab472420059b3b213`.

No donor identifier, individual protein value, CLint value or prediction is
stored in GitHub artifacts or committed to the repository.

## Source Table S3 is now qualified

The 21-page supporting PDF is native text; OCR or figure digitisation is not
required. A strict source parser identifies Table S3 over pages 5–7:

- **75 source rows = 15 donors × 5 probe drugs**;
- the exact same 15 donor/batch identifiers occur in every drug block;
- drugs: midazolam, omeprazole, diclofenac, bupropion and bufuralol;
- columns contain corresponding CYP concentrations in HLM/HH/HL, CYP amount in
  incubation, HH Kpuu and five equation-derived clearance columns;
- `ND` is retained as source missingness rather than converted to zero.

The equation-column mapping is independently checked against the paper:

- `eq.1` = source HLM `CLint,mic`, mL/min/kg body weight;
- `eq.4` = source HH `CLint,hep`, mL/min/kg body weight.

For midazolam, parsed aggregate `eq.1` is **374.7 (203.4–818.9)** and
`eq.4` is **134.0 (45.5–413.1)**, reproducing the rounded values reported in
the article. Bufuralol has exactly **2 ND HLM CLint** rows, also matching the
published description.

Therefore the former status
`donor_level_clint_table_verified = false` is obsolete.

## Prespecified held-donor benchmark

For each drug and assay system separately:

1. target = log source CLint;
2. baseline = training-donor mean log CLint;
3. biological predictor = **one corresponding CYP protein concentration**
   measured in the same system;
4. fixed Ridge alpha = 1 after train-only standardisation;
5. leave-one-donor-out;
6. no variable selection or hyperparameter search;
7. paired donor bootstrap of the loss difference is descriptive only and is
   not multiplicity-adjusted.

Primary metric below is
`Δ log-MAE = MAE(CYP protein) − MAE(training mean)`; negative is better.

### Human hepatocytes (HH)

| Probe drug | CYP | n | Baseline log-MAE | CYP protein log-MAE | Δ log-MAE | Descriptive 95% paired bootstrap |
|---|---|---:|---:|---:|---:|---:|
| midazolam | CYP3A | 15 | 0.4350 | **0.1626** | **−0.2724** | **[−0.4603, −0.1115]** |
| omeprazole | CYP2C19 | 15 | 0.9011 | **0.4289** | **−0.4721** | **[−0.7624, −0.1496]** |
| diclofenac | CYP2C9 | 15 | 0.2344 | **0.1439** | **−0.0905** | **[−0.1624, −0.0110]** |
| bupropion | CYP2B6 | 15 | 0.2826 | 0.2702 | −0.0124 | [−0.0589, 0.0322] |
| bufuralol | CYP2D6 | 15 | 0.3405 | 0.3751 | +0.0347 | [−0.0024, 0.0794] |

Descriptive pooled log-MAE across the repeated drug observations is
**0.4387 → 0.2762** (Δ −0.1626). This pooled number is **not an independent
75-person analysis**: the same 15 donors contribute up to five drugs.

### Human liver microsomes (HLM)

| Probe drug | CYP | n | Baseline log-MAE | CYP protein log-MAE | Δ log-MAE | Descriptive 95% paired bootstrap |
|---|---|---:|---:|---:|---:|---:|
| midazolam | CYP3A | 15 | 0.3390 | 0.2199 | −0.1191 | [−0.2571, 0.0116] |
| omeprazole | CYP2C19 | 15 | 0.7936 | 0.7673 | −0.0263 | [−0.1704, 0.1648] |
| diclofenac | CYP2C9 | 15 | 0.4812 | 0.4882 | +0.0070 | [−0.0369, 0.0572] |
| bupropion | CYP2B6 | 15 | 0.4866 | 0.5184 | +0.0317 | [−0.0435, 0.1181] |
| bufuralol | CYP2D6 | 13 | 0.4172 | 0.4932 | +0.0760 | [−0.0073, 0.2142] |

Descriptive pooled log-MAE is **0.5059 → 0.4975** (Δ −0.0084). None of the
five per-drug descriptive bootstrap intervals excludes zero.

## Scientific interpretation

This is a genuine same-donor demonstration that measured molecular protein
abundance can carry predictive information for measured drug metabolism.
However the effect is **system- and substrate-dependent**.

The strongest signal occurs in intact hepatocytes for midazolam,
omeprazole and diclofenac. It is weak or absent for other substrate/system
pairs, and the HLM pooled improvement is negligible. Protein abundance is
therefore **not a universal surrogate for CLint**.

This result supports Omics2PK's central concept of testing a functional
molecular phenotype against a pharmacological outcome, while also showing why
the platform must preserve negative and heterogeneous results.

It does **not** establish:

- RNA superiority or RNA+protein benefit — Wegler has no qualified matched RNA;
- systemic patient clearance, AUC or concentration prediction;
- external donor/site validation;
- causal mediation;
- a clinical dose recommendation.

## Remaining optional Wegler extension

Data S1 contains the global HH/HLM proteomic matrices. A future high-dimensional
or multi-protein benchmark would require a separately prespecified donor-column
mapping and nested feature handling. It is **not required** for the primary
targeted-CYP result above and should not replace it with a p≫n search on 15
donors.
