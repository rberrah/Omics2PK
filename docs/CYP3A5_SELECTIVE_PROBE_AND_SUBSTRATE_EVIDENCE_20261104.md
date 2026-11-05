# Public same-donor CYP3A5 selectivity and substrate-specific metabolism evidence (2026-11-04)

## Scope

This note separates: **genotype → molecular expression → isoform-selective catalytic phenotype**, **shared CYP3A catalytic phenotype**, and **drug-specific metabolic outcome**. A drug's *in-vitro* metabolite/depletion signal after drug incubation is NOT independently measured baseline endogenous metabolomics. True tacrolimus PK remains untested.

## Newly qualified sources

| Study | DOI / primary source | Real study components | Use and fundamental limitation |
|---|---|---|---|
| Li et al. 2014 | [10.1124/dmd.113.054726](https://pmc.ncbi.nlm.nih.gov/articles/PMC3935135/) | T-5 N-oxidation identified as a **highly CYP3A5-selective** probe reaction (reported >100-fold selectivity over CYP3A4) in recombinant systems and liver preparations | Functional E2 enzyme-specific **reference assay**; selectivity is NOT absolute specificity; the measured T-5 reaction is an outcome, not independent predose metabolomics |
| Bissada et al. 2019 | [10.1124/dmd.119.088823](https://pmc.ncbi.nlm.nih.gov/articles/PMC6815940/) | n12 genotyped HLM, n15 separately genotyped single-donor hepatocytes; **both** midazolam 1′-hydroxylation (shared CYP3A) and **T-5 N-oxidation** (CYP3A5-selective); lapatinib M1 and reactive-metabolite assays | Ideal **mechanistic** counterexample; original per-donor activity numeric **matrix not yet directly recovered**, so do not estimate patient-level heldout MAE from published group summaries |
| Topletz et al. 2013 | [10.1124/dmd.113.051094](https://pmc.ncbi.nlm.nih.gov/articles/PMC3876810/) | **20 individual HLM** Table 2, CYP3A5 genotype (*3/*3, *1/*3, *1/*6, *1/*7 etc), quantified CYP3A4 and CYP3A5 proteins, **observed** vinorelbine metabolism rate, published derived CLint/CLh | Ideal **substrate-specific negative control**: high CYP3A5 expression is NOT sufficient to imply higher vinorelbine microsomal clearance; calculated CLint and CLh are derived from rate/assumptions and are NOT independent outcomes. Still needs a DOI-pinned automated original-source replay |
| Grangeon et al. 2025 | [10.3390/ph18121789](https://pmc.ncbi.nlm.nih.gov/articles/PMC12735975/) | Original DOI-pinned EuropePMC JATS audit retrieved tables 1–3. Table 1 includes 24 individual HLM IDs and separate pools; CYP3A5 genotypes reported for **HLM1–HLM16** in supplementary **Table S1**, while Table 2 reports large-scale protein quantities and additional figures show enzyme activities. | Genotype/protein/activity **per-person numeric crosswalk remains UNVERIFIED** until source PDF supplement is inspected, exact HLM ID matching is proved and numeric enzyme rate is not inferred from correlation figures |

## Bissada mechanistic negative control

In **n15** independent donor hepatocytes, published T-5 N-oxide activity (pmol/min/million cells):
- CYP3A5*3/*3, n8: **1.7 ± 0.7**;
- CYP3A5*1/*3, n5: **16.4 ± 9.3**;
- CYP3A5*1/*1, n2: **9.4 ± 1.5**.

The heterozygote mean exceeds the homozygous *1/*1 mean in this **tiny sample**; do NOT force a strict per-genotype rank or interpret these as universal percent activity scores. In the same n15 hepatocytes, lapatinib M1 formation correlated with **shared CYP3A midazolam** phenotype (published R²=0.87), NOT with CYP3A5-selective T-5 activity (R²=0.07). This is critical evidence that drug-specific CYP3A4 vs CYP3A5 contribution differs from enzyme-specific activity and why metabolite formation cannot be silently treated as a CYP3A5 activity biomarker. None of these values are original per-person vectors for a fresh supervised CV.

## Topletz genotype/protein discordance

Original Table 2 reports CYP3A5*1/*3 donor **IUL-57** among the *low measured expression* group (CYP3A5 protein 0.8 pmol/mg), and *1/*7 donor **IUL-71** (1.1 pmol/mg). This is not an excuse to overwrite CPIC's qualitative allele function assignments: alleles describe a biological genetic expectation, not a fixed observed abundance per tissue or per donor. It is equally invalid to classify *1/*3 as a no-function *diplotype* merely because a paper's measured protein was low.

**Prespecified retrospective Topletz comparator**, once DOI-locked original rows are machine-readable:
- observed primary outcome = **vinorelbine rate at published substrate concentration**, not derived physiological CLh;
- identical 20 real HLM donors and LODO groups;
- models: intercept → CYP3A5 PGx alone → CYP3A4 protein alone → CYP3A4 protein + PGx → CYP3A4+CYP3A5 protein (+/- PGx);
- fixed Ridge alpha1, train-only normalization, no selected features or retrospective graph tuning;
- negative mechanism control: CYP3A5 protein or PGx should not be assigned an effect **by construction**, and additional features may worsen MAE; measured source units must be explicit;
- *1/*6 and *1/*7 are compound alleles that require correct functional interpretation rather than simply assigning number of `*1` without haplotype consideration.

Published interpretation remains that vinorelbine CL is mainly associated with CYP3A4 phenotype. This is not proof for other CYP3A substrates.

## Exploratory vinorelbine negative-control computation (20 original article rows, MANUALLY TRANSCRIBED)

The 20 original Topletz Table-2 rows were **manually transcribed in an ephemeral local analysis** from the publicly printed article. Same 20 source donor identities and exactly one held-out donor for every arm. Each fitted model uses fixed Ridge(alpha=1) on log observed `v_[1 µM vinorelbine]` with strictly train-only z-score standardization. The unfitted comparator is the mean observed donor rate in the other 19 donors. The table is original human liver research, but transcription was **NOT independently machine-replayed from the DOI-pinned article**; metrics are therefore preliminary and should not be claimed as an independently reproduced public dataset benchmark.

| Predeclared exploratory arm | LODO MAE (pmol/min/mg measured rate) |
|---|---:|
| Other donors' mean | 16.412632 |
| CYP3A5*1 copy count only | 16.647751 |
| CYP3A5 protein only | 16.475507 |
| CYP3A4 protein only | **8.767444** |
| CYP3A4 protein + CYP3A5*1 copy count | 9.199472 |
| CYP3A4 + CYP3A5 protein | 9.038718 |
| CYP3A4 + CYP3A5 proteins + CYP3A5*1 copy count | 9.254531 |

These results are an **exploratory negative control on the source study**: adding CYP3A5-related variables did not improve the heldout MAE beyond CYP3A4 protein alone in these 20 microsomal donors. This is entirely consistent with the original authors' conclusion that CYP3A4 dominated vinorelbine metabolism, but cannot establish causality, extrapolate to clinical vinorelbine exposure or infer anything about tacrolimus clearance. No CYP3A5-specific activity was independently measured in this table, and the reported derived CLint/CLh must not be treated as additional outcomes.

A reusable source-qualification and donor-heldout model `src/omics2pk/topletz_vinorelbine_source_gate.py` was added, with **synthetic** contract tests. It requires the caller to supply the authentic original 20 source rows with DOI and exact unique donor keys; raw source table and donor records are not committed. This function does **not** independently authenticate the original source HTML, which remains inaccessible to the current automated acquisition route.

## Source acquisition status (automated audit)

The DOI-locked audit `scripts/audit_public_cyp3a5_followup_sources.py` retrieves the original articles from Europe PMC JATS first and PMC HTML second, and reports only table structure/supplement attachment metadata and SHA256. On 2026-11-04 original JATS for **Grangeon succeeded** (tables 1–3). At that time the public Bissada and Topletz endpoints were unreachable/structure-blocked by the current runner: **ABSTAIN_OFFICIAL_ARTICLE_UNAVAILABLE** is an *access state*, not proof article data do not exist. The inventory does NOT verify Grangeon supplementary donor values and does NOT perform numerical CYP3A5-vs-PK prediction. Explicitly keep the abstention rather than using web-rendered article rows or figures as if the original table had been independently parsed.

## Implemented assay-identity validation

`src/omics2pk/cyp3a_activity_phenotype.py` records CYP3A5-selective T-5, shared CYP3A midazolam/testosterone probes, endogenous 4β-hydroxycholesterol, and post-incubation lapatinib/vinorelbine outcomes with explicit primary references. `benchmark_cyp3a5_activity` now **rejects original-source CYP3A5 activity inputs without a validated source-attested probe code and its original assay DOI**. Shared midazolam, drug metabolite and unspecified activity cannot serve as a CYP3A5-selective E2 label. This is a conservative software gate, not independent wet-lab assay qualification.

No real human original CYP3A5-selective activity + RNA + protein + PGx + baseline metabolomics cohort was downloaded or numerically benchmarked. Cohorts remain small and molecular/PK links incomplete.
