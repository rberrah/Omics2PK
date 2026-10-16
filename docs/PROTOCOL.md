# Omics2PK — scientific protocol v0.1 (exploratory)

## Core hypothesis

Does a drug- and tissue-contextualized patient molecular network add **independent** predictive value for PK/PD beyond baseline clinical/PGx covariates, molecular expression alone, conventional pathway scores and population PK models?

NetPK provides a pharmacological reference network and reported covariate modules. SWEET provides statistical gene co-expression edges **per sample**, not biochemical or causal relations. Never silently merge these edge types.

## Comparison matrix

Evaluate, on **identical subjects** for each paired contrast, the following data availability regimes:
- Clinical covariates only, clinical + DNA genotype, clinical + RNA phenotype, clinical + protein/metabolite phenotype, all available omics. Do not call RNA a genotype.
- For RNA, compare gene expression, PCA, GSVA or similar pathway scores, SWEET connectivity/edge summaries and partial network embeddings.
- For multiomics compare single-view and joint-view methods (e.g. MOFA2), with matched training subjects, missingness and sample count reported.

| Model | Features | Algorithm / purpose |
|---|---|---|
| C0 | Clinical + PGx | Published NetPK or reference NLME comparator |
| M1 | Clinical + preprocessed full omics | XGBoost, tabular ML |
| M1bis | Clinical + SWEET features only | XGBoost, network-only arm |
| M2 | Clinical + omics + SWEET | XGBoost, incremental network effect |
| M3-basic | Clinical + expression/PCA/GSVA/MOFA2 scores | PopPK NLME |
| M3-S | Clinical + SWEET scores **without** Pillar A | PopPK NLME |
| M3-A | Clinical + drug/tissue SWEET scores **with** Pillar A | PopPK NLME, main hypothesis |
| M4 | Clinical + graph kernel/embedding with/without Pillar A | Hybrid PopPK, exploratory |

Keep structural PK model and parameter estimation procedure fixed across NLME representation arms unless structural change is separately prespecified.

## Analysis rules

1. Predefine one measured primary endpoint per dataset. Microsomal CLint is a different task from patient-level apparent CL/F or sparse-sample AUC. Do not combine them as if interchangeable.
2. All subjects from a donor/patient must be kept in a single outer partition. Nested feature selection, normalization, SWEET reference construction, pathway priors, score weights, tuning and model estimation must use **training only**. A truly external cohort stays locked.
3. Distinguish baseline prediction (covariates without patient concentration) from clinical sparse-sampling MAP-BE with C0 or C0+C1+C3; an improved C0-trained XGBoost is not necessarily an improved mechanistic PK model.
4. M3/M4 must use actual dosing and concentration-time observations fitted jointly by a verified NLME engine. Ridge on synthetic individual CL does not count.
5. Test genotype-effect double-counting, normalization/tissue shift, correlated pathways, patient-swapped graphs, randomized networks and irrelevant drug knowledge priors. Report null/negative results.
6. Report held-out subject-level MAE/bias, distribution calibration, uncertainty, complexity/identifiability and uncertainty intervals. An in-sample BIC improvement alone is insufficient.
7. SWEET network transductive behavior and large genome-wide edge counts require special independent assessment; avoid passing a whole dense N-gene squared matrix to SAEM.

## Scientific data priority

- [Official SWEET](https://github.com/SysMednet/SWEET): public algorithm and demo RNA, not observed PK.
- [Yin / Battail, iScience 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12309972/): sample network → treatment response study.
- [Wegler et al. 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8041379/): 15 donor matched proteome + measured liver CLint, not clinical transcriptome.
- BIOMARGIN omics and potentially separate ABIS/ISBA PK require explicit lawful patient/date linkage by P&T before any clinical plan is feasible.

Simulated outcome benchmarks are for software sanity only. Clinical, biological and PK conclusions require real measured endpoints.
