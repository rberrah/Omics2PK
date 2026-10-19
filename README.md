# Omics2PK

**General-purpose research framework for evaluating how transcriptomic, proteomic and metabolomic states improve individual pharmacokinetic prediction beyond clinical information and pharmacogenetics.** Tacrolimus/CYP3A is a current flagship validation case, **not a restriction of the framework to one drug or pathway**.

![General Omics2PK graphical abstract: measured clinical/genetic and single/multiomics data pass through direct-feature, functional-score, RNA-network or latent representations, compared with regression, XGBoost and PopPK on held-out patients and genuine measured PK.](docs/assets/omics2pk-generated-optimized.webp)

[Original generated image (PNG, full resolution)](docs/assets/omics2pk-generated-original.png) · [Optimized WebP](docs/assets/omics2pk-generated-optimized.webp) · [Editable vector adaptation (SVG)](docs/assets/omics2pk-graphical-abstract.svg) · [Scientific tacrolimus protocol](docs/TACROLIMUS_FIRST_CYP3A4_CYP3A5_RESEARCH_PROTOCOL_20261019.md) · [Real-cohort source audit](docs/CYP3A_TACROLIMUS_OMICS_PK_PUBLIC_SOURCE_ACQUISITION_AUDIT_20261019.md) · [MultiOmics producer](https://github.com/rberrah/MultiOmics)

> **Status (19 October 2026): exploratory, source-qualified research prototype; not clinical software.** The data-ingestion and many research comparators work, but **no independently validated gain in tacrolimus PK prediction over clinical + authentic CYP3A4/5 genotype has yet been demonstrated**. Omics2PK is **not a medical device**, a validated dose recommendation system or an automatically fitted clinical PopPK engine.

## The question

For the **same people, measured PK outcomes, prediction landmark and held-out folds**, do individual molecular states improve prediction of drug AUC, CL/F or observed concentrations **beyond the clinical and authentic genotype reference appropriate to that drug**? Tacrolimus with CYP3A4/CYP3A5 is the initial flagship, but the contract permits other medicines and relevant CYP, UGT, transporter or biological pathways.

The comparison is designed to be falsifiable: an omic signal may be redundant, unstable, tissue-inappropriate, or harmful for predictive performance. Links between RNA, enzyme activity and ex-vivo intrinsic clearance are **not** equivalent to predicting in-vivo tacrolimus disposition.

### Research design

| Comparison axis | Candidate arms | Primary control |
|---|---|---|
| **Source information** | Clinical only; clinical + authentic PGx; + one RNA / protein / metabolite layer; + matched multi-omic layers | Clinical + PGx, on the same eligible people |
| **Representation** | Prespecified CYP3A functional score; direct low-dimensional values; fold-fitted PCA; latent factors; experimental inductive RNA network | Unsupervised PCA and fixed low-dimensional score |
| **Estimator** | Penalized regression; XGBoost where sample size supports it; eventual dose–concentration PopPK likelihood | Identical training/held-out subjects and drug PK endpoints |
| **Drug-specific processes** | CYP3A4/5 metabolism; intestinal first pass; ABCB1 transport; whole-blood partition; post-transplant state | Prespecified, measured or explicitly unavailable covariates |

**No one omic layer, SWEET graph or ML method is assumed superior.** Original SWEET concerns RNA networks; MS1 peaks must not be relabelled as RNA features. Train/test separation includes all molecular preprocessing, imputation, signature weights, graph reference fitting, XGBoost tuning and latent factors. A true PopPK comparison requires original doses, times and concentration observations; regressing an already estimated CL/F is only a preliminary benchmark.

### Original-cohort PK/PGx readiness — available without model fitting

Before searching for an apparent omic effect, assess the original biological
person/sample overlap and clinical PK/PGx availability. The new standard-library
`omics2pk audit-pk-readiness /path/to/exchange` command checks a genuine
MultiOmics V1 source and reports **unique subjects versus assay rows**, one
specimen per modality, repeated samples, and all mono-/multiomic intersections.
It deliberately reports **unknown** clinical PK/genotype eligibility if the
original clinical crosswalk was not provided. With authorized private in-memory
availability metadata, the API also checks **subject + original assay +
specimen** and *original* pre-dose collection clocks before deeming an
individual potentially eligible for PK benchmarking. It never exposes
individual measurements or identifiers.

The [BIOMARGIN–ABIS readiness protocol](docs/MULTIOMICS_PK_READINESS_AND_BIOMARGIN_ABIS_JOIN_20261019.md)
specifies the minimal source-key and availability counts to request before
patient-level transfer, and why a large public RNA dataset by itself does not
establish drug PK predictive benefit.

### New fair single- versus multi-omic benchmark (research, synthetic tests)

The library now provides a strictly donor-/patient-held-out **five-fold** ablation that compares a clinical reference and genuine PGx, individual modality PCA2 and direct Ridge, optional **externally fixed biological scores**, unsupervised PCA **early fusion** and a **fixed late fusion** of single-omic predictions. Optional XGBoost on low-dimensional projections is intentionally blocked when fewer than 80 original people are available. **No ranking or performance claim has been established on a new original clinical cohort by the software addition alone.**

For private authorized originals, `src/omics2pk/multimodal_exchange_pk.py` consumes versioned and SHA-checked source `DOWNSTREAM_EXCHANGE_V1` matrices plus a separately supplied **in-memory person→assay→collection time→PK/PGx** sidecar, with explicit complete-case eligibility. It refuses row-order joins, repeated participant assays, source-unknown patients and invented `T0` timestamps. The benchmark `src/omics2pk/multimodal_pk_ablation.py` fits modality transforms only inside outer training folds and emits **only aggregate error metrics**, with descriptive paired bootstrap intervals explicitly lacking model refitting. Requires `pip install -e '.[benchmark]'` and suitable original measured PK outcomes; no clinical/PK data belong in the public repository.

See [matched-modalities research protocol](docs/MULTIMODAL_PK_ABLATION_SOFTWARE_QUALIFICATION_20261019.md) for API, assumptions and leakage guards.

## What is already established?

| Evidence gate | Original observations and result | What it does **not** establish |
|---|---|---|
| **Technical exchange** | Real-source-compatible `MultiOmics → DOWNSTREAM_EXCHANGE_V1 → Omics2PK` interface with SHA-256, sample/assay identity checks, private-source support and synthetic CI | End-to-end public-browser export parity in every setting, or biological validity of scores |
| **Enzyme phenotype (E2)** | Genuine donor-matched human hepatic RNA → measured CYP activity; positive gains for several enzymes in independent donor-heldout tests | In-vivo tacrolimus PK or superiority over complete patient PGx |
| **Measured clinical CYP3A probe PK (E3)** | Original Jarmusch plasma n=13 at predose with CYP3A4/5 and **midazolam** CL/F: clinical + PGx MAE **7.343** vs + global MS1 PCA2 **9.888** on day 1 (worse); day 9 n=12 **4.601 → 4.749** (worse) | Tacrolimus efficacy, omic superiority, independent cohort validation |
| **Another real clinical PK cohort (E3)** | Rotroff metformin n=32, original predose metabolomics + measured AUC/renal CL: global PCA offered no convincing incremental gain | CYP3A/tacrolimus response |
| **Transplant multi-view omics** | BIOMARGIN/ROCKET original author notebook: **131** unique people, six RNA/miRNA views (35 six-view complete); public GEO blood RNA-seq **58,828 genes × 384 assays** confirmed | Six distinct omic modalities; GEO ↔ BIOS ↔ ABIS person crosswalk or matched original TAC PK/PGx |

Source reports: [real midazolam pilot](docs/JARMUSCH14_ORIGINAL_GNPS_BIOM_TO_REAL_MIDAZOLAM_PK_20261019.md) · [private MultiOmics numerical parity](docs/JARMUSCH13_PRIVATE_MULTIOMICS_NATIVE_EXPORT_AND_PK_PARITY_20261019.md) · [metformin](docs/E3_REAL_METFORMIN33_PREDOSE_384_OMICS_AUC_CLR_HELDOUT_20261018.md) · [BIOMARGIN public original-data audit](docs/BIOMARGIN_131_REAL_PUBLIC_COUNTS_AND_GEO_BIOS_CROSSWALK_20261019.md).

The exploratory midazolam errors are in original source numerical units not yet independently established as physical CL/F units. These small cohorts do not support comparative XGBoost versus SWEET versus multi-omic claims.

### The current *scientific* data bottleneck

To demonstrate the central tacrolimus claim we still need, **within authorized original patients**:

1. Original molecular measurements **before the PK prediction landmark** and exact person/specimen/visit identity. For longitudinally treated graft recipients, distinguish contemporaneous state from truly pre-treatment omics.
2. Authentic **CYP3A5/CYP3A4** variants (with allele coverage and missingness), clinical covariates, formulation, dose history, exact concentration sampling times, hematocrit and time after transplant.
3. Measured tacrolimus PK (dose-adjusted trough and/or dense AUC / concentration–time profiles), external or center-disjoint validation and permission for a research-only sample-level join.

**Candidate data requests:** BIOMARGIN ↔ ABIS (the original **GEO title → BIOS → ABIS** crosswalk is presently missing), Vanhove n147 tacrolimus + CYP3A phenotype/PGx (full person table not publicly released), Zhu n137 longitudinal tacrolimus + metabolomics (original matrix and genotypes missing), Phapale n29 tacrolimus predose metabolomics (complete numeric individual arrays missing). See the [qualified real source catalogue](docs/CYP3A_TACROLIMUS_OMICS_PK_PUBLIC_SOURCE_ACQUISITION_AUDIT_20261019.md). MIMIC-IV is suitable for certain clinical PK workflows but **does not provide original patient-linked omics in its standard released schema**.

## Data flow and repository boundaries

```text
MultiOmics (independent repository)
  └─ real individual assays + metadata + QC / explicitly qualified transforms
       └─ source exchange V1 (manifest.json, samples.csv, features.csv, matrices/)
           └─ Omics2PK: strict source/person/time/units gates
               ├─ clinical + genuine PGx reference
               ├─ + single-omic direct / prespecified functional score
               ├─ + matched multiomic factors or RNA network (only when supported)
               └─ held-patient PK comparison → optional future PopPK integration
                    └─ Omics2PD: future separately qualified project
```

No patient-level original material belongs in this public repository, its PRs, Actions, artifacts or documentation. MultiOmics, PMx Explain and Omics2PK are separate projects; this repository does **not** vendor their applications. Protected thesis Pillar A/B/C/D implementations are **not** published here.

## Installation and basic use

Requires **Python ≥3.11**. Standard source-exchange validation uses the Python standard library.

```bash
git clone https://github.com/rberrah/Omics2PK.git
cd Omics2PK
python -m pip install -e .
python -m unittest discover -s tests -v
omics2pk validate-exchange /path/to/extracted/exchange
```

`validate-exchange` checks hashes, sample IDs, matrix schema and numeric content; it is **not** a PK model fit. Use only data you have permission to process. Do not silently infer missing samples or clinical timestamps.

For qualified **normalized individual transcriptomic expression** (not raw RNA counts), an optional exploratory SWEET preparation path is:

```bash
omics2pk prepare-sweet /path/to/exchange /private/sweet_inputs
```

For optional Python research estimators:

```bash
python -m pip install -e '.[benchmark]'
```

This installs benchmarking dependencies; it does **not** unlock source-restricted datasets or guarantee biomedical superiority. The real Jarmusch/BIOMARGIN source qualification scripts require their original authorized sources or retrieve verified public inputs ephemerally. See the dated study reports before running them. The RNA-only SWEET-inspired experimental network is **not** a verified reproduction of original SWEET.

## Reproducibility, release and use boundaries

- **Versioning:** public evidence scripts, sample-level contract checks and synthetic CI; source-specific files have SHA-256 and documented source limitations.
- **Scientific validation:** positive and negative experiments are retained. Technical green CI, an exploratory metric or an in-vitro phenotype is **not** a validated PK model.
- **Privacy:** preserve institutional consent and data-use restrictions; never commit per-patient original omics, genotypes, PK measurements, consent metadata or linkage keys.
- **Licensing:** [PolyForm Noncommercial 1.0.0](LICENSE); **source-available, not OSI open source**. Third-party sources, data and licensed components retain their own terms. Scientific citation metadata: [CITATION.cff](CITATION.cff).
- **History:** the former chronological README with dated discovery logs is preserved intact in [docs/history](docs/history/README_CHRONOLOGICAL_ARCHIVE_20261019.md), rather than allowing outdated intermediate source-access states to dominate the project home page.

This software is for **research and teaching only**, not clinical dose selection or automated therapeutic recommendations.
