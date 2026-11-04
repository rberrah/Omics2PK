# Omics2PK

**General-purpose research framework for evaluating how transcriptomic, proteomic and metabolomic states improve individual pharmacokinetic prediction beyond clinical information and pharmacogenetics.** Tacrolimus/CYP3A is a current flagship validation case, **not a restriction of the framework to one drug or pathway**.

![General Omics2PK graphical abstract: measured clinical/genetic and single/multiomics data pass through direct-feature, functional-score, RNA-network or latent representations, compared with regression, XGBoost and PopPK on held-out patients and genuine measured PK.](docs/assets/omics2pk-generated-optimized.webp)

[Original generated image (PNG, full resolution)](docs/assets/omics2pk-generated-original.png) · [Optimized WebP](docs/assets/omics2pk-generated-optimized.webp) · [Editable vector adaptation (SVG)](docs/assets/omics2pk-graphical-abstract.svg) · [Scientific tacrolimus protocol](docs/TACROLIMUS_FIRST_CYP3A4_CYP3A5_RESEARCH_PROTOCOL_20261019.md) · [Real-cohort source audit](docs/CYP3A_TACROLIMUS_OMICS_PK_PUBLIC_SOURCE_ACQUISITION_AUDIT_20261019.md) · [MultiOmics producer](https://github.com/rberrah/MultiOmics)

> **Status (22 October 2026): source-qualified research framework with real positive and negative E2/E3 benchmarks; not clinical software.** The repository now includes leakage-controlled source contracts, a qualified experimental inductive RNA-network workflow, real hepatic RNA/protein functional benchmarks and authentic clinical Omics→PK examples. **No independently validated gain in tacrolimus PK prediction over clinical + authentic CYP3A4/5 genotype has yet been demonstrated.** The public **Analyze** interface therefore remains explicitly **in development**. Omics2PK is **not a medical device**, a validated dose recommendation system or an automatically fitted clinical PopPK engine.

## The question

For the **same people, measured PK outcomes, prediction landmark and held-out folds**, do individual molecular states improve prediction of drug AUC, CL/F or observed concentrations **beyond the clinical and authentic genotype reference appropriate to that drug**? Tacrolimus with CYP3A4/CYP3A5 is the initial flagship, but the contract permits other medicines and relevant CYP, UGT, transporter or biological pathways.

The comparison is designed to be falsifiable: an omic signal may be redundant, unstable, tissue-inappropriate, or harmful for predictive performance. Links between RNA, enzyme activity and ex-vivo intrinsic clearance are **not** equivalent to predicting in-vivo tacrolimus disposition.

### Original human-donor source replay: genotype, proteins, enzyme function (2026-11-04)

An [independently reproduced published Latham PHH source pilot](docs/LATHAM_PHH_SOURCE_REPLAY_REAL_20261104.md) now reads **original public Table 4 in memory** (14 donor hepatocyte rows; the same 13 protein-complete donors for all four comparison arms). It verifies CYP3A5 *3/*6/*7 genotype classes and compares CYP2C8+CYP3A5 PGx against additional measured CYP proteins and CYP2C8/shared-CYP3A probe activities in fixed donor-heldout Ridge models. Replayed imatinib **in-vitro** source-table MAE: PGx 2.906576, +proteins 3.003704, +probe activities 2.762944, +both 3.262812; [original-source GitHub Actions replay](https://github.com/rberrah/Omics2PK/actions/runs/41804929408) passed. Among CPIC-classified nonexpressers, 9/9 with quantitative CYP3A5 peptide data had a positive signal and 10/10 had positive *shared* CYP3A probe activity—**neither establishes CYP3A5-specific catalytic function**. These are very small in-vitro original data and do **not** demonstrate clinical PK gains, multi-omic benefit, or tacrolimus dose utility. Original donor-level data never enter this repository.

### CYP3A5-selective versus shared CYP3A activities (source audit, 2026-11-04)

A [new published human-evidence and assay-specificity report](docs/CYP3A5_SELECTIVE_PROBE_AND_SUBSTRATE_EVIDENCE_20261104.md) distinguishes validated selective T-5 N-oxidation from shared CYP3A midazolam/testosterone activity, endogenous 4β-hydroxycholesterol and **post-drug-incubation** metabolites. The original-source activity benchmark now refuses an unqualified or shared CYP3A probe as a **CYP3A5-selective** training target. Literature controls include genotyped livers with drug-substrate-specific metabolic outcomes (Topletz vinorelbine, Bissada lapatinib) and Grangeon human-liver protein/genotype source inventory. **No same-donor complete PGx+RNA+protein+metabolomics+selective activity+clinical PK dataset is yet verified.** Live external-source audits can be abstained by remote host availability; offline contracts are the mandatory PR CI.

### Original GEO liver miRNA–CYP3A activity: acquisition PASS, biological source-parity ABSTAIN (2026-11-04)

The [GSE49012 original GEO source-qualification report](docs/GSE49012_ORIGINAL_MIRNA_CYP3A_SOURCE_PARITY_20261104.md) documents **34 same-liver original public microRNA arrays and source-declared CYP3A activities**: 22 cirrhotic and 12 normal, 683 biological miRNA probes, identical subject-heldout partitions and strictly train-only molecular PCA. Although the original GEO sample activities can be reproduced, their **22-cirrhotic-donor activity sum (12.533) is incompatible with the publication's whole n28 cirrhotic activity mean (0.28)** under the same-endpoint-and-subset assumptions. The normal n12 GEO mean (0.4708) matches the published 0.47. Thus any computed miRNA MAEs are **source-data diagnostics only**, not biological performance evidence; no benefit/negative omics claim is qualified. An explicit source-checking gate reports `ABSTAIN_GEO_ACTIVITY_PUBLICATION_GROUP_PARITY_FAILURE`. Provenance resolution is tracked in [issue #138](https://github.com/rberrah/Omics2PK/issues/138). The external source replay can be invoked separately; mandatory offline source tests remain in CI.

### Biological knowledge reconciliation (research-only prototype, 2026-11-03)

A [CYP3A5 functional-evidence pilot](docs/CYP3A5_FUNCTIONAL_EVIDENCE_PILOT_20261103.md) and [held-donor phenotype-proximity protocol](docs/PHENOTYPE_PROXIMITY_AND_FUNCTIONAL_ACTIVITY_20261103.md) and [split-conformal interval/assay-leakage qualification](docs/CYP3A5_ACTIVITY_UNCERTAINTY_AND_LEAKAGE_20261103.md) now triage potential contradictions between qualified *1/*3/*6/*7 allele calls and transcript/protein/activity/metabolite observations. Total-gene RNA detection is never treated as proof of active protein; ambiguous metabolites and incomplete PGx cause explicit warnings/abstentions. The module assigns **no fixed percent-of-function weights**, does **not** fit a PK model, and has **no demonstrated predictive gain**. MultiOmics V1 producer artifacts remain unchanged. Synthetic rule checks run in the offline Python tests.

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

**Multicenter validation:** an additional `validation_mode="leave_one_site_out"` option holds out each complete original research center in turn (source-attested `site_ids` or clinical-sidecar `site_id`). Site labels cannot silently enter patient-random folds. This is still **within-parent-study site-disjoint validation**, not an independent external cohort.

For private authorized originals, `src/omics2pk/multimodal_exchange_pk.py` consumes versioned and SHA-checked source `DOWNSTREAM_EXCHANGE_V1` matrices plus a separately supplied **in-memory person→assay→collection time→PK/PGx** sidecar, with explicit complete-case eligibility. It refuses row-order joins, repeated participant assays, source-unknown patients and invented `T0` timestamps. The benchmark `src/omics2pk/multimodal_pk_ablation.py` fits modality transforms only inside outer training folds and emits **only aggregate error metrics**, with descriptive paired bootstrap intervals explicitly lacking model refitting. Requires `pip install -e '.[benchmark]'` and suitable original measured PK outcomes; no clinical/PK data belong in the public repository.

See [matched-modalities research protocol](docs/MULTIMODAL_PK_ABLATION_SOFTWARE_QUALIFICATION_20261019.md) for API, assumptions and leakage guards.

## What is already established?

| Evidence gate | Original observations and result | What it does **not** establish |
|---|---|---|
| **Technical exchange** | Real-source-compatible `MultiOmics → DOWNSTREAM_EXCHANGE_V1 → Omics2PK` interface with SHA-256, sample/assay identity checks, private-source support and synthetic CI | End-to-end public-browser export parity in every setting, or biological validity of scores |
| **RNA → enzyme phenotype (E2)** | Genuine U.S. EPA n≈100 donor-matched human liver RNA → measured CYP activity; several enzymes improve held-donor prediction, with negative CYP2E1 retained | In-vivo tacrolimus PK, protein benefit or superiority over complete patient PGx |
| **Protein → measured CLint (E2)** | Wegler n=15 same liver donors × five probe drugs: targeted CYP protein improves held-donor HH CLint prediction strongly for midazolam/omeprazole/diclofenac, while HLM gain is nearly null overall | RNA-vs-protein superiority, systemic clearance, clinical dose benefit or external validation |
| **Measured clinical CYP3A probe PK (E3)** | Original Jarmusch plasma n=13 at predose with CYP3A4/5 and **midazolam** CL/F: clinical + PGx MAE **7.343** vs + global MS1 PCA2 **9.888** on day 1 (worse); day 9 n=12 **4.601 → 4.749** (worse) | Tacrolimus efficacy, omic superiority, independent cohort validation |
| **Another real clinical PK cohort (E3)** | Rotroff metformin n=32, original predose metabolomics + measured AUC/renal CL: global PCA offered no convincing incremental gain | CYP3A/tacrolimus response |
| **Transplant multi-view omics** | BIOMARGIN/ROCKET original author notebook: **131** unique people, six RNA/miRNA views (35 six-view complete); public GEO blood RNA-seq **58,828 genes × 384 assays** confirmed | Six distinct omic modalities; GEO ↔ BIOS ↔ ABIS person crosswalk or matched original TAC PK/PGx |

Source reports: [Wegler protein→CLint](docs/WEGLER_SOURCE_STRUCTURE_RESULTS_20261016.md) · [real midazolam pilot](docs/JARMUSCH14_ORIGINAL_GNPS_BIOM_TO_REAL_MIDAZOLAM_PK_20261019.md) · [private MultiOmics numerical parity](docs/JARMUSCH13_PRIVATE_MULTIOMICS_NATIVE_EXPORT_AND_PK_PARITY_20261019.md) · [metformin](docs/E3_REAL_METFORMIN33_PREDOSE_384_OMICS_AUC_CLR_HELDOUT_20261018.md) · [BIOMARGIN public original-data audit](docs/BIOMARGIN_131_REAL_PUBLIC_COUNTS_AND_GEO_BIOS_CROSSWALK_20261019.md).

The exploratory midazolam errors are in original source numerical units not yet independently established as physical CL/F units. These small cohorts do not support broad comparative XGBoost-versus-SWEET-versus-multi-omic claims. Wegler is deliberately kept at the **in-vitro E2** level: its heterogeneous HH/HLM result is evidence that a measured molecular phenotype can matter, not evidence that protein universally predicts patient PK.

### The current *scientific* data bottleneck

The remaining high-priority blockers are now primarily **source acquisition or authorized linkage**, not missing generic modeling code:

- **#63:** authenticated Zhang n100 or Ohtsuki n17 same-donor RNA + protein + CYP-activity matrix;
- **#50:** original rosuvastatin n40 processed author CSVs with true AUC/Cmax ↔ omics subject crosswalk;
- **#31:** original Phapale n29 / He n109 patient-level molecular ↔ tacrolimus PK matrices;
- **#14:** lawful deterministic BIOMARGIN ↔ ABIS patient/date/dose/concentration crosswalk.

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

For qualified **normalized individual transcriptomic expression** (not raw RNA counts), an optional exploratory RNA-network preparation path is:

```bash
omics2pk prepare-sweet /path/to/exchange /private/sweet_inputs
```

For optional Python research estimators:

```bash
python -m pip install -e '.[benchmark]'
```

This installs benchmarking dependencies; it does **not** unlock source-restricted datasets or guarantee biomedical superiority. The real Jarmusch/BIOMARGIN source qualification scripts require their original authorized sources or retrieve verified public inputs ephemerally. See the dated study reports before running them.

The RNA-only network remains explicitly **SWEET-inspired and inductive**, not a reproduction of original transductive SWEET. The historical upstream z-score archive mismatch is now explained and current-source normalization is numerically qualified. On authentic GSE9588, train-only per-probe z-scoring was leakage-safe but increased sample-weight clamps in all five deterministic n=150 reference sensitivities, so source-normalized values remain the prespecified representation for that already-normalized array source rather than automatically selecting a new preprocessing transform.

## Reproducibility, release and use boundaries

- **Versioning:** public evidence scripts, sample-level contract checks and synthetic CI; source-specific files have SHA-256 and documented source limitations.
- **Scientific validation:** positive and negative experiments are retained. Technical green CI, an exploratory metric or an in-vitro phenotype is **not** a validated PK model.
- **Privacy:** preserve institutional consent and data-use restrictions; never commit per-patient original omics, genotypes, PK measurements, consent metadata or linkage keys.
- **Licensing:** [PolyForm Noncommercial 1.0.0](LICENSE); **source-available, not OSI open source**. Third-party sources, data and licensed components retain their own terms. Scientific citation metadata: [CITATION.cff](CITATION.cff).
- **History:** the former chronological README with dated discovery logs is preserved intact in [docs/history](docs/history/README_CHRONOLOGICAL_ARCHIVE_20261019.md), rather than allowing outdated intermediate source-access states to dominate the project home page.

## Public website boundary

The public site follows the same navigation thread as MultiOmics: **Apprendre / Analyze / Écosystème / Ressources**. Under **Apprendre**, **Objectif** explains *why* Omics2PK is scientifically relevant and **Démo guidée** teaches *how* to design the comparison. **Analyze remains explicitly “en cours”**: there is no public patient-data upload or pretend one-click clinical model while the external data gates above remain unresolved.

This software is for **research and teaching only**, not clinical dose selection or automated therapeutic recommendations.
