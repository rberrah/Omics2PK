# Omics2PK

**Research-only alpha. Not clinically validated.**

Omics2PK connects the reproducible outputs of [MultiOmics](https://github.com/rberrah/MultiOmics) to research on patient-specific molecular networks (SWEET), drug-specific knowledge from Thesis Pillar A and PopPK/PKPD covariate discovery in Pillar B (route R2).

## New: transition to longitudinal PK (October 2026)

The immediate goal is now **real original patient concentration–time profiles and dosing histories**, followed by **time-varying RNA/protein/metabolite covariates** on matched subjects, rather than additional cross-sectional enzyme correlations. The source-gated PK_LONG_V1 contract verifies observed versus simulated events, predose sampling, BLQ/LOQ, full dose history and prospective molecular sampling timestamps. Original **Theoph**: 12 humans, 12 doses, 132 actual C(t) observations, Cmax/Tmax/AUC0-last descriptive NCA. Separate pinned original RDA audit covers 120 humans, 198 distinct PK occasions and 2,480 measured concentration rows for **mavoglurant**, 32 subjects with warfarin PK/PD mixed records, and explicitly simulated theophylline multidose data. **None of these original PK series contains omics**; they are safe PK architecture/clinical-series qualifications, not proof of transcriptome-driven precision dosing.

See [Longitudinal Omics→PK roadmap](docs/LONGITUDINAL_PK_AND_DYNAMIC_OMICS_ROADMAP_20261016.md) and [original measured PK dataset audit](docs/ORIGINAL_LONGITUDINAL_PK_SOURCE_AUDIT_20261016.md).

## MultiOmics source-interoperability audit: E2 molecular evidence → E3 PK/NCA

We have reproduced the **actual** `rberrah/MultiOmics` V1 browser exporter against its current pinned public source. The native 8-subject, 16-specimen, 49-assay demo **loses 100% of the original categorical T0/T12 visit metadata** on downstream export, because the producer copies `collection_time` but not the input `timepoint`. Raw RNA, unknown provenance scales and missing PK dosing histories remain correctly blocked for clinical interpretation. Omics2PK now offers `python -m scripts.audit_multiomics_pk_readiness <local_exchange_dir>` for locally verifying the original V1 output, plus a strict optional `individual_signatures` contract and measured `Cmax/Tmax/AUC0last` NCA target contract. No attempt is made to derive patient PK from group fold changes, an omics-only ZIP, or unverified chronological labels.

[Full native output audit and producer change request](docs/REAL_MULTIOMICS_OUTPUT_AUDIT_AND_PK_ADAPTER_GAPS_20261016.md). The agreed cross-repo changes are filed as [MultiOmics issue #20](https://github.com/rberrah/MultiOmics/issues/20); this repository does not alter MultiOmics source. Omics→NCA→PopPK/PBPK evidence remains future until genuinely observed, subject-linked data and time-validated covariates exist.

## Public patient-linked omics and PK original-source hunt (2026-10-16)

A first **real individual clinical tacrolimus exposure target** is now available from Wang et al. 2026 PLOS ONE [clinical S1 original supplement](https://doi.org/10.1371/journal.pone.0343817): **36 liver transplant recipients with dose, actual whole-blood C0, C/D, age, sex, height and weight**, including **three aged <18**. The source spreadsheet SHA-256 and original C/D physical dose/weight consistency are checked in CI. A **clinical-only** 36-fold held-patient-out reference (age+sex+BMI fixed ridge vs constant original-source baseline) showed nearly identical MAE; **no omic measurements were joined or fitted**. NCBI SRA metagenome accession PRJNA1416570 is public, but an exact source-defined crosswalk to this clinical table and the original metabolite matrix **remain unverified**.

Original **Fannin 2010 GEO GSE15158 (88 RNA assays)** and **Bushel/Fannin 2016 GSE70784 (432 RNA assays)** were live audited; both have authentic serial blood gene expression but neither currently supplies source-verified individual clinical PK concentration-time tables. Fannin's published early serum APAP graph is not a primary numeric curve; its post-dose RNA is a downstream response, not an admissible baseline predictor of early AUC. Other source candidates include Degraeve et al. 2023 n351 genuine mouse TAC PK blood samples + gut ABCB1 (unverified animal RNA/PK crosswalk), Zhu et al. 2022 n137/693 clinical TAC C0+metabolomics (original patient matrix not openly verified), Zenodo 24-participant 18-visit metabolomics (no original drug PK), plus source-qualified negative HuggingFace training corpora.

**Full audit and strict evidence decisions:** [Public Omics×PK source triage](docs/PUBLIC_SAME_SUBJECT_OMICS_PK_SOURCE_TRIAGE_20261016.md). The repository will not stitch unrelated omics and PK cohorts together or treat a C0/D endpoint as PopPK CL/F/AUC.

## Public PK/omics evidence gate: 93 human TAC microbiomes; parent-metabolite PBPK

In the original [Degraeve et al. 2024 tacrolimus renal cohort](https://doi.org/10.1002/cpt.3077), the genuine ENA `PRJNA929904` source contains **93 distinct original 16S samples** for 93 clinical CYP3A5-nonexpressing transplant recipients. Actual **C0 was observed**, while the paper's AUC0-24 was **machine-learning estimated** from sparse concentrations (with missing-point cohort-median substitution). Original individual C0/AUC ↔ bacterial specimen crosswalk has NOT been obtained, and automated official supplementary ZIP requests returned HTTP 403. Thus these are **not** ready for a patient-validated Omics→AUC study.

Separately, the public [R pksensi APAP dataset](https://search.r-project.org/CRAN/refmans/pksensi/html/APAP.html) contains **32 source rows from three literature studies, including 30 complete parent-acetaminophen / glucuronide / sulfate observations**, with real drug metabolite concentration-time signals. They can support **PBPK multi-analyte population-curve checks**, but there are NO original patient identifiers or independent RNA/metabolomic feature matrices in this RDA.

See [strict original-source audit and E3/PBPK limitations](docs/PUBLIC_93_KIDNEY_TAC_AND_APAP_PBPK_SOURCE_AUDIT_20261016.md). No independent E3 clinical Omics→PopPK gain is claimed.

## What is executable today

- Offline inspection of the **MultiOmics DOWNSTREAM_EXCHANGE_V1** directory: checks SHA-256, study-data access declaration, assay/sample/subject and feature-ID matching, finite numeric matrices and recorded processing scale.
- Conversion of suitable **normalized transcriptomic** matrices into the input files expected by the original SWEET implementation: expression.txt, patient.txt, gene.txt. The converter rejects unaddressed missing data, repeated subjects, incompatible namespaces and mixed tissues.
- Strict quality control of **already-computed SWEET** edge lists and extraction of descriptive network features for the proposed M1bis arm (network-wide edge statistics, per-gene signed/absolute strength). No mechanistic or CYP activity interpretation.
- Experimental training-frozen **SWEET-inspired** network proposal (distinct from upstream SWEET), with subject-level LOO/cross-fitted and held-out features, optionally separating genome-wide sample weights from targeted edge genes. Public 80-sample expression-only diagnostic, but no PK endpoint or validated causal interpretation.
- CI that executes the **real native MultiOmics exporter** on its bundled three-omic synthetic demo: 8 longitudinal subjects, 49 assays, raw RNA strictly blocked from SWEET.
- Research-only **real observed pharmacological benchmark** using 60 NCI-60 cell lines and predeclared doxorubicin GI50, with raw/PCA/network and fixed XGBoost (M1/M1bis/M2), negative shuffled-network control and tumor-type-adjusted baselines. **No observed primary MAE improvement; tumor-type-only is strongest.** Not clinical PK or verified biological mechanism.
- Additional cross-study observed GI50 study: **47 tumor-lineage-concordant NCI-60 cell lines**, 17,987 RNA features and 280 metabolites, with fixed RNA-only, metabolite-only, joint and tumor-adjusted Ridge/XGBoost comparisons. **No measured benefit of adding metabolomics to RNA in primary MAE**; one reclassified cell line explicitly excluded. These assays are NOT same-time or same-specimen multiomics.
- Public Wegler 2021/ Yin-Battail 2025 source-format inventory: bounded retrieval of PMC supporting metadata, workbook sheet names, and Zenodo record information. **No fitted CLint, donor-level link or reproduced SWEET clinical result is claimed**.
- **Thesis-locked Omics→PK evidence gate:** local, strictly fail-closed subject/time/drug/dose/concentration pairing for **RNA alone, proteomics alone, metabolomics alone, or 2–3 modalities together**; explicit cross-assay window. Strict pre-dose covariate timing prevents post-dose biological measurement leakage. Only aggregate QC is returned — **no PopPK fit or proof of clinical improvement**.
- **Wegler source-schema audit:** reads open 2021 supporting PDF and compartment-specific proteomics workbook temporarily in GitHub Actions with official PMC AWS URLs and fixed SHA-256, emits only structure/keyword counts, no donor data or invented CLint rows.
- Unit tests using synthetic data; no private patient files or cloud PK inference in GitHub Actions.

**Not implemented yet:** a *scientifically validated* inductive SWEET method (the implemented experimental method is explicitly *SWEET-inspired*), drug-specific mechanistic patient scores, Pillar-A overlays, XGBoost/PopPK fits, and clinical validation. The pinned upstream regression check reproduces raw example networks but flags an official z-score archive mismatch (issue #4). Treat the outputs as input preparation only.

## Run locally

Requires Python >=3.11.

    python -m pip install -e .
    omics2pk validate-exchange /path/to/MultiOmics-archive/exchange
    omics2pk prepare-sweet /path/to/MultiOmics-archive/exchange ./sweet-input
    # After independently running official SWEET locally:
    omics2pk summarize-sweet ./sweet-edge-outputs ./sweet-input/patient.txt ./sweet-features
    python -m unittest discover -s tests -v

Never commit the generated SWEET inputs for controlled or clinical subjects: they contain molecular measurements and subject-linked assay IDs.

## Protocol and milestones

- [SWEET network outputs and scientific caveats](docs/SWEET_OUTPUTS.md)
- [Experimental frozen-reference networks and leakage controls](docs/INDUCTIVE_SWEET_RESEARCH.md)
- [Native MultiOmics producer integration validation](docs/NATIVE_MULTIOMICS_INTEGRATION.md)
- [Public measured pharmacology evidence matrix](docs/EVIDENCE_MATRIX.md)
- [Doxorubicin observed-response negative benchmark results](docs/RESULTS_NCI60_DOXORUBICIN_20261016.md)
- [NCI-60 cross-study 47-line RNA + metabolite GI50 protocol (PR #9 merged)](docs/NCI60_TWOOMICS_BENCHMARK.md)
- [Observed cross-study NCI-60 two-omics negative benchmark results (47 cell lines)](docs/RESULTS_NCI60_TWOOMICS_20261016.md)
- [Wegler measured CLint source qualification and Yin SWEET source gates](docs/PUBLIC_PK_SOURCE_QUALIFICATION.md)
- [Thesis-aligned Omics→PK evidence ladder and future site plan](docs/OMICS_TO_PK_THESIS_PROOF_AND_SITE_PLAN.md)
- [Wegler 2021 PDF/proteome source-schema results and data availability blocker](docs/WEGLER_SOURCE_STRUCTURE_RESULTS_20261016.md)
- [**First real transcript→drug-metabolizing enzyme activity validation, 100 donor-matched livers (8 CYP outcomes)**](docs/EPA_100_DONOR_TRANSCRIPT_CYP_VALIDATION_20261016.md)
- [**First real CYP proteomics→microsomal tucatinib clearance validation, 21 donor-held-out folds**](docs/TUCATINIB_21_DONOR_MEASURED_VALIDATION_20261016.md)
- [Independent original-source audit and actual first validation gates (including mandatory RNA ABSTAIN)](docs/SOURCE_AUDIT_AND_FIRST_VALIDATION_GATES_20261016.md)
- [Qualified public Omics×PK data catalogue: in-vitro CLint, RNA functional assays, clinical TDM](docs/PUBLIC_OMICS_PK_DATASET_CATALOGUE_20261016.md)
- [Research roadmap and decision gates](docs/ROADMAP.md)
- [Comparative methods protocol: ML vs PopPK, SWEET and evidence graphs](docs/PROTOCOL.md)
- [Data security, local-only computation and publication](docs/DEPLOYMENT_AND_PUBLICATION.md)

## Publishing / licensing

The GitHub repository can be public without hosting any clinical records. New original source is offered under PolyForm Noncommercial 1.0.0 with the mandatory NOTICE attribution; authored texts are CC BY-NC-SA 4.0. See LICENSE, LICENSE-CONTENT.md and CITATION.cff. Commercial use requires separate permission from rights holders. This is source-available, **not OSI open-source**. We do not bundle private Thesis code or the upstream SWEET implementation.
