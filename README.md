# Omics2PK

**Research-only alpha. Not clinically validated.**

Omics2PK connects the reproducible outputs of [MultiOmics](https://github.com/rberrah/MultiOmics) to research on patient-specific molecular networks (SWEET), drug-specific knowledge from Thesis Pillar A and PopPK/PKPD covariate discovery in Pillar B (route R2).

## What is executable today

- Offline inspection of the **MultiOmics DOWNSTREAM_EXCHANGE_V1** directory: checks SHA-256, study-data access declaration, assay/sample/subject and feature-ID matching, finite numeric matrices and recorded processing scale.
- Conversion of suitable **normalized transcriptomic** matrices into the input files expected by the original SWEET implementation: expression.txt, patient.txt, gene.txt. The converter rejects unaddressed missing data, repeated subjects, incompatible namespaces and mixed tissues.
- Strict quality control of **already-computed SWEET** edge lists and extraction of descriptive network features for the proposed M1bis arm (network-wide edge statistics, per-gene signed/absolute strength). No mechanistic or CYP activity interpretation.
- Experimental training-frozen **SWEET-inspired** network proposal (distinct from upstream SWEET), with subject-level LOO/cross-fitted and held-out features, optionally separating genome-wide sample weights from targeted edge genes. Public 80-sample expression-only diagnostic, but no PK endpoint or validated causal interpretation.
- CI that executes the **real native MultiOmics exporter** on its bundled three-omic synthetic demo: 8 longitudinal subjects, 49 assays, raw RNA strictly blocked from SWEET.
- Research-only **real observed pharmacological benchmark** using 60 NCI-60 cell lines and predeclared doxorubicin GI50, with raw/PCA/network and fixed XGBoost (M1/M1bis/M2), negative shuffled-network control and tumor-type-adjusted baselines. **No observed primary MAE improvement; tumor-type-only is strongest.** Not clinical PK or verified biological mechanism.
- Additional cross-study observed GI50 study: **47 tumor-lineage-concordant NCI-60 cell lines**, 17,987 RNA features and 280 metabolites, with fixed RNA-only, metabolite-only, joint and tumor-adjusted Ridge/XGBoost comparisons. **No measured benefit of adding metabolomics to RNA in primary MAE**; one reclassified cell line explicitly excluded. These assays are NOT same-time or same-specimen multiomics.
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
- [NCI-60 cross-study 47-line RNA + metabolite GI50 protocol (exploratory, PR #9)](docs/NCI60_TWOOMICS_BENCHMARK.md)
- [Observed cross-study NCI-60 two-omics negative benchmark results (47 cell lines)](docs/RESULTS_NCI60_TWOOMICS_20261016.md)
- [Research roadmap and decision gates](docs/ROADMAP.md)
- [Comparative methods protocol: ML vs PopPK, SWEET and evidence graphs](docs/PROTOCOL.md)
- [Data security, local-only computation and publication](docs/DEPLOYMENT_AND_PUBLICATION.md)

## Publishing / licensing

The GitHub repository can be public without hosting any clinical records. New original source is offered under PolyForm Noncommercial 1.0.0 with the mandatory NOTICE attribution; authored texts are CC BY-NC-SA 4.0. See LICENSE, LICENSE-CONTENT.md and CITATION.cff. Commercial use requires separate permission from rights holders. This is source-available, **not OSI open-source**. We do not bundle private Thesis code or the upstream SWEET implementation.
