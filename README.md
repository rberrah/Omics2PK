# Omics2PK

**Research-only alpha. Not clinically validated.**

Omics2PK connects the reproducible outputs of [MultiOmics](https://github.com/rberrah/MultiOmics) to research on patient-specific molecular networks (SWEET), drug-specific knowledge from Thesis Pillar A and PopPK/PKPD covariate discovery in Pillar B (route R2).

## What is executable today

- Offline inspection of the **MultiOmics DOWNSTREAM_EXCHANGE_V1** directory: checks SHA-256, study-data access declaration, assay/sample/subject and feature-ID matching, finite numeric matrices and recorded processing scale.
- Conversion of suitable **normalized transcriptomic** matrices into the input files expected by the original SWEET implementation: expression.txt, patient.txt, gene.txt. The converter rejects unaddressed missing data, repeated subjects, incompatible namespaces and mixed tissues.
- Strict quality control of **already-computed SWEET** edge lists and extraction of descriptive network features for the proposed M1bis arm (network-wide edge statistics, per-gene signed/absolute strength). No mechanistic or CYP activity interpretation.
- Unit tests using entirely synthetic fixture files. No remote uploads or patient files in GitHub Actions.

**Not implemented yet:** SWEET network inference or exact published-output replication, drug-specific mechanistic patient scores, Pillar-A overlays, XGBoost/PopPK fits, and clinical validation. Treat the outputs as input preparation only.

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
- [Research roadmap and decision gates](docs/ROADMAP.md)
- [Comparative methods protocol: ML vs PopPK, SWEET and evidence graphs](docs/PROTOCOL.md)
- [Data security, local-only computation and publication](docs/DEPLOYMENT_AND_PUBLICATION.md)

## Publishing / licensing

The GitHub repository can be public without hosting any clinical records. New original source is offered under PolyForm Noncommercial 1.0.0 with the mandatory NOTICE attribution; authored texts are CC BY-NC-SA 4.0. See LICENSE, LICENSE-CONTENT.md and CITATION.cff. Commercial use requires separate permission from rights holders. This is source-available, **not OSI open-source**. We do not bundle private Thesis code or the upstream SWEET implementation.
