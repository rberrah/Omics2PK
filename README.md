# Omics2PK

**Research-only alpha. Not clinically validated.**

Omics2PK connects the reproducible outputs of [MultiOmics](https://github.com/rberrah/MultiOmics) to research on patient-specific molecular networks (SWEET), drug-specific knowledge from Thesis Pillar A and PopPK/PKPD covariate discovery in Pillar B (route R2).

## What is executable today

- Offline inspection of the **MultiOmics DOWNSTREAM_EXCHANGE_V1** directory: checks SHA-256, study-data access declaration, assay/sample/subject and feature-ID matching, finite numeric matrices and recorded processing scale.
- Conversion of suitable **normalized transcriptomic** matrices into the input files expected by the original SWEET implementation: expression.txt, patient.txt, gene.txt. The converter rejects unaddressed missing data, repeated subjects, incompatible namespaces and mixed tissues.
- Unit tests using entirely synthetic fixture files. No remote uploads or patient files in GitHub Actions.

**Not implemented yet:** SWEET network inference or exact published-output replication, patient network scores, Pillar-A overlays, XGBoost/PopPK fits, and clinical validation. Treat the outputs as input preparation only.

## Run locally

Requires Python >=3.11.

    python -m pip install -e .
    omics2pk validate-exchange /path/to/MultiOmics-archive/exchange
    omics2pk prepare-sweet /path/to/MultiOmics-archive/exchange ./sweet-input
    python -m unittest discover -s tests -v

Never commit the generated SWEET inputs for controlled or clinical subjects: they contain molecular measurements and subject-linked assay IDs.

## Protocol and milestones

- [Research roadmap and decision gates](docs/ROADMAP.md)
- [Comparative methods protocol: ML vs PopPK, SWEET and evidence graphs](docs/PROTOCOL.md)
- [Data security, publishing and optional Oracle compute](docs/DEPLOYMENT_AND_PUBLICATION.md)

## Publishing / licensing

The GitHub repository can be public without hosting any clinical records. Software licensing and ownership review remain pending; visibility alone does not grant redistribution rights. We do not bundle the private Thesis repository or upstream SWEET source.
