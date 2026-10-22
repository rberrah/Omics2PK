# Hepatic same-donor RNA/protein/multi-view readiness — 22 Oct 2026

## Scientific question

Can measured hepatic RNA predict drug-specific metabolism sufficiently, or does
quantified protein add held-donor-out information on the **same original human
livers**?

This question cannot be answered by combining different studies merely because
they measure complementary layers. The donor identity intersection is part of
the scientific evidence.

## Current qualified sources

| Source | RNA | Protein | Drug-specific observed outcome | Same-donor fit status |
|---|---|---|---|---|
| U.S. EPA n≈100 | yes, CYP qPCR | no matched quantitative protein layer | measured CYP microsomal activity | `ABSTAIN_NO_MATCHED_PROTEIN` |
| Tucatinib/imatinib HLM n21/n17 | no matched RNA | yes, CYP proteins | measured CLint / substrate depletion and functional activity | `ABSTAIN_NO_MATCHED_RNA` |
| Wegler n15 liver donors | no qualified matched RNA in current source contract | yes, proteomics | measured in-vitro CLint across compounds | `ABSTAIN_NO_MATCHED_RNA` |
| Jamwal n74 human livers | study reports CYP3A4 mRNA | study reports CYP3A4 protein | measured midazolam 1'-hydroxylation / CLint biology | `ABSTAIN_SOURCE_ROWS_NOT_ACQUIRED` |

The Jamwal 2018 article and its freely described supporting information confirm
that RNA, CYP3A4 protein and midazolam metabolism were measured in the same
74-donor study. However, the currently accessible supporting material is
organized as histology/gender summaries and correlation figures rather than an
authenticated 74-row individual RNA+protein+activity matrix. Figure
digitization is not an acceptable substitute.

Source:
- DOI: https://doi.org/10.1021/acs.molpharmaceut.8b00159
- PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC7065849/

## Software gate

`src/omics2pk/hepatic_multiview_readiness.py` now blocks model fitting before
any numerical comparison unless:

1. the source is human and authenticated;
2. individual source rows have actually been acquired;
3. RNA, protein and drug-specific outcome donor identifiers are explicit;
4. the complete-case intersection is the same predeclared donor cohort across
   all compared arms.

Explicit statuses include:

- `ABSTAIN_SOURCE_ROWS_NOT_ACQUIRED`
- `ABSTAIN_NO_MATCHED_RNA`
- `ABSTAIN_NO_MATCHED_PROTEIN`
- `ABSTAIN_NO_MATCHED_DRUG_OUTCOME`
- `ABSTAIN_NO_SAME_DONOR_INTERSECTION`
- `ABSTAIN_INCOMPLETE_CASE_MISMATCH`
- `READY_SAME_DONOR_RNA_PROTEIN_DRUG_OUTCOME`

The existing Jamwal local benchmark
`src/omics2pk/cyp3a4_rna_protein_mdz_source_gate.py` remains the numerical
ablation once an authorized original 74-donor matrix is available.

## Interpretation

At present Omics2PK has:
- genuine RNA → enzyme-activity evidence (EPA);
- genuine protein/activity → drug-specific CLint evidence (TKI/Wegler);
- a study that biologically contains the desired RNA+protein+midazolam bridge
  (Jamwal), but not a qualified individual matrix in the repository workflow.

Therefore **RNA vs protein vs RNA+protein superiority is not yet numerically
established**. This is an acquisition limitation, not permission to fuse
unmatched studies.
