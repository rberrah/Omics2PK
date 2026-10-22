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
| **Zhang 2016 n=100** | authors report 10 CYP mRNA measurements | authors report matching CYP protein content | measured probe-substrate intrinsic activities | **P0 acquisition · `ABSTAIN_SOURCE_ROWS_NOT_ACQUIRED`** |
| **Ohtsuki 2012 n=17** | authors report matched mRNA | absolute CYP protein quantification | measured CYP functional activities | **P1 replication · `ABSTAIN_SOURCE_ROWS_NOT_ACQUIRED`** |
| U.S. EPA n≈100 | yes, CYP qPCR | no matched quantitative protein layer | measured CYP microsomal activity | `ABSTAIN_NO_MATCHED_PROTEIN` |
| Tucatinib/imatinib HLM n21/n17 | no matched RNA | yes, CYP proteins | measured CLint / substrate depletion and functional activity | `ABSTAIN_NO_MATCHED_RNA` |
| Wegler n15 liver donors | no qualified matched RNA in current source contract | yes, proteomics | measured in-vitro CLint across compounds | `ABSTAIN_NO_MATCHED_RNA` |
| Jamwal n74 human livers | study reports CYP3A4 mRNA | study reports CYP3A4 protein | measured midazolam 1'-hydroxylation / CLint biology | `ABSTAIN_SOURCE_ROWS_NOT_ACQUIRED` |

The highest-priority source is now **Zhang 2016** (DOI
https://doi.org/10.1124/jpet.116.233635): the paper reports 100 normal human
liver samples with ten CYP molecular/functional axes. The source article
establishes the experiment, but Omics2PK still lacks an authenticated original
donor-level table with stable donor IDs, mRNA, protein, probe activity and
units. Midazolam must be treated as a measured CYP3A probe activity rather than
silently relabelled as perfectly isoform-specific CYP3A4 activity.

**Ohtsuki 2012** (DOI https://doi.org/10.1124/dmd.111.042259) is the secondary
independent acquisition target: 17 human livers with absolute protein,
corresponding mRNA and functional activities. The paper explicitly places
individual-donor protein data in Supplemental Tables S5/S6, CYP assay conditions
in S1, and absolute enzyme activities in the supplemental data; analytical
under-LOQ protein values are also explicitly defined. A fresh source audit on
22 October 2026 probed both the exact historical ASPET DOI-supplement route and
the current DMD article route: **both returned HTTP 403**. Thus the experiment
is strongly publication-supported, but the original donor-level three-view
file is still not legitimately acquired.

Jamwal 2018 remains a CYP3A4-focused candidate. Its accessible material has not
provided an authenticated 74-row individual RNA+protein+activity matrix.
Figure digitization is not an acceptable substitute.

Detailed source acquisition contract:
`MULTICYP_SAME_DONOR_SOURCE_ACQUISITION_GATE_20261022.md`.

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

Two local numerical gates are available without bundling original donor data:

- `multicyp_rna_protein_activity_source_gate.py`: preregistered multi-CYP
  donor-LODO RNA / protein / RNA+protein comparison for Zhang n=100 or
  Ohtsuki n=17, requiring every published source donor to remain represented
  for every included CYP endpoint while BLQ/missing values remain explicit;
- `cyp3a4_rna_protein_mdz_source_gate.py`: Jamwal-specific n=74 comparator.

The multi-CYP gate refuses arbitrary donor subsets, mixed assay units, multiple
post-hoc probe activities for one CYP endpoint, pooled rows and unknown DOIs.

## Interpretation

At present Omics2PK has:
- genuine RNA → enzyme-activity evidence (EPA);
- genuine protein/activity → drug-specific CLint evidence (TKI/Wegler);
- several publications whose experiments contain the desired same-liver
  RNA+protein+activity bridge, with **Zhang n=100 now the primary acquisition
  target**, but no qualified original individual matrix in the workflow.

Therefore **RNA vs protein vs RNA+protein superiority is not yet numerically
established**. The software gate is ready; the remaining blocker is legitimate
original-source acquisition and donor-level provenance.
