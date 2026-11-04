# Omics2PK — Omic context and mutation-dependent information gain

**Research protocol, 2026-11-04. Not a clinical PK decision rule.**
Related research branch: PR #137 (draft). Original human-source material must never be
committed or stored as workflow artifacts. No patient-level predictions output.

## Central biological question

The most informative omics view need not be the same for every genetic background.
A loss-of-function variant may prevent transcript abundance from tracking intact
protein or true enzymatic activity. Conversely, a valid transcript can be informative
where measured protein is noisy, or a metabolite can track a realized phenotype more
closely than expression. But metabolite abundance alone is often **not enzyme- or
transporter-specific**, and protein abundance alone is not functional activity.

Primary comparison:
1. matched baseline characteristics, + **authentic observed mutation information**;
2. + original **RNA** only / **protein** only / **metabolites** only;
3. all three combinations of two, + all three together;
4. prespecified **somatic TP53 status × molecular-view interactions**;
5. nested 80/90%-retained-gain **minimum-view selection separately by TP53
   stratum**, with `k=0`/no additional omic as a valid outcome.

The baseline and every view share the **same exact biological units, holdout
folds, outcome, and preprocessing safeguards**. Outer test targets do not
select a model or a molecular panel; if a genotype stratum is too small,
selection falls back to the global inner-training result and the fallback
count is reported explicitly. Relative performance must be compared on
original heldout MAE and uncertainty, **not a percentage of SHAP mass**.

## Dataset 1 — actual genomic mutations and three omics, but PD only

Candidate authentic source: [Bioconductor rcellminerData 2.34.0](https://bioconductor.org/packages/release/data/experiment/html/rcellminerData.html)
and [NCI CellMiner mutation documentation](https://discover.nci.nih.gov/cellminer/).
The `mut` layer reports **binary gene-level predicted-deleterious mutations**
(presumed predominantly somatic, not validated complete germline diplotypes).
They include event types such as frameshift, nonsense, splicing or selected
missense. A binary call is **not** allele dose, direct enzymatic measurement,
CYP3A5 star-allele genotype, or proof of complete protein inactivation.

Join to the already qualified **47** NCI-60 biological cell-line identities:
- original doxorubicin GI50 (cell-line **PD in vitro**), tumor lineage, RNA
  (17,987), original candidate metabolite signals (280);
- original Guo et al. SWATH protein intensities (3,171), independently
  publisher-SHA-pinned;
- original CellMiner mutation matrix, exact line name **and tissue-prefix**
  crosswalk, preserving missing gene calls as unavailable, never as WT;
- matching is by identical cell-line identity, NOT a same-day matched
  biological aliquot across independent laboratories.

Code: `scripts/export_nci60_cellminer_mutations.R`,
`src/omics2pk/mutation_multiview.py`,
`scripts/benchmark_original_nci60_mutation_multiview.py`,
`tests/test_mutation_multiview.py`.
Original-source workflow `.github/workflows/nci60-somatic-mutations-threeomics.yml`.
Mutation flags include prespecified TP53 and those direct doxorubicin-gene
mutation calls genuinely present, plus source-qualified mutation burden.
TP53 grouping is *a priori* (DNA damage response hypothesis) and is NOT a
pharmacogenetic CYP star-allele interpretation. All interactions are
regularized Ridge alpha10; original 3×5 repeated disjoint outer folds and
3-fold inner training selection. The three repetitions are **47**, not 141,
independent biological samples.

## Dataset 1 — second test: mutation-aware biological priors vs negative controls

Use exact source feature identities (not post hoc ranked picks):

| View | Fixed mechanism-inspired panel | Biological limit |
|---|---|---|
| **RNA** | TOP2A, TOP2B, ABCB1, ABCC1, ABCG2, CBR1, CBR3 | Gene expression may not represent functional protein, particularly with variants |
| **SWATH proteins** | Actually measured TOP2B, ABCB1, CBR1 | Missing TOP2A and other direct-prior proteins must not be imputed as present; peptides do not prove activity |
| **Metabolites** | Six original glutathione/redox source-name candidates | Hypothesis only: source chemical IDs/MSI level-1 identity and GSSG/NADPH/flux absent |

Each panel is evaluated with exactly the same **observed source somatic
mutation baseline** (tumor lineage, observed genotype flags, assay-derived
mutation burden) and optional TP53×panel interactions. The null control
draws **32 random feature panels of the same size**, source-matched on
*training-only* measurement detection rate and variance, independently
within each assay. This is a matched feature selection control, **not**
a network degree-conditioned full-omics STRING test or causal biological
function assay. Compare full-sample and TP53-specific MAE descriptively.

Code: `src/omics2pk/mutation_mechanism_controls.py`,
`scripts/benchmark_original_nci60_mutation_mechanism.py`,
`tests/test_mutation_mechanism_controls.py`.

## Dataset 2 — true pharmacogenetics, fewer available omics

Original Latham PHH Table 4, DOI `10.1016/j.dmd.2025.100196`: 14 human
hepatocyte donors, **13** complete common cases with CYP2C8 diplotypes and
CYP3A5 *1/*3/*6/*7, measured CYP proteins, independent CYP2C8 and
**shared CYP3A4/5** probe activities, and original *in-vitro* imatinib
clearance. Four existing original donor-LOO arms: PGx; PGx+protein;
PGx+probe activities; PGx+both. New descriptive breakdown by CYP3A5
nonexpressers versus other genotypes and CYP2C8 reference versus variant
diplotype. Groups with fewer than **five** donors **abstain** from
stratum-specific performance estimates; allele *6/*7 individual effects
are not estimated. No separate RNA or endogenous untargeted metabolomics
exists in this original table; adding a surrogate would fabricate a view.
This is not tacrolimus plasma PK and the shared midazolam/CYP3A assay is
**not** isoform-selective CYP3A5 function.

Existing reproducible human-source pilot:
`src/omics2pk/latham_phh_source_gate.py`,
`scripts/benchmark_latham_phh_original.py`,
`tests/test_latham_phh_source_gate.py`.

## Mandatory acceptance criteria

- Source authenticity and same-biological-identity alignment (fail closed),
  gene/variant assay provenance, mutation carrier counts, no mutation
  missing-as-WT, no cross-fold preprocessing leakage.
- Every fixed arm and nested selector compared on the same outer observations;
  no best modality chosen after viewing test performance.
- Subgroup gains require enough independent carriers/noncarriers;
  estimate interaction uncertainty without treating 3 repeated predictions
  as 3 patients.
- Mechanistic links cannot be inferred from SHAP/STRING alone.
- For the thesis primary goal, a **separate true human cohort with complete
  germline CYP3A5/CYP3A4 PGx + omics + measured enzyme activity/PK** is
  still mandatory. BIOMARGIN/ABIS tacrolimus linking is the priority source.

**Status:** workflows are source-qualification gates. An offline Python
test passing is *software qualification only*, not evidence of a real
mutation-modality superiority. Record original numerical results here
only after the actual CellMiner + 47-line + all assay joins succeed.
