# NCI-60 cross-study RNA + metabolite GI50 comparison (exploratory)

## Question
For a frozen set of NCI-60 cancer cell-line identities with **observed**
doxorubicin response, does transcriptomics plus metabolomics predict
in-vitro GI50 more accurately than transcriptomics alone or a tumor-type
baseline? This is **not** a human PK study.

## Data sources (pinned by full Git SHA)
- mixOmics-org/mixOmics 9f9fd66b69105fdb618b54588b638b567e85c0d6:
  60 NCI-60 cell lines, measured response to doxorubicin.
- Mathelab/NCI60_GeneMetabolite_Data
  4c140708c18c41910458414a626b1a51d3dd2a85:
  IntLIM public normalized expression (Affymetrix microarray; 17,987
  genes across 57 lines) and metabolite abundance (Metabolon; 280
  metabolites across 58 lines).

The IntLIM vignette explicitly describes upstream log2 transformation,
gene probe annotation/filtering, and technical-metabolite repeat
preprocessing. They cannot be treated as newly measured paired assays.
The two repositories represent **distinct original investigations and
laboratory conditions**. Cell-line identity is not the same as
same-sample, same-time, or experimentally paired multiomics. No
causal, CYP activity, mechanistic target, or PK claim is possible here.

## Prespecified eligibility
- Exactly **48** initially aligned lines must have collision-free case/punctuation keys
  across NCI-60 measured-drug response, IntLIM metadata, RNA and
  metabolite matrices.
- Independent source tumor-lineage labels must agree after a
  transparent small abbreviation harmonization (e.g., LEUK/LEUKEMIA,
  MELAN/MELANOMA, OVAR/OVARIAN). **MDA-MB-435 must be excluded**:
  mixOmics labels it BREAST, whereas IntLIM labels it Melanoma.
  The prespecified modeling cohort therefore has **47** lines.
- 12 syntactically unresolved names and the one histology-discordant
  MDA-MB-435 stay excluded; no fuzzy or historical aliasing.
- Input row orders are checked **exactly** and never inferred from
  numeric matrix position. Temporary crosswalk keeps source labels.
- All source-derived molecular matrices are temporary on the runner;
  only aggregate errors and feature counts are downloadable.

## Models and experimental integrity
- Outcome: measured doxorubicin GI50 (one drug declared before
  evaluation). Exactly the same 47 cell lines are used in every arm.
- Five outer subject-disjoint folds, three prespecified fold seeds
  20261016, 20261017, 20261018: repeated observations of **the same 47
  lines**, not 144 independent lines.
- Baselines: training-fold mean; cancer lineage alone.
- Fold-only preprocessing: >=50% nonmissing training values per feature,
  training-only median imputation, top 300 RNA / 120 metabolite
  TRAINING-variance features, STANDARD SCALE/PCA-3 separately for each
  view. Upstream full-study normalization cannot be undone and
  remains a specific limitation.
- Fixed Ridge alpha 10; comparisons RNA, metabolites, both, each with
  and without cancer lineage; raw selected-expression comparators.
  Secondary XGBoost on 3-PC views, fixed parameters. A shuffled
  TRAIN-only metabolite representation is a sanity negative control.
- Metrics: same-fold out-of-fold MAE and RMSE; differences versus
  lineage-only comparator. A descriptive paired bootstrap resamples
  47 independent line identities, averaging repeated-fold error
  within line. No inferential superiority claim or external validation.

## Interpretation / publication guard
Positive differences in one response/drug/dataset must not be presented
as a clinical PK result, as proof of multimodal mechanism, or as
independent confirmation of SWEET. The source transcriptome is not
genotype, drug transporter activity, or a new sample-specific gene
network. SWEET/Pillar-A analysis remains a **separate, not-yet-validated**
phase. Even a positive multiomics result here requires truly paired
specimens and external replication before mechanism or MIPD relevance.

## Execution

    Rscript scripts/export_public_nci60_twoomics.R <mixOmics-dir> <IntLIM-dir> ./temp_source
    python scripts/benchmark_public_nci60_twoomics.py --data ./temp_source --out ./aggregate.json --repeats 3
    rm -rf ./temp_source

GitHub Actions executes this with pinned public source commits, checks
the synthetic leakage tests, and retains **aggregate JSON only**.
Measured negative-increment results from actual successful CI are recorded
in docs/RESULTS_NCI60_TWOOMICS_20261016.md. This is an exploratory
no-go for claiming any incremental metabolomics benefit or clinical PK
relevance, not a no-go for the wider Omics2PK research program.
