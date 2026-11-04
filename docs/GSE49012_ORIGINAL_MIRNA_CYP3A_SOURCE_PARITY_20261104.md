# Original GEO GSE49012: liver miRNA versus measured CYP3A activity — source-parity ABSTENTION
**Date:** 2026-11-04. **Primary original source:** [NCBI GEO GSE49012](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49012). **Publication:** Vuppalanchi et al., *PLoS ONE* 2013, DOI [10.1371/journal.pone.0074471](https://pmc.ncbi.nlm.nih.gov/articles/PMC3772944/), PMID 24058572.

## Objective and what was actually downloaded

A genuine independently retrieved, original GEO family SOFT source was downloaded in-memory from:

`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE49nnn/GSE49012/soft/GSE49012_family.soft.gz`

SHA-256 of the original compressed SOFT bytes:
`d653ef3aef9fa9e237d0cd3ce5842e50c88ccdbf28418919e5f4deff9e5ed170`

The source-identity gate verified exactly **34 original biological liver samples** (22 status DISEASE/cirrhotic and 12 NORMAL), each with its own original *GSM* ID, original tissue/platform GPL17470, individually linked **normalized miRNA microarray intensities** and a positive *GEO-supplied* field `hepatic cyp3a activity`. The source contains 738 probe identifiers, of which **683** are human miRNA/let feature names; control probes are excluded before any biological PCA. Sample rows/labels, original per-person activity measurements and individual predictions are never committed or exported.

This is **shared CYP3A midazolam microsomal probe activity**, not CYP3A5-selective T-5 N-oxidation and not original midazolam plasma CL/F, tacrolimus AUC or PopPK. No actual CYP3A5 genotypes or per-person CYP proteins were identified in the SOFT.

Verified reproducible original-source workflow: [GitHub Actions #41818473268](https://github.com/rberrah/Omics2PK/actions/runs/41818473268). The paper's total functional cohort is **28 cirrhosis plus 12 normal**, while the miRNA GEO array is **22 cirrhosis plus 12 normal**. This difference was handled explicitly; 22 is not silently made into 28.

## Decisive publication-versus-GEO phenotype discrepancy

The original paper reports hepatic shared-CYP3A activity **0.28 ± 0.06 (mean ± SEM)** in **28 cirrhotic livers**, and **0.47 ± 0.07** in **12 normal livers**. The original linked per-sample GEO activity fields give these aggregate numbers:

| Original phenotype measure | GEO same-sample n | GEO sample mean | GEO sample sum | Paper summary |
|---|---:|---:|---:|---|
| Cirrhosis | 22 | **0.569682** | **12.533** | n28, mean **0.28 ± 0.06** |
| Normal | 12 | **0.470833** | 5.650 | n12, mean **0.47 ± 0.07** |

Under the plausible, explicitly unverified interpretation that the 22 GEO cirrhotic cases are members of the published 28 and that the paper summarizes the *same* nonnegative activity endpoint in the *same* units, the 28 patients' total activity is bounded by `28 × (0.28 + 0.005) = 7.98` allowing rounding to two decimals. **The observed 22-person subset already sums to 12.533 > 7.98.** The two source claims cannot simultaneously be true under those assumptions.

The normal group does reproduce the original paper mean. Thus, simply transforming all values by an arbitrary common unit factor is scientifically unsupported and would break the normal-group concordance.

**SOURCE QUALIFICATION STATUS: `ABSTAIN_GEO_ACTIVITY_PUBLICATION_GROUP_PARITY_FAILURE`.** This is a **source identity/phenotype-level discrepancy**, not model underperformance. It can reflect metadata errors in some GEO samples, mismatched source outcome definitions, mistaken sample subset assumptions or a publication summary problem; the present evidence does not uniquely identify the cause. No values are corrected, removed, inferred from figures, winsorized to force agreement, or overwritten. Original authors/source holders must resolve the sample-identity and phenotype provenance if this is to support a biological inference.

Programmatic gate: `src/omics2pk/geo_cyp3a_mirna.py` now publishes only group aggregates and a fail-closed `source_publication_parity_audit`, explicitly setting `biological_superiority_inference_permitted=False`. The original article publication means and group n's are declared constants for this source-specific audit, not claims of general biological CYP3A activity reference ranges.

## Descriptive *diagnostics only*: heldout regression results on deposited fields

These results were computed on the authentic GEO per-sample numerical fields, with exact same-sample 5fold disease-stratified heldout partitions, fixed `Ridge(alpha=1)`, **train-only** feature standardization. The original source already provides normalized intensities; original upstream study-wide normalization was not retrospectively re-fitted. For the PCA arm, median imputation, scaling and PCA are all trained within each outer training partition. Disease is **always** included for each miRNA comparator to prevent crediting disease-related signals to specific molecular features.

| Study-data diagnostic arm | GEO-supplied measured-field heldout MAE | Median absolute error |
|---|---:|---:|
| Disease status only | **0.36701085** | **0.15924637** |
| Disease + hsa-miR-142-5p | 0.37637518 | 0.18335158 |
| Disease + hsa-miR-142-3p | 0.37355470 | 0.17189699 |
| Disease + hsa-miR-130a* | 0.37589540 | 0.16419292 |
| Disease + hsa-miR-21* | 0.37723702 | 0.17605478 |
| Disease + hsa-miR-27b* | 0.37855820 | 0.18048298 |
| Disease + unsupervised 683-miRNA PCA2 | 0.38586772 | 0.21110742 |
| Disease + unsupervised 683-miRNA PCA4 | 0.38246172 | 0.17633043 |

**DO NOT interpret this as a negative biological validation of the omic hypothesis.** The unresolved source-phenotype discrepancy prevents inferring anything reliably about the activity actually analyzed in the paper.

Exploratory 95% disease-stratified **paired** percentile bootstrap of already-heldout absolute-error differences, 1,999 resamples (not refitted models):
- PCA2 minus disease ΔMAE = +0.01885687, interval **[−0.00605829, +0.04457683]**;
- PCA4 ΔMAE = +0.01545087, interval **[−0.02235770, +0.05131000]**;
- miR142-5p ΔMAE = +0.00936433, interval **[−0.01004107, +0.02902459]**.

All three intervals include zero; they are also exploratory *conditional on the deposited source values* and do not fix the source discrepancy. The result is not an independent external validation or a calibrated treatment decision.

## Crucial microRNA identity issue: mature versus star strand

Some source probes are `hsa-miR-27b*`, `hsa-miR-21*`, `hsa-miR-130a*`; source also reports both `hsa-miR-142-3p` and `hsa-miR-142-5p`. Do not silently equate a historical `*` strand with a miRNA named without the star, or substitute miRNA-3p for -5p. A modern mature-miRNA release/miRBase sequence-level mapping with probe sequences and documented accession is needed before merging with independent Liu n55 mRNA/miRNA findings. The original Liu 2016 supplementary file contains group/correlation tables but **not an authenticated full n55 individual response matrix**.

## Implementation and action needed

- Parser/benchmark: `src/omics2pk/geo_cyp3a_mirna.py`
- Public original SOFT replay: `scripts/benchmark_original_geo_gse49012_mirna.py`
- Synthetic-only contracts: `tests/test_geo_cyp3a_mirna.py`
- CI: `.github/workflows/gse49012-original-activity.yml`

Minimum original-source correction request to data authors:
1. Confirm whether 22 GEO cirrhotic GSM entries are exactly a subset of the n28 liver phenotype cohort used for Table 1.
2. Confirm per-GSM original assay units, substrate concentration, model and transformation of `hepatic cyp3a activity`, and whether paper Table 1 reports this same measure.
3. Explain the discordant large cirrhotic activities; provide provenance-preserving corrected GEO annotation or a verified donor-to-assay crosswalk, with author-approved release if appropriate.
4. Re-run the fixed same-donor models only after source original phenotype has passed this qualification; do not select favorable omics by examining the disputed outcomes.

**Final scientific judgment:** an original accessible miRNA-to-CYP3A assay dataset exists, but **as of this audit its activity labels cannot be reconciled with the published cirrhotic-cohort summary**. It is a valid software/source-QC stress test, **not an E2 biological efficacy proof**; it says nothing about genotype vs omics or individual tacrolimus PK.
