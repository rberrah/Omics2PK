# Hepatic CYP metabolized-drug protein abundance + enzyme phenotype fusion, original-source report

**Date:** 18 October 2026. Original data, source identity checked in independent public studies with the **SAME human liver microsome donors**.

## Study aim and underlying data

Rather than using non-metabolized metformin to demonstrate hepatic metabolic clearance, we tested **tucatinib** and **imatinib**, both CYP2C8/CYP3A substrates. Original research articles:

- Tucatinib HLM, Latham et al., DOI https://doi.org/10.1016/j.dmd.2025.100061, official PMC https://pmc.ncbi.nlm.nih.gov/articles/PMC13095627/
- Imatinib HLM, Latham et al., DOI https://doi.org/10.1016/j.dmd.2025.100196, official PMC https://pmc.ncbi.nlm.nih.gov/articles/PMC12990092/

The **actual 22 printed Table-1 rows** (21 unique individual donor livers plus one excluded 150-donor pooled HLM) were automatically read *in memory only* from both official public PMC article pages, with source DOI checks and original donor-ID matching. All four published imatinib Spearman CYP associations independently reproduced (CYP2C8 abundance 0.909314, CYP3A abundance 0.901961, CYP2C8 probe activity 0.664216, CYP3A probe activity 0.862745). Same 21 donor IDs and CYP abundances/activities across the two articles. The four imatinib clearance nonquantifiable observations were **kept missing** (not zero). No original Table 1, donor identifier, individual protein measurements, clearances, nor model predictions are redistributed in Omics2PK source/CI artifacts.

Source-verified CI: https://github.com/rberrah/Omics2PK/actions/runs/39147377727 . The full-text Europe PMC XML was unavailable for these two manuscripts, so the downloader used each **official PMC HTML Table 1**, pinned as source SHA hashes in the run; no fabricated data, credentials or commercial source access.

## Predeclared model contrasts and same-donor methodology

Identical leave-one-original-biological-donor-out outer folds, Ridge α=1 fixed, train-only standardization of all features, physical negative predictions clipped to 0, no feature search or tuning:

1. no-omic benchmark — remaining donor mean clearance or metabolite peak area ratio;
2. **two CYP protein concentrations** (CYP2C8 and total CYP3A4+CYP3A5) from quantitative targeted proteomics;
3. **two functional enzyme probe activities** (CYP2C8 amodiaquine probe and CYP3A midazolam probe);
4. four predictors **CYP proteins + CYP probe activities**.

The 2-protein and 2-activity arms must reproduce existing source-gated LODO models to 0.000002 absolute MAE before 4-predictor fusion can be reported. Descriptive uncertainty is paired 1,999-replicate percentile bootstrap of already heldout donor-specific loss differences (no model refitting). The model hyperparameters are fixed, so this is internal CV, NOT true external validation.

### Results actually computed on original independently measured human HLM assays

MAE is given in each source-defined outcome unit. Smaller is better.

| Actual source in-vitro observed endpoint | Real independent donor n | Original donor-heldout mean MAE | CYP proteins MAE | CYP activities MAE | Combined protein+activity MAE | Combined - proteins |
|---|---:|---:|---:|---:|---:|---:|
| Tucatinib apparent HLM metabolic clearance, µL/min/mg | 21 | 6.245714 | 2.041447 | 1.940258 | **1.783073** | **−0.258374** |
| Imatinib apparent HLM metabolic clearance, µL/min/mg | 17 | 28.869894 | 14.166086 | 14.026988 | **13.620372** | −0.545714 |
| N-desmethyl imatinib formation peak area ratio (not CL), dimensionless | 21 | 0.146683 | 0.082222 | **0.058313** | 0.059077 | −0.023145 |

Descriptive paired 95% bootstrap intervals for the combined-minus-proteins MAE difference:
- Tucatinib: **[−0.891842, +0.243625]** µL/min/mg; combination better in 10/21 individuals and worse in 11/21 despite lower mean.
- Imatinib quantified clearance: **[−3.793775, +2.435248]** µL/min/mg; combination better in 8/17, worse in 9/17.
- N-desmethyl imatinib peak ratio: **[−0.047171, −0.001215]** PAR; combination better in 12/21, worse in 9/21. The functional activity-only arm has the numerically lowest MAE, so this is not an advantage of 4 features over 2 activities.

No multiplicity-controlled inference, no independent donor centre, no true population PK/plasma AUC measured here, and these two drugs are not biological external validation of each other because **the 21 donor livers are the same**.

## Crucial distinction: this is not transcriptomics + proteomics

CYP enzyme probe activity is a **quantitative functional phenotype**, not measured mRNA and not a distinct omic measurement. This experiment compares **one molecular omics (proteomics) alone, one functional phenotyping assay alone, and their fusion**. In these 21 donors there are no corresponding transcriptomic RNA and metabolomic matrices available in the official source Table 1. It would be scientifically false to conclude that multiomics beats (or fails to beat) transcriptomics alone for these drugs.

The existing actual matched **RNA + metabolome** comparison in Omics2PK on 47 NCI60 lines showed *no* extra benefit of metabolites compared with RNA in that **doxorubicin GI50 drug response**, but GI50 is **not liver metabolic clearance**. See [NCI-60 original data report](RESULTS_NCI60_TWOOMICS_20261016.md).

For a true RNA-vs-multiomics-for-CLint experiment, the prerequisite is **same donor, same drug and original measured PK-functional outcome** with contemporaneous RNA abundance, protein quantities and (optionally) metabolites. Required common cases and frozen train-only dimension reduction are nonnegotiable. Without those data, label the missing comparison `ABSTAIN_MISSING_SAME_DONOR_RNA`, not an efficacy result.

## Practical experimental priorities

1. Pursue **same-donor human liver RNA + CYP-protein + matched xenobiotic in-vitro CLint** with a *drug metabolized by hepatic enzymes*. Current qualified hepatic evidence supports CYP2C8/CYP3A imatinib/tucatinib (protein+activity), but no RNA. Other target options: CYP3A midazolam/tacrolimus with in vivo C(t)/AUC and original genotype and quantitated protein/RNA samples.
2. Once available, prespecify matched donor folds with clinical or mean control, RNA-only, proteins-only, RNA+proteins, metabolomics-only and all-measured layers; include a within-donor anti-leakage/randomized layer control and compare **paired** delta MAEs with cautious small-n intervals. Do not tune the omics layer count by looking at original outer test outcomes.
3. Keep metformin as a **renal OCT2/MATE transport and absorption example**, not hepatic CYP clearance. Do not assume ABCB1/P-gp controls its clinical intestinal absorption; PMAT/OCT transporters have much stronger source support.
4. Continue to keep the genuinely measured E3 metformin n32 study as a validated example of an original-source matched clinical data pipeline **with neutral results** when adding three metabolomics PCs, rather than trying to re-label renal CLr as hepatic CYP metabolism.

## Reproduction

- Engine: `src/omics2pk/hepatic_protein_activity_fusion.py`
- Ephemeral true original public-source table parser/runner: `scripts/benchmark_original_hepatic_tki_protein_activity.py`
- Tests: `tests/test_hepatic_protein_activity_fusion.py`
- Workflow: `.github/workflows/hepatic-cyp-original-protein-activity-fusion.yml`

Run: `python -m scripts.benchmark_original_hepatic_tki_protein_activity` with `beautifulsoup4==4.12.3` and internet access. The official paper rows are fetched in memory, not saved. All displayed numerical results are **aggregated CV scores** and source provenance checks only.
