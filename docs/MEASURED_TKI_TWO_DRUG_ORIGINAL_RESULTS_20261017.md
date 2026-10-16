# Real original Table-1 CYP → in-vitro clearance: tucatinib and imatinib, 16–17 October 2026

**Scope: exploratory within-source reanalysis, n=21 ORIGINAL biological HLM donors shared across BOTH articles. Not clinical PK, not independent external validation, not validated IVIVE.**

Sources:
- Tucatinib, Latham et al., *Drug Metab Dispos* 2025, DOI [10.1016/j.dmd.2025.100061](https://doi.org/10.1016/j.dmd.2025.100061), [original PMC Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC13095627/), PRIDE PXD057282.
- Imatinib, Latham et al., *Drug Metab Dispos* 2025, DOI [10.1016/j.dmd.2025.100196](https://doi.org/10.1016/j.dmd.2025.100196), [original PMC Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC12990092/), PRIDE PXD066482.

## Original-source identity validation

**The actual 21 individual HLM donor IDs and source lot IDs are identical in both published Tables 1.** The original CYP2C8 abundance, CYP3A4+CYP3A5 abundance, amodiaquine activity and midazolam activity coincide by donor. The additional 150-donor pooled preparation is excluded in both. Four imatinib source donors (H0509, H0502, H0028, H0525) show *no measurable imatinib substrate depletion*, and original Table 1 has no imatinib clearance entry; they are missing/non-quantifiable outcomes, **not observed clearance=0**. Original 17 observed HLM imatinib CL: 2.5168 to 112.58 μL/min/mg. N-desmethyl-imatinib *peak area ratio at 10 min* is available in 21 biological donors. It is not a measured CL, and must remain a distinct target.

Original Table 1 positive CYP phenotypes and target columns were transcribed temporarily for private numerical reanalysis; no raw study rows, reconstructions of publication tables, identifiers, publisher archives or protected proteomic files were committed. Independent reanalysis reproduced both articles' **eight** published CYP↔clearance Spearman coefficients at the precision printed:

| Predictor | Tucatinib n=21 | Imatinib quantified n=17 |
|---|---:|---:|
| CYP2C8 abundance | 0.904190 (article 0.904) | 0.909314 (article 0.909) |
| CYP3A abundance | 0.881455 (article 0.882) | 0.901961 (article 0.902) |
| CYP2C8 probe activity | 0.781423 (article 0.781) | 0.664216 (article 0.664) |
| CYP3A probe activity | 0.907438 (article 0.907) | 0.862745 (article 0.863) |

## Actual measured-source held-donor-out exploratory comparisons

For every independently held-out donor, keep their entire profile out of the training set; mean reference, protein ridge and activity ridge are independently refit on the remaining donors. All ridge models have fixed α=1, use only two corresponding measured phenotypes, standardize on training donors only, and truncate nonphysical negative predictions at zero. Source values retain original published rounding. All scores are mean absolute error (MAE), in **endpoint-specific units**.

| Endpoint | n original biological donors | No-omics training mean MAE | CYP protein two-covariate MAE | CYP activity two-covariate MAE |
|---|---:|---:|---:|---:|
| Tucatinib apparent microsomal CL, μL/min/mg | 21 | 6.245714 | **2.041447** | 1.940258 |
| Tucatinib CL in **exact imatinib-observed donor subset** | 17 | 6.6684 | **2.1903** | 2.2867 |
| Imatinib apparent microsomal CL, μL/min/mg, quantifiable only | 17 | 28.869894 | **14.166086** | 14.026988 |
| N-desmethyl-imatinib formation peak-area ratio (PAR), dimensionless | 21 | 0.146683 | **0.082222** | 0.058313 |

Original 21-donor tucatinib estimates reproduce the previously recorded Omics2PK local evaluation (mean 6.246, proteins 2.041, activity 1.940). Relative MAE changes due to the two measured CYP protein covariates: **−67.3%** for tucatinib (21 donors), **−50.9%** for imatinib measured cases (17 donors), and **−43.9%** for metabolite formation PAR (21 donors). These percentages are improvements *within the same small study panel*, not external gains or true patient pharmacokinetic exposure predictions.

### Exploratory uncertainty and eligibility limitations

Illustrative 10,000-resample percentile bootstrap intervals of paired **already held-out per-donor absolute error differences** (protein minus no-omic mean, random generator seed 20261017):

- Tucatinib n21: ΔMAE −4.2043 μL/min/mg, approximate percentile interval **[−6.1711, −2.4684]**; 17/21 donors improved.
- Tucatinib restricted to same imatinib CL n17 donor IDs: ΔMAE −4.4781, interval **[−6.5498, −2.6469]**; 17/17 improved.
- Imatinib original measured n17: ΔMAE −14.7038 μL/min/mg, interval **[−24.5500, −4.4802]**; 15/17 improved.
- N-desmethyl-imatinib n21: ΔMAE −0.0645 PAR, interval **[−0.1024, −0.0244]**; 18/21 improved.

**Critical:** these are paired bootstraps of already generated LODO error values, *not* nested bootstrap refitting or confidence bounds on truly new donor populations; endpoints are correlated and multiplicity not controlled. Imatinib CL analysis excludes four donors with poor measurable depletion, which probably favors high-metabolism donors: the outcome-availability mechanism is not missing completely at random. A clinical or independent-donor study could reverse these findings. A censored original imatinib regression requires an assay-specific depletion detection threshold and source time data; neither should be invented from absent entries. The same shared 21 donors do not constitute independent cross-study replication, despite two independent molecules.

## Reproduction, source rights and next gates

The benchmark implementation is [src/omics2pk/matched_tki_metabolism.py](../src/omics2pk/matched_tki_metabolism.py), local CLI [scripts/validate_matched_tki_local.py](../scripts/validate_matched_tki_local.py), full source contract [here](TWO_DRUG_SHARED_DONOR_PUBLIC_SOURCES_20261016.md). To replicate the original values independently, curate the genuinely published Table-1 data yourself into local CSVs with explicit DOI, original donor IDs and NA for four imatinib clearances; execute the CLI and inspect published Spearman checks before reading scores. No original full table is bundled in GitHub due to the publisher's explicit text/data-mining and reuse reservation; only aggregate outcome analysis and tested reusable code are public.

**Next:** source-gated independent n26 NAT1/NAT2 sulfamethoxazole hepatocyte metabolism cohort; compare NAT2 genotype, two NAT enzymatic probes, NAT2 peptide abundance and clinical demographics **on arm-matched complete cases**. NAT/SMX endpoints are also *in vitro* metabolic formation, NOT observed in-vivo drug clearance, AUC, trough nor DILI.

**Clinical Omics→PK E3/E4 proof remains pending truly same-human linked dose/time/C(t)+pre-dose omics.** The same-donor mechanistic evidence is a preliminary biological bridge, not validation of omics-based individualized dosing.
