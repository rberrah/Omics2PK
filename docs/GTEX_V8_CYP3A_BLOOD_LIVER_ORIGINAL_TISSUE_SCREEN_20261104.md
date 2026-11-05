# GTEx V8 — empirical tissue-source qualification for CYP3A4/CYP3A5
**2026-11-04, verified ORIGINAL DCC public source, SUCCESS:** [GTEx source-replayed CI](https://github.com/rberrah/Omics2PK/actions/runs/41899982546).

## Source and provenance
Original GTEx Adult V8 official DCC median TPM GCT (public, RNA-seq bulk tissue, not patients treated with tacrolimus): `https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz`.

Source-only Python script `scripts/audit_gtex_v8_cyp3a_blood_liver_medians.py` downloads **the entire original compressed published median GCT** and verifies the GCT structure, exact official tissue column labels (not API names), unique gene symbols, finite nonnegative median TPMs and original row count. The workflow stores only **aggregate source gene/tissue medians and file hash**.

**Do not relabel `Whole Blood` as PBMC, `Liver` as purified hepatocytes, or bulk values as enzyme activities.** The Small Intestine - Terminal Ileum segment is not necessarily proximal intestinal CYP3A first-pass metabolism; it is a separate anatomical region.

## Original results, real public GTEx V8 median TPMs

| Source-identified gene | Liver | Whole Blood | Small Intestine - Terminal Ileum |
|---|---:|---:|---:|
| CYP3A4 (`ENSG00000160868.14`) | **335.297** | **0.0501182** | 48.2465 |
| CYP3A5 (`ENSG00000106258.13`) | **155.853** | **0.0870020** | 61.4648 |
| FMO3 (hepatic context; `ENSG00000007933.12`) | 152.507 | 0.0374734 | 1.88056 |
| ALB (liver-context positive control; `ENSG00000163631.16`) | 25201.3 | 2.46845 | 1.21139 |
| PTPRC (hematopoietic marker; `ENSG00000081237.18`) | 3.85455 | **218.423** | 28.2409 |

**Tissue-median comparison only.** CYP3A4/CYP3A5 RNA are almost undetectable in median bulk whole blood despite prominent liver and ileal medians. The blood compartment captures immune-lineage material, a different biology. This supports tissue-specific qualification and challenges a simplistic « blood CYP3A4 RNA = liver CYP3A4 metabolic capacity » assumption.

**Non-inferable from the actual table:** whether *within-donor* whole-blood expression is correlated with liver expression, whether an individual blood inflammatory signature predicts liver CYP3A activity, isoform-specific CYP3A5*3 expression, protein-level active enzyme fraction, intestinal contribution to oral bioavailability or **clinical tacrolimus clearance/AUC**. The median-tissue GCT holds no individual sample IDs, genomic diplotypes or PK endpoints.

For blood-to-liver correlation, require original same-donor matched sample-level **GTEx gene TPM GCT** (much larger), de-identified original tissue sample attributes and donor-level crosswalk; quantify detection floor and measurement variance. Never infer correlation from two group medians. Full donor genotype data are controlled-access; avoid implying complete genotype is public.

## PK implications
1. Tissue-local CYP gene RNA is a **mechanistic gold-standard explanatory variable only for that tissue**, not a clinically accessible blood biomarker.
2. Compare accessible predose blood-derived metabolomic, immune and EV markers to genotype + real clinical PopPK, **not** to liver-expressed RNA model in independent unrelated people.
3. First existing authentic source: Rotroff 2016 metformin n33 (n32 measured AUC/CLr and predose GC-TOF), already executed in [original CI](https://github.com/rberrah/Omics2PK/actions/runs/39145479441): clinical reference versus + predose PCs showed **no convincing added benefit**. No source genotype or PopPK comparator has yet been fitted, so still an *initial E3* proof, not the full hypothesis.
4. Rosuvastatin MTBLS3128 n40 is very close to the desired human design, but **original paired quantitative AUC/Cmax and predose processed intensity CSV have not been acquired from authentic linked source rows** despite exhaustive source audit. No rosuvastatin Omics→PK gain can be asserted in Omics2PK.

References: [GTEx consortium 2020 Science](https://doi.org/10.1126/science.aaz1776), [original GTEx DCC open-access files](https://gtexportal.org/home/downloads/adult-gtex/), [Rotroff 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4906013/), [Silveira 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8633954/).
