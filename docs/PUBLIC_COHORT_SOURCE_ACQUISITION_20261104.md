# 2026-11-04: public E2 human-liver acquisition decision — Zhang, Liu, and NEW PXD085173

## Result: source discovery is NOT a donor-level data release

Original article and public repository records were cross-checked, independently of the already-frozen same-donor comparison code. No unverified original rows were downloaded, reconstructed from plots, or added to GitHub.

| Lead | Publication and authentic source | Matched measurements **asserted by authors** | Individual numeric crosswalk independently verified now? | Scientific use and current status |
|---|---|---|---|---|
| **P0 Zhang n100** | Zhang et al. 2016, [10.1124/jpet.116.233635](https://pubmed.ncbi.nlm.nih.gov/27189963/), PMID 27189963. Same-liver-bank corroboration [10.1124/dmd.116.069849](https://pubmed.ncbi.nlm.nih.gov/27271371/). | 100 normal human livers: RNA, absolute/targeted CYP proteins, measured CYP probe activities (10 isoforms), some PGx. Some measurements belong to different denominators (per mg microsomal protein vs per pmol CYP isoform). | **NO**. Primary publisher web abstracts and scholarly indexes give no source-authenticated 100×CYP numerical table in current access. | First choice for the actual RNA-versus-protein functional E2 benchmark; original source request required. |
| **P0 Liu n55** | Liu et al. 2016, [10.1038/srep26544](https://www.nature.com/articles/srep26544), [PMC4876377](https://pmc.ncbi.nlm.nih.gov/articles/PMC4876377/), PMID 27211076. | 55 human livers: CYP3A4*1G, CYP3A5*3, 13 miRNAs, CYP3A4/5 mRNA, shared CYP3A testosterone 6β-hydroxylation activity, **in-vitro** atorvastatin depletion plus ortho/para-hydroxy products. | **NO**. Public supplementary PDF provides aggregate association, regression and genotype-stratum summaries but no authenticated **55-row donor-matched cross-layer table**. | Uniquely attractive for genotype + regulatory miRNA + RNA → shared function → drug-specific in-vitro metabolism; no source-supported CYP3A5*6/*7 complete panel, no CYP3A5-selective catalytic assay and no human plasma atorvastatin PK. |
| **P1 NEW PXD085173 n61** | University of Washington Human Liver Bank S9 quantitative DIA, [ProteomeXchange/PanoramaPublic PXD085173](https://www.omicsdi.org/dataset/panorama/PXD085173), **announced 2026-10-23**; associated preprint Wen et al. [10.64898/2026.10.12.758325](https://doi.org/10.64898/2026.10.12.758325) **published 2026-10-19**. | **61 donor-derived human liver S9** fractions: quantitative DIA proteomics, processed peptide/protein abundances, with the **same 61-HLS9 panel** assayed for CYP2C9-mediated diclofenac and tolbutamide hydroxylation. Preprint claims functional contribution of FABP1 beyond CYP2C9 protein. | **PARTIAL**: the *proteomics dataset itself* is publicly indexed as available at multiple processing levels, but a 61×{FABP1,CYP2C9 protein, two observed metabolite-formation rates} **same-ID table and functional endpoint download have NOT been verified**. PanoramaWeb source deep link returned HTTP 403 during this audit. | Best new near-term mechanistic **protein+context→activity** replication; two-drug same-donor falsification. No proven paired mRNA, genotypes, metabolomics or **observed clinical PK**, so never market it as a full multomic or PGx superiority experiment. |
| **P1 Ohtsuki n17** | Ohtsuki et al. 2012, [10.1124/dmd.111.042259](https://pubmed.ncbi.nlm.nih.gov/21994437/), PMID 21994437 | n17 donor liver biopsies with CYP proteins, matched transcripts and functional probe activities. Original supplementary tables S5/S6 reportedly contain individual protein measurements; other layers must be joined via original specimen IDs. | **NO** source-authenticated all-layer crosswalk in current environment. | Independent, smaller replication if complete original supplementary material or permission to use a deidentified row-level extract can be obtained. |

## Biological interpretation: important positive example of a *context protein* (not an omics shortcut)

FABP1 is not a CYP gene or an enzyme-specific metabolic activity measure. Its quantitative abundance may modify substrate availability and CYP2C9 catalytic kinetics, potentially via FABP1–CYP2C9 protein interactions. The n61 source offers a prospectively specified test of **measured CYP2C9 protein alone versus measured CYP2C9+FABP1 proteins**, applied to **both diclofenac and tolbutamide** using exactly the same donors and fold partitions.

This is *one proteomic modality with two proteins*, NOT transcriptomics + proteomics or de novo metabolomic functional phenotype. Measured drug-specific hydroxylation after incubation is the **outcome**, not an independent baseline metabolite feature. The published models' performance must not be described as our independent validation.

CYP3A5 *3/*6/*7 knowledge rules must **not** be reused for this CYP2C9 experiment. If PGx becomes available, a separately authenticated CYP2C9 star allele/phenotype mapping would be required.

## Original source access gates — NO fabricated cohort values

### P0 Zhang

Request exact originally measured values and sample crosswalk, retaining **all 100 biological donors**, CYP/probe identity, and for each CYP donor RNA, protein (physical units), microsomal/isoform activity (physical units), genotype where originally measured, source analytical/BLQ statuses, experimental replicates and external preprocessing information. Keep original sample rows private; only aggregate results to repository.

**Verified publication contact:** Hai-Ling Qiao, `qiaohl@zzu.edu.cn` (explicitly listed by PubMed for PMID 27189963; historical address, delivery not verified).

### P0 Liu

Request exactly **55 individual source liver donors** and columns: donor pseudonym, genotype CYP3A5*3 and CYP3A4*1G with missing statuses, each of the 13 **named/strand-specific** miRNAs, CYP3A4/CYP3A5 qPCR values, independently measured shared testosterone-6β-CYP3A activity, measured atorvastatin parent depletion and ortho/para-hydroxylation, units/LOD and experimental source identifiers. Do **not** substitute plotted summary ranks/correlations or convert CYP3A activity to CYP3A5-only function.

**Verified original corresponding author:** Shi-Long Zhong, `zhongsl@hotmail.com` (explicitly published on page 1 of the open Scientific Reports article PDF; historical address, delivery not verified).

### P1 NEW FABP1/CYP2C9 n61

The ProteomeXchange record describes downloadable **DIA proteomic processing levels**, but a drug-metabolism phenotype file cannot be assumed from this description. First attempt official PanoramaPublic processed *protein* and sample manifest access. Only then request **the original 61 biological liver IDs linking processed CYP2C9 and FABP1 protein amounts** to *directly measured* diclofenac 4'-hydroxylation and tolbutamide 4-hydroxylation (specify source assays, initial substrate concentrations, measured units, BLQ values, replicate-aggregation and clinical liver disease covariates). Request whether any donor-level CYP2C9 genotype / RNA / baseline metabolomics was actually measured; absence is acceptable.

**Verified public dataset submitter:** Alex Zelter, `azelter@uw.edu` (OmicsDI PXD record). Original research team corresponding author Nina Isoherranen, `ni2@uw.edu` (University of Washington publications). Submission does **not** demonstrate permission to share the 61 drug metabolism observations.

## Pre-registered comparative assays after acquisition

### Zhang and Ohtsuki (RNA/protein quantitative function)

Use the pre-existing conservative `src/omics2pk/multicyp_rna_protein_activity_source_gate.py` and [source request manifest](ZHANG_OHTSUKI_AUTHOR_SOURCE_REQUEST_MANIFEST_20261022.md). Within each genuine source and each original enzyme/probe, benchmark training mean, RNA alone, protein alone, RNA + protein with same complete-case donors and no cross-source joins or rank reconstruction. Treat CYP3A midazolam as **shared CYP3A**.

### Liu (regulatory PGx/RNA/miRNA/drug outcome)

Only after n55 authentic joined rows: fixed held-donor comparison **CYP3A4/5 genotype only → +mRNA → +miRNA → +both**, prespecify each miRNA mature/star annotation and abstain on missing. Secondary endpoint: independent, observed atorvastatin parent depletion and OAT/PAT product formation *in-vitro*, with explanatory enzyme activity treated as a separate intervention/arm and **never used as a predictor of the same-assay reaction**. Unlike clinical E3, drug substrate measurements are contemporaneous assays, not plasma PK.

### PXD085173 (one-modality plus contextual mechanism, two substrates)

Original 61 verified donors only. At the **same original donor splits**, compare 4 arms:
1. other donors' mean;
2. CYP2C9 protein alone;
3. FABP1 protein alone (explicit negative/explanatory control);
4. **CYP2C9 + FABP1 protein**.

Do this **separately** for two observed drug hydroxylation rates, with secondary paired comparison of whether gain transfers across the two substrates. Predeclare Ridge alpha 1, train-only scaling, held-out donor MAE and paired loss contrasts. If a single donor lacks outcome/protein, preserve the full source manifest and use the **same declared complete-case subset for all arms of that drug**; do not silently drop donors from only the losing arm. Clinical pathologies may influence both protein and activity; if an original disease field exists use a prespecified baseline and disease-stratified diagnostics. **No benefit claimed until a 61-row original source crosswalk is validated.**

## Acquisition-ready request language (English, for the user to send, NOT sent)

**Zhang / Liu:** `Subject: Request for de-identified same-donor CYP multi-layer data — [DOI]`

Dear Dr [published corresponding author],

I am conducting academic research on whether molecular measurements (RNA, microRNA, protein) improve out-of-sample prediction of measured human hepatic enzyme activity beyond genotype alone. Your study [DOI] reports these measurements in the same liver donor cohort. Would it be possible to obtain the original supplementary source table or a de-identified donor-aligned extract containing the original assay values, sample keys, units, missingness/BLQ markers and probe-method descriptions? We do not require identifiable clinical information; the raw source would be held outside the public GitHub repository. We would cite your work appropriately and are open to any reasonable data-use requirements.

Sincerely, [researcher, institutional signature]

**Wen / Zelter / Isoherranen:** `Subject: PXD085173 — matching the 61 donor-level proteomics with CYP2C9 metabolism rates`

Dear Dr [recipient],

I am evaluating held-donor prediction of CYP2C9 functional activity from independent hepatic proteomic measurements. I located your new PXD085173 DIA dataset and the associated 2026 preprint (10.64898/2026.10.12.758325). Could you confirm whether the processed 61-donor protein abundance matrix includes an anonymized key that matches the individually measured diclofenac and tolbutamide hydroxylation rates? If the activity table is not in the public deposit, could you share an approved de-identified phenotype/key file, including assay units, BLQ status, original liver-disease status if available and sample-level measurement provenance? No clinical identifiers are requested, and the raw data would not be published on GitHub. We would acknowledge the source and follow the applicable permissions.

Sincerely, [researcher, institutional signature]

## Current acquisition outcome

- `ZHANG_100=ABSTAIN_AUTHENTIC_DONOR_ROWS_NOT_ACQUIRED`
- `LIU_55=ABSTAIN_AUTHENTIC_DONOR_ROWS_NOT_ACQUIRED`
- `OHTSUKI_17=ABSTAIN_AUTHENTIC_DONOR_ROWS_NOT_ACQUIRED`
- `PXD085173_61_PROTEOMICS=PUBLIC_DATASET_INDEX_CONFIRMED_BUT_LOCAL_ROWS_NOT_VALIDATED`
- `PXD085173_61_ACTIVITY_AND_PROTEOME_CROSSWALK=ABSTAIN_ORIGINAL_ACTIVITY_KEY_NOT_ACQUIRED`

Keep the previous `GSE49012` source-phenotype/publication parity blocker separate; **do not treat PXD085173 as a correction or replacement for that cohort**.
