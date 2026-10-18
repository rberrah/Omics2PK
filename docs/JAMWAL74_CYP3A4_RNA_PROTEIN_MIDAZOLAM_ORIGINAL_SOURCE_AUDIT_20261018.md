# CYP3A4 hepatic metabolism: original n74 RNA vs protein + midazolam enzymatic endpoint
**Source audit date:** 2026-10-18. **Evidence classification: measured experiment published, true individual paired numerical dataset NOT independently acquired; quantitative fusion trial ABSTAIN.**

## Study and scientific relevance

Jamwal R et al., 2018, *Molecular Pharmaceutics*, DOI [10.1021/acs.molpharmaceut.8b00159](https://doi.org/10.1021/acs.molpharmaceut.8b00159), PMCID [PMC7065849](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065849/), PMID 29792708.

The authors analyzed CYP3A4 *mRNA*, CYP3A4 *protein* by targeted proteomics and *in-vitro* **midazolam 1′-hydroxylation/CYP3A4 activity** in **74 actual human liver donor samples** (normal n24, NAFL n26, NASH n24), plus PBPK simulations of theoretical midazolam exposure. In this study the ORIGINAL experimental endpoint is *in-vitro* enzyme activity, and PBPK predictions are NOT observed individual clinical midazolam AUC.

This is the desired biochemical bridge between measured RNA and quantified hepatic drug metabolism. However, **targeted CYP3A4 qPCR is not full transcriptomics** and one targeted CYP3A4 proteomic measurement is not the full proteome; if this study becomes usable it supports a conservative *two-assay modality ablation*, not unrestricted omics universality.

## Original supplemental source verification

Official manuscript supplementary URL advertised by the original PMC paper:
`NIHMS1566110-supplement-Supplementary_material.docx`. The publicly described supplement consists of:
- Table I: histopathology group summary (24/26/24 donor classes);
- Table II: original PCR primers;
- Table III: sex-stratified CYP activity/mRNA/protein summary;
- Tables IV/V: **PBPK PREDICTED** mean pharmacokinetic parameters (not observed individual clinical concentrations);
- Figures I–V: various summary trends and RNA–protein–activity correlations.

**Machine audit actually run on CI**, two official source methods:
1. [Europe PMC REST `PMC7065849/supplementaryFiles`](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7065849/supplementaryFiles) returned XML, **296 bytes**, `Content-Type: application/xml`, not a DOCX nor ZIP of measured participant records.
2. [Direct PMC original manuscript DOCX](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065849/bin/NIHMS1566110-supplement-Supplementary_material.docx) returned an HTTP error from the GitHub Actions runner.

These results prove **source retrieval failure in these endpoints**, not the absence of individual RNA/protein/activity measurements in all possible investigator repositories. The original publisher supplementary material is described publicly as primarily group summaries and correlation FIGURES. Thus we cannot assume the 74 genuine source tuples `(donor key, mRNA, protein, activity)` are available, cannot claim a numeric regression, and must NOT infer raw values by digitizing correlations/figures or reusing published R².

Original source audit: `scripts/audit_jamwal74_original_cyp3a4_rna_protein_activity.py`; tests: `tests/test_jamwal74_original_cyp3a4_source.py`; reproducible opt-in workflow: `.github/workflows/jamwal74-cyp3a4-original-source.yml`. The parser examines DOCX OOXML tables, embedded spreadsheets and direct XLSX original ZIP contents if the official source becomes accessible; **it never emits participant identifiers or raw clinical/molecular observations**. Its output explicitly reports `OFFICIAL_PMC7065849_SUPPLEMENT_UNAVAILABLE_ABSTAIN` and keeps numeric evidence flags **false**. This source access failure is not treated as a validation of absence or as a statistical null result.

## New original-data-only RNA/protein comparator — coded, unexecuted on Jamwal donors

Source-locked implementation: `src/omics2pk/cyp3a4_rna_protein_mdz_source_gate.py`.
Local runner: `scripts/benchmark_local_jamwal74_cyp3a4_rna_protein.py`.
Synthetic-only scientific/anti-leakage tests: `tests/test_cyp3a4_rna_protein_mdz_source_gate.py`.

A real authorized per-donor table must include **exactly 74 independent donor IDs, explicit source DOI and original row type, biological donor age/sex, original normal/NAFL/NASH histology group, actual positive numeric original CYP3A4 mRNA normalized to 18S, CYP3A4 protein concentration (pmol/mg), measured in-vitro midazolam 1′-hydroxylation (pmol/min/mg)**. This strict original study field schema is a private curation contract; it is not a claim that any released supplement already has those exact header names. Correct before original fit if article reports a different physical activity unit.

With an authentic and authorized original source, prespecified fixed leave-one-donor-out Ridge `α=1` and TRAIN-ONLY scaling compare **training-only mean**, **age+sex+NAFL/NASH**, **RNA only**, **protein only**, **RNA+protein**, **clinical+RNA**, **clinical+protein**, and **clinical+RNA+protein** on **identical original heldout people**. Primary contrast is *clinical+RNA+protein minus clinical+RNA MAE*; secondary is *RNA+protein minus RNA only*. The descriptive bootstrap resamples already heldout per-donor losses, **not new prospective donor/model-refit uncertainty**, so does not certify external performance. The numeric test will ABSTAIN on duplicates, missing values, group/source mismatches, pooled tissue, simulated PBPK, or wrong original donor number.

Run locally only after authorized original data are acquired:

```bash
python -m pip install -e .
python -m scripts.benchmark_local_jamwal74_cyp3a4_rna_protein /absolute/private/path/verified_74_donors.csv
```

The original donor CSV must be stored outside the public repo. Original identifiers, full data, fold predictions and rows are not printed, committed or uploaded. **The tests use fabricated values only; no source-original Jamwal model performance exists yet.**

## Current evidence and next highest-value source

- **EPA n100 human liver**: original qPCR CYP → observed CYP enzyme activities; genuine, independent original donor-heldout single RNA modality; **no same-person proteomics**.
- **Tucatinib/imatinib n21 HLM human donors**: original CYP proteins + measured functional phenotypes → observed metabolized drug clearance, genuine two-view source on same individuals; **no same-person RNA**.
- **Jamwal n74 human livers**: original study has both RNA and protein and observed substrate-specific midazolam activity in the experimental laboratory; **individual raw n74 paired numeric assay rows not obtained through official public source endpoints**, so no direct comparative claim.
- **NCI-60 47 matched lines**: original RNA and metabolite feature matrices with measured doxorubicin GI50; molecular fusion did not improve RNA; **GI50 is antitumor drug response, not hepatic CLint or systemic plasma PK**.
- **Rotroff n32 metformin clinical E3**: real predose metabolome and individual measured original AUC/renal CLr; no convincing incremental metabolomic gain; **metformin is not a hepatic CYP metabolic-clearance example**.

**Next scientific priority:** obtain the actual Jamwal 74-donor individual molecular+activity table or identify an alternative **genuinely same-donor human liver RNA+protein+metabolized-drug CLint** release. Manuscript summaries alone are insufficient; no extrapolation across unmatched EPA/TKI donors.

**Formal research status:** `ANALYSIS_ENGINE_UNIT_TESTED; ORIGINAL_JAMWAL74_INDIVIDUAL_DATA_NOT_OBTAINED; ORIGINAL_RNA_VS_PROTEIN_CYP_CLEARANCE_MODEL_ABSTAIN`.
