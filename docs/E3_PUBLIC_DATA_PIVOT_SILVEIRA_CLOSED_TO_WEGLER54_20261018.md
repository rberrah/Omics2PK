# E3 public-source pivot after five private author Google Drive CSVs (18 October 2026)

## Stop condition: original Silveira/Duarte n40 rosuvastatin predose LC-MS

Original primary study: DOI 10.3389/fphar.2021.752960, MetaboLights MTBLS3128. The exact author GitHub tree (notebook, README, requirements) contains no processed CSV. The original notebook points to five Google Drive CSVs. **The user has directly confirmed the Google Drive is private/inaccessible. Stop attempts to download or treat these CSV links as public.** Independent official EBI inventory confirms 145 raw mzML and 156 instrument ZIPs with 3,174 internal instrument-format files; there are zero named CSV/XLSX/tabular-format files inside the ZIPs. The publisher supplement did not provide a recognized original 40-person numeric AUC+Cmax table.

The raw LC-MS are genuine and available for a future independent processing analysis, but this cannot independently reconstruct individual clinical PK outcomes or establish biological subject IDs. Thus Silveira is **NOT a reproduced Omics→PK E3 validation**. Status: AUTHOR_PROCESSED_CSV_CLOSED; RAW_SPECTRA_AVAILABLE; CLINICAL_E3_ABSTAIN until direct authorization and complete data from original investigators. Never invent source PK values, infer clinical ID pairing from row order, or digitize mean curves into per-person values.

## New highest-priority candidate: Wegler et al. rosuvastatin n54

Wegler et al. Clinical Pharmacology & Therapeutics (2021), DOI 10.1002/cpt.2056, PMCID PMC7984432:
- **54 independent clinical patients**, oral 20 mg rosuvastatin; the authors measured actual plasma profiles (0–24 h) and individual NCA AUC and Cmax, plus matched liver biopsies and measured OATP1B1, OATP1B3, OATP2B1, NTCP transporter protein concentrations.
- The primary article explicitly documents that original measured clinical pharmacokinetics and liver proteomics were paired in the same people. The reported AUC twofold prediction success (78%) is the AUTHORS' original result, not an independent Omics2PK result.
- The published Supplementary Data S1 is an ~12 KB original XLSX described as concentrations of transporter proteins and donor allometrics; it is NOT automatically proof that individual measured PK AUC/Cmax columns are also released.
- A read-only schema auditor examines original Europe PMC publisher supplementary ZIP in memory without logging participant identifiers or original values. This source gate must identify actual numeric individual AUC/Cmax and matched same-person protein rows **before** running any genuine E3 modeling.
- This original n54 Wegler *in vivo rosuvastatin transport and proteomics cohort* is DIFFERENT from the n15 Wegler in-vitro CLint and protein dataset already discussed in Omics2PK issue #11.

## Fallback protocol if 54-subject numeric PK source cannot be recovered

1. Mark n54 E3 as SOURCE_BLOCKED, not as a negative result.
2. Continue existing genuine E2 donor-matched sources: 100 independent human livers paired qPCR CYP transcript levels and microsomal enzyme activities; 21 individual human liver microsomes with tucatinib measured in-vitro CLint and proteomics; other source-qualified donor data.
3. For actual clinical E3, prefer fully public measured original patient omics and PK numeric target data; metformin n17 has 255 actual PK concentrations but source predose molecular matrices need proof; Phapale n29 tacrolimus also needs original participant-specific mapping.
4. For future n54 E3, differentiate hepatic protein measurement obtained during biopsy after drug dosing from a pre-dose noninvasive prediction biomarker: real transporter-based bottom-up PK scaling can be mechanistically valid, but it is not prospective pre-dose blood metabolomics.
5. Keep subject-level observed data, source identifiers and profiles out of public GitHub and CI; source audits output aggregates only. No external R², no invented confidence intervals or simulated biological ground truth.

## Primary sources

- Silveira original: https://doi.org/10.3389/fphar.2021.752960
- MetaboLights: https://www.ebi.ac.uk/metabolights/editor/MTBLS3128/overview
- Author GitHub: https://github.com/GustavoHBDuarte/Metabolomics_PK_proj
- Wegler n54 original: https://doi.org/10.1002/cpt.2056
- Wegler original published supplementary: https://ascpt.onlinelibrary.wiley.com/doi/10.1002/cpt.2056
- Europe PMC original article: https://europepmc.org/articles/PMC7984432
