# Public source search: CYP3A5 genotype → omics → activity → drug metabolism (2026-11-03)

## Search outcome and scientific qualification

The full same-donor `CYP3A5*3/*6/*7 PGx + isoform-resolved RNA + specific CYP3A5 protein + independent pre-landmark metabolomics + CYP3A5-specific activity + clinical tacrolimus PK` **has NOT been found as an unrestricted, source-verified participant-level dataset**. Do not treat generic CYP3A probe activity as CYP3A5-specific or post-incubation drug metabolite formation as independent baseline metabolomics.

| Priority | Original source | Real same-donor measurements | Public access / required qualification |
|---|---|---|---|
| P0 immediate small **PHH imatinib** | Latham et al. 2025, DOI 10.1016/j.dmd.2025.100196, PMC12990092 **Table 4**, PRIDE PXD066482 | 14 PHH donor rows: CYP2C8/3A4/3A5 diplotypes, CYP2C8 and CYP3A4/5 probe activities, targeted CYP2C8/3A4/3A5 protein, imatinib in-vitro clearance and N-desmethyl formation. 13 donors have complete protein; one protein-missing row must remain missing | **14 directly visible public individual rows, no application needed for table**. Raw targeted proteomics deposited PRIDE. No RNA, no native untargeted predose metabolome, no genuine in-vivo plasma PK; CYP3A4/5 probe not CYP3A5-specific. n13 is *below* the standard benchmark's predeclared n30 threshold |
| P0 already research-owner acquired **HLC** | Sage syn4499, `syn88644` RNA, `syn89614` genotype, `syn4629` activity; Yang et al. 2010 DOI 10.1101/gr.103341.109; GEO GSE9588 | Prior original-source audit: 466 liver microarrays, ~387–397 source-matched probe activities; overlapping real genotypes n195 total, 171–172 complete for selected CYP endpoints | **Already qualified in Omics2PK via owner-supplied archives**, not newly obtained. CYP3A4 midazolam/testosterone activity; no separately proven CYP3A5-specific activity, no proteomics, no independently matched metabolome. Genotyping array provides rs776746 but does not document complete *3/*6/*7 coverage; do not relabel as a full diplotype |
| P1 new genotype–protein–activity small liver | Grangeon et al. 2025, DOI 10.3390/ph18121789, PMC12735975 | 16 genotyped liver microsome donors; CYP3A5 PGx in supplemental Table S1; CYP-specific protein and shared CYP3A4/5 activity provided as a study | Publication/supplement accessible; **individual quantitative protein+activity-by-PGx crosswalk not yet independently established**; must inspect all supplementary tables and source vendor assay-level records before any numerical regression |
| P1 broader genotype–protein | Pedersen et al. 2024, DOI 10.1124/dmd.124.001692 | 250 genotyped postmortem livers, 116 hepatic CYP protein assays including CYP3A5 | Paper confirmed; **downloadable unified n116 per-person genotype/protein table not yet confirmed**, enzyme activity/metabolome/PK not established |
| P1 preclinical CYP pathway | Latham et al. 2025, DOI 10.1016/j.dmd.2025.100061, PRIDE PXD057282 | 21 HLM CYP2C8/3A proteins+activities and real tucatinib metabolic depletion | Already analyzed via original article same-donor table in Omics2PK; no original public matching RNA or genotype reported in verified table |
| P1 source-available RNA–function | US EPA ORD-037038, DOI 10.23719/1519197 | Original public spreadsheet links real hepatic RNA transcripts and CYP probe activities | **Official public XLSX downloadable**; prior results already in repo. Does not establish CYP3A5-specific PGx/metabolomics/clinical PK |
| P1 true multi-omic fusion benchmark | EATRIS-Plus latest **v4**, Zenodo 19136013 (2026-03-20) | 125 healthy individuals, blood cell RNA, plasma proteome, metabolome, lipidome, multi-omic sample links | **v4 open**: R MultiAssayExperiment and Python MuData (141.7 MB) with metadata. Full germline WGS and expanded clinical phenotypes in ClinData are **restricted**. No drug-specific CYP3A activity or drug PK established; blood is not liver. Biological integration control, NOT CYP clearance superiority |
| P2 larger genotype/phenotype RNA | Wang/Collins et al. 2022, DOI 10.1111/cts.13398, PMC9652438 | 246 livers with variants including CYP3A5*3/*6/*7 and promoter variants, CYP3A4/5 RNA, relevant CYP3A4 protein associations | Article confirms same-cohort analyses; **full n246 individual matrix publicly downloadable not verified**. GEO GSE207590/207791 is ATAC/4C on small hepatocyte experiments, NOT a substitute for the original n246 donor matrix |
| P2 regulated molecular | EGA EGAS00001003426 | 150 matched liver samples, targeted 340 ADME genes DNA variants and 303 expression genes | **Access controlled**, no direct CYP function or original clinical drug PK proved from the indexed data |
| P2 regulated in-vivo biomarker | Izat et al. 2026, DOI 10.1002/cpt.70369, COCKTAIL | n54 4β-hydroxycholesterol/cholesterol, liver CYP3A4 protein, measured midazolam IV/oral PK (n30 longitudinal follow-up) | Strong potential real PK link but patient matrices require study governance approval. 4β-OHC is NOT CYP3A5-specific |

## Public PHH imatinib Table 4: preliminary source-reading computation

Using the **13 exact complete protein rows** in Latham Table 4 (not mixing with the 21 HLM donors), a private one-off script transcribed the published values with CYP3A5 diplotype classes (two *1 / one *1 / zero *1) and CYP2C8 carrier flag. No individual data or donor identifiers are committed. For each arm, a **same-donor LODO** standardize-on-train → fixed Ridge(alpha=1) on log(original source imatinib-clearance value) → exp transform. No hyperparameter search.

| Predictors (all n13, same LODO folds) | Descriptive heldout clearance MAE on printed Table-4 numeric scale |
|---|---:|
| CYP2C8 carrier status + CYP3A5*1 copy count | 2.90658 |
| Same PGx + CYP2C8/3A4/3A5 protein concentration | 3.00370 |
| Same PGx + CYP2C8 amodiaquine and CYP3A midazolam probe activities | **2.76294** |
| Same PGx + all above proteins and activities | 3.26281 |

This is **NOT a prospectively independent evaluation**: n13, class balance n9 no-CYP3A5*1/n3 one copy/n1 two copies; imatinib disposition in PHH depends strongly on CYP2C8; all units and original source crosswalk assumptions must be verified. The paper inconsistently describes clearance units as per million cells in Table 4 and per mg protein in part of the text, so preserve the printed numeric scale and **do not assert an absolute physical unit** until resolved. Absence of a gain for proteomics here is not evidence that proteins have no biological utility.

Do **not** use the measured N-desmethyl-imatinib formation product as a baseline metabolomic predictor of imatinib clearance: this is a contemporaneous drug-incubation outcome and would confound or leak function.

Original source: https://pmc.ncbi.nlm.nih.gov/articles/PMC12990092/ (Table 4); PRIDE: https://www.omicsdi.org/dataset/pride/PXD066482 .

## Next predeclared steps

1. **Reproduce the PHH n13 exercise** by a source-pinned table HTML parser with tests for expected 14 original IDs, 13 complete proteins, units, *3/*6/*7, no row-order mismatch. Keep it separate from the main n≥30 functional benchmark, with `PILOT_ONLY_NOT_GENERALIZABLE`.
2. **Use already-held HLC originals** for exact genotype+RNA vs activity assays, but explicitly call outcomes CYP3A4/CYP3A-shared and genotype a partial rs776746 call, not CYP3A5-specific (do not apply the complete-diplotype model).
3. **Inspect Grangeon v2025 n16** supplements for an independent individual phenotype, and full Genotype/protein/enzymatic 16× table. No fabricated values or mixed donors.
4. **Test EATRIS v4** with verified subject-level multiomic overlap and a held-donor identity/assay-provenance gate. This validates fusion, **not pharmacologic efficacy**.
5. For definitive comparative CYP3A5 and observed drug PK, seek original specific enzyme-activity data, genuinely baseline omics, and full CYP3A5*3/*6/*7 panel and CYP3A4*22 genotypes in the **same individuals**. Until then keep the central E2/E3 superiority endpoint `NOT_TESTABLE_FULL_COHORT`.

All 2026-11-03 public-source facts are research leads, not blanket license to redistribute third-party raw data.
