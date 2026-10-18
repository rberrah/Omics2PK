> **Qualification finale (18 octobre 2026) : VERROU LEVÉ.** Les sections ci-dessous conservent l'historique du *premier audit provisoire* qui ne reconnaissait pas les échantillons codés. Le contrôle approfondi a établi **33/33 appariements exacts** entre les identifiants cliniques numériques et les 33 triplets A/B/C métabolomiques, et 32 sujets avec AUClast0–24 et CLr réellement mesurées. Le premier benchmark E3 clinique a été exécuté dans [run 39145479441](https://github.com/rberrah/Omics2PK/actions/runs/39145479441) : aucun gain convaincant de l'ajout de trois composantes omiques, ni pour AUC ni pour CLr. Voir le [rapport E3 final](E3_REAL_METFORMIN33_PREDOSE_384_OMICS_AUC_CLR_HELDOUT_20261018.md). Le statut provisoire de non-appariement exposé plus bas est **obsolète** et conservé seulement pour la traçabilité de l'audit.

# Original Rotroff metformin n33 — data-level clinical Omics→PK E3 qualification (18 October 2026)

## Literature evidence

Original peer-reviewed publication: Rotroff DM et al., *Frontiers in Pharmacology* (2016), DOI https://doi.org/10.3389/fphar.2016.00135, original Europe PMC study PMCID **PMC4906013** and [publisher supplementary ZIP source](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4906013/supplementaryFiles).

Article reports **33 originally dosed people, 384 GC-TOF metabolomic features** and original individual metformin plasma concentration/exposure PK endpoints, with sampling at true predose timepoint A, and after dose at timepoints B and C. This makes it a genuine **E3 biological-study candidate**, unlike in-vitro E2 or simulated source data.

## Actual public-source inspections

Source verification performed against the **official real publisher supplement** in the opt-in GitHub Actions [original source run 39144061209](https://github.com/rberrah/Omics2PK/actions/runs/39144061209) and [independent cross-sheet follow-up run 39144421104](https://github.com/rberrah/Omics2PK/actions/runs/39144421104).

The original supplementary ZIP contained exactly **five .xlsx books, one .docx and ten images**. Examined all original Excel sheets and their actual numeric cell counts in memory, without uploading or emitting individual original rows.

- One original sheet: **392 rows**, up to **103 numeric measurements per row** (consistent with many metabolomic features arranged along rows, but sample/timepoint mapping must be demonstrated).
- Another original sheet: **118 rows**, up to **384 numeric measurements per row** (an original individual-sample candidate feature matrix, with biological/pre-dose status still unqualified).
- Other source tables: **34-row × up to 17-cell** and **34-row × up to 12-cell** numeric clinical/covariate tables, a 58-row sheet, and annotated statistical tables. Presence of 34 rows does NOT identify 33 biological donors without labels; do not infer person pairing from row count.
- In one 34-row source table, original cells contain **AUC tokens** (twice); in another, a **metformin token**, but no independently verified original header directly representing individual patient **metformin AUC AND Cmax** in one source sheet. AUC may refer to *glucose-response AUC* rather than *metformin exposure AUC*, requiring original manuscript/codebook attribution. This distinction is fundamental.
- An expanded crosssheet diagnostic tested exact matching of nonnumeric pseudonym candidates appearing in the first three columns of each original possible clinical 34-row table against original wide metabolomic matrices. **ZERO exact casefolded key overlaps across the four candidate sheet pairs**. No individual measured human AUC/Cmax ↔ original predose metabolite row mapping was verified.
- The source supplementary ZIP outer SHA varied across successive downloads (possible dynamic archive assembly). **Pin the SHA of each original actual .xlsx member separately, not the dynamically generated outer transport ZIP**, and do not infer scientific change solely from outer digest.

### Decision / qualification gate

**Source status:** ORIGINAL_MTBLS3128_AUTHOR_CSV_CLOSED; NEW_PUBLIC_ROTROFF33_ORIGINAL_OMICS_FEATURES_VISIBLE; PUBLIC_NUMERIC_CLINICAL_TABLES_VISIBLE; NO_ORIGINAL_UNAMBIGUOUS_SAME_PERSON_PREDOSE_OMICS_METFORMIN_AUC_CMAX_JOIN; **CLINICAL_E3_ABSTAIN**.

This is stronger than a generic literature match but **NOT** a valid original 33-person individual Omics→PK benchmark. It would be scientifically incorrect to sort 33 metabolomic profiles and 33 clinical PK rows into the same order, to equate drug AUC to glucose AUC, or to train/test with postdose metabolites as allegedly predose predictive factors.

Source scripts:
- scripts/audit_rotroff33_metformin_original_pk_omics.py
- scripts/audit_rotroff33_crosssheet_pk_omics.py
- scripts/audit_rotroff33_original_sample_naming.py
- tests/test_rotroff33_metformin_original_pk_omics.py
- tests/test_rotroff33_crosssheet_pk_omics.py
- tests/test_rotroff33_original_sample_naming.py

## Next defensible actions

1. Verify the original data dictionary for the omic specimen source names, timepoint A/B/C, independent patient identifiers and the clinical table's metabolic-vs-drug AUC wording; count actual 33-person intersections using only explicitly documented source keys.
2. If this exact source joint map is available, stage the legitimate patient rows **locally only** (do not publish original person-level omics/PK through GitHub); test source hashes, measured units and input causal ordering, then fit training-only nested feature selection and donor-heldout benchmark with a clinical baseline and negative controls.
3. If no original numeric drug PK targets and explicit crosswalk can be established, record E3 missingness and **continue genuinely paired E2 EPA n100 human donor transcript→functional CYP measurements** and tucatinib n21 HLM proteomics→CLint; never re-label as patient PK.
4. Request donor-to-sample crosswalk and independent drug Cmax/AUC from original authors only with research permission, without assuming those fields exist in supplement. Avoid closed Google Drive of another rosuvastatin study.

No original human subject values, sample labels or original proprietary data are part of this public repository. This audit is reproducible and fail-closed.
