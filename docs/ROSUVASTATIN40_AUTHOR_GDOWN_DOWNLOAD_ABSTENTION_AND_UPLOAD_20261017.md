# Rosuvastatin n40 — authenticated or owner-mediated original CSV transfer required

**Qualification performed 2026-10-17.** Primary manuscript Silveira et al., *Frontiers in Pharmacology*, DOI https://doi.org/10.3389/fphar.2021.752960, genuine 40-human predose LC-MS metabolomics and targeted rosuvastatin PK (AUC/Cmax). Original data accession MetaboLights **MTBLS3128**, original author repository https://github.com/GustavoHBDuarte/Metabolomics_PK_proj.

## Sources genuinely obtained

1. MetaboLights public original ISA sample table: **142 LC-MS source sample acquisitions**, two assay tables of **69 and 73**, strict sample-name matches (PR #49); raw `mzML` are deposited with the original accession. The original MAF metabolite assignments match those source sample names but have **zero numeric sample abundance cells**, so they are *not* ready-to-fit matrices (PR #51).
2. Original author 21 MB `Metabolomics_Elastic_net.ipynb` at immutable commit `dbf9182b0f72f8b290f623f7cc6560c5a5f99d73`. Author README explicitly mentions processed CSV feature matrices used for study Elastic Net; **five distinct original Google Drive links** embedded in five `pd.read_csv` code cells, not in the GitHub tree (PR #52).
3. Original five URL accesses by Python urllib returned only Google Drive HTML pages (PR #52).
4. **Second route actually executed with pinned gdown 5.2.0**, public original notebook IDs (extracted only in process memory), each file attempted independently. [Actual original-source action run #39026878636](https://github.com/rberrah/Omics2PK/actions/runs/39026878636): **5/5 `PUBLIC_GDRIVE_DOWNLOAD_UNAVAILABLE`**, no CSV recovered. Every original volunteer name, opaque Drive token/URL, code cell output and numeric PK target is intentionally omitted from GitHub Actions logs. The outcome was not invented, it was evaluated in the original source runner.
5. Public MetaboLights edit URL `https://www.ebi.ac.uk/metabolights/editor/MTBLS3128/files` was supplied by the user but is an interactive editor; official public `/ws/studies/MTBLS3128/files?include_raw_data=false` independently supplies ISA metadata. An interactive /editor/ page being visible to an authenticated user does not imply an HTTP API token for the ChatGPT assistant.

## Status, not a model score

**`AUTHOR_FIVE_PROCESSED_CSV_LINKS_CONFIRMED; ALL_FIVE_AUTOMATED_DOWNLOADS_BLOCKED; ORIGINAL_NUMERIC_PK_X_OMICS_SAME_PERSON_CROSSWALK_UNVERIFIED; CLINICAL_E3_ABSTAIN.`**

Author's reported R² (AUC ~0.93 and Cmax ~0.94 ESI+) under LOOCV are **reported only**, not replicated externally; original best-feature/hyperparameter search could bias estimates in n40 unless all tuning/selection and transformations occur inside folds.

## Minimal requested original file upload

The user has offered to download original data via a browser if automated link transfer fails. **Request all five original CSV files** referenced in the pinned notebook. At minimum, upload **ESI positive predose processed matrix** and **the true person-level quantitative rosuvastatin AUC and Cmax labels**, ideally also ESI negative and 44 clinical laboratory attributes. Keep any author-provided original subject IDs *untouched* (do not reorder/relabel datasets or sort and join by row index); originals should be uploaded to the private chat, never committed to this public repository. If authors only share CSV omics without person-level AUC/Cmax, state explicitly that E3 cannot yet be reproduced.

After successful upload: load originals privately, hashes/allowed-use recorded, classify measured PK units, study human biological ID group per original demographics, sample predose timing and modality, baseline-only/covariates, then prespecified nested held-out-by-human models/negative controls with no target leakage.

Script: `scripts/audit_rosuvastatin_author_gdown_csv.py`; dedicated source CI and offline synthetic tests detect access limitations without publishing patient data.
