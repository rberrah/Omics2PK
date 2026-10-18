# MTBLS3128 vs original author GitHub — checked 18 October 2026

## Exact same source study, distinct data products

Silveira et al., Frontiers in Pharmacology 2021, DOI https://doi.org/10.3389/fphar.2021.752960, original predose serum LC-MS and independently measured rosuvastatin AUC/Cmax among 40 humans.

The original author code is https://github.com/GustavoHBDuarte/Metabolomics_PK_proj at immutable commit dbf9182b0f72f8b290f623f7cc6560c5a5f99d73. The entire git tree contains exactly three files: 21,266,591-byte Jupyter notebook, README, requirements; **zero CSV files**. The notebook refers to five Google Drive processed CSVs via five genuine pandas read_csv calls. Google Drive HTTPS/gdown original attempts did not retrieve any CSV rows. This is a machine-download failure, not proof of inaccessible or absent scientific data.

The author README explicitly describes public original raw spectral files at EBI https://www.ebi.ac.uk/metabolights/editor/MTBLS3128/overview. The original experimental MS samples underwent ProgenesisQI peak-picking/alignment, internal standard normalization, pooled-QC normalization, cube-root transformation and QC RSD selection before generating five derived modeling CSVs. **Raw mzML is not an already normalized 40-person biomarker matrix.**

## Sources newly qualified

- Public EBI original six ISA/MAF metadata files, 145 mzML + 156 ZIP in the FTP folder (verified official source census). The visible FTP names include no direct processed CSV/XLSX, and all 156 ZIP central directories were subsequently checked via bounded EBI HTTPS Range requests in the independent original-source audit (below).
- Publisher supplement https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8633954/supplementaryFiles: three original Excel workbooks, four worksheets. Previous inspection examined only the first-row headers. Original 18 Oct source-scoped CI https://github.com/rberrah/Omics2PK/actions/runs/39138162728 now checks every cell for AUC, Cmax, subject labels and numerical data alignment: **zero candidate 40-row original quantitative patient AUC+Cmax tables detected** using this prespecified recognition pattern. Some summary labels contain AUC or Cmax, but do not supply a validated 40-subject target matrix. Unconventional layouts cannot be logically excluded by header-pattern screening alone.
- Added bounded ZIP central-directory metadata probe with synthetic-only tests and filenames redacted. An initial FTP `MLSD` operation abstained and was fixed using source-verified `NLST` fallback. **The final original-data CI inspected all 156/156 ZIP central directories successfully, 3,174 internal source members and zero failures**: 1,526 `.dat`, 448 `.idx`, 604 `.inf`, 447 `.sts`, 149 `.txt`; **zero CSV, TSV, XLS/XLSX, Parquet, RDS or HDF5-format members**. Inspect [successful exhaustive original EBI ZIP metadata source run](https://github.com/rberrah/Omics2PK/actions/runs/39139008662). No instrument payloads downloaded. This is an extension/name audit, not proof that arbitrary text contents cannot contain numerical data.

## Evidence verdict

| Source question | Result |
|---|---|
| EBI contains genuine study original MS observations | YES, source data deposited |
| Author GitHub directly includes five processed modeling CSVs | NO, 3 committed non-CSV files only |
| Author code has five external CSV references | YES, pinned Google Drive links |
| Numeric rows and feature matrices from these five CSVs recovered | NOT YET |
| True 40-person AUC/Cmax table in publisher XLSX | NOT IDENTIFIED after all-cell structural scan |
| ZIP central-directory formats for all 156 official archives inspected | YES, 156/156 directories; 3,174 members; 0 CSV/XLSX/other named tabular datasets (no raw contents read) |
| Unique biological specimen-to-patient AUC/Cmax crosswalk | NOT QUALIFIED |
| Clinical E3 omics→PK gain independently reproduced | NO |

**Formal status:** RAW_EBI_SAME_STUDY_CONFIRMED; FIVE_PROCESSED_CSV_EXTERNAL_ONLY; SUPPLEMENT_NO_IDENTIFIED_40_SUBJECT_PK_TARGET_TABLE; ALL_156_ZIP_MEMBER_FORMATS_VERIFIED_NO_NAMED_STRUCTURED_CSV; CLINICAL_E3_ABSTAIN.

## Next gate

Priority is original five processed author CSVs from the URLs inside the notebook, downloaded through an authorized browser or author-provided file access, with original IDs retained; use locally only. If they remain unavailable, request original person-linked NCA AUC, Cmax, sample ID crosswalk and predose Progenesis intensity matrices from the authors. Re-processing all original mzML is a legitimate methodological fallback only if the actual PK targets and biological matching also exist. Then use the private source gate and 40-subject LODO scripts already in Omics2PK; external/patient-held-out claims require further independent qualification and nested feature tuning. Do not transfer individual PK/omics data to public CI.
