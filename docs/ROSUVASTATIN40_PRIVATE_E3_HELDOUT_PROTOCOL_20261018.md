# E3 clinical rosuvastatin: private source-gated held-human-out reference (2026-10-18)

**Purpose:** prepare a *scientifically constrained* first independent human metabolome→quantitative PK benchmark **without claiming the original n40 matched matrix has been obtained**.

Primary original study: Silveira et al. 2021, DOI [10.3389/fphar.2021.752960](https://doi.org/10.3389/fphar.2021.752960), [MTBLS3128](https://www.ebi.ac.uk/metabolights/MTBLS3128). The original five processed author CSVs / individual subject-linked AUC and Cmax remain a source-retrieval blocker. The original MAF files do **not** contain numerical abundance values. See [local original CSV acquisition and structural gate](ROSUVASTATIN40_LOCAL_ORIGINAL_CSV_PRIVATE_UPLOAD_PROTOCOL_20261017.md).

## Workflow (only after obtaining the actual original author CSVs)

1. Run the earlier \`scripts/qualify_local_original_rosuvastatin_csvs.py\` inventory and strict join. Do not match by row number or assume 142 LC-MS acquisitions = 142 humans. Verify the source is original, **each of 40 biological human subjects** has a measured quantitative AUC/Cmax, selected metabolite intensities were measured **before drug exposure**, and original units and plasma methods are known.
2. Prepare a **private** author-verified *combined* CSV containing exactly one row per original volunteer. The merge must use original biological IDs or an independently authenticated one-to-one author crosswalk. Do not put the combined file in GitHub, public Actions, website assets or CI.
3. Independently record \`sha256sum /path/to/original_verified_combined.csv\` (or \`shasum -a 256\`) with the source documentation and units in the approved local environment. Freeze molecular feature names **before examining PK associations**, avoiding any target-based column screen.
4. Run the script locally with explicit column names and four manual attestations (do not assert these unless supported by inspected source files):

\`\`\`bash
python -m pip install -e .
python scripts/benchmark_private_original_rosuvastatin.py \
  --original-combined /PRIVATE/PATH/original_verified_combined.csv \
  --original-sha256 THE_ACTUAL_LOWERCASE_SHA256 \
  --person-id-column ORIGINAL_AUTHOR_SUBJECT_ID \
  --auc-column ORIGINAL_AUTHOR_AUC \
  --cmax-column ORIGINAL_AUTHOR_CMAX \
  --clinical-column AGE_NUMERIC \
  --clinical-column BMI_NUMERIC \
  --omic-column ORIGINAL_PREDOSE_INTENSITY_A \
  --omic-column ORIGINAL_PREDOSE_INTENSITY_B \
  --attest-original-author-csv-authenticated \
  --attest-same-biological-subject-ids-verified \
  --attest-all-selected-omic-signals-predose-verified \
  --attest-auc-cmax-assay-units-verified
\`\`\`

The uppercase values are placeholders, **not guessed names**. Do not give attestation flags unless a researcher has positively verified each condition. At least two named genuinely numeric predose molecular variables and an original 40-row file are required. Numeric clinical covariates are optional; when absent the comparator is a train-only constant, **not** a clinical+PGx reference.

## Fixed evaluation / protected statistics

- Outer **40 leave-one-biological-volunteer-out** folds. Original ID controls grouping. Sorting by patient ID (never emitted) makes the aggregate result invariant to CSV row order. Neither outcome is a predictor, and any extra postdose PK descriptor in source columns prevents running.
- **Predeclared ridge with no tuned hyperparameters**, on \`log(AUC)\` and \`log(Cmax)\` separately. The feature matrix is imputed by **training-subject** median (omics only), log1p-transformed (omics only), and z-scaled on the training fold. Feature dimensionality must be bounded. The fixed dual ridge penalty equals the number of included predictors; the clinical/constant baseline and clinical+omics arm are fit on exactly the **same original 40 held-out people**.
- Report original-unit held-out MAE, RMSE, bias and within-person paired delta MAE. Report a 2,000-replicate subject-resampled exploratory percentile 95% interval, **not** an external validation confidence claim. Each outcome tested separately; no multiplicity-adjusted significance or mechanism attribution is implied.
- **No train/test feature selection based on PK outcomes is allowed**. If the researcher chooses among the published 4,800 configurations after inspecting outcomes, this protocol is not confirmatory; source provenance, assay batch, missingness, analyte identity and calibration remain to be audited separately.
- Input subject identities, metabolite values, measured PK values, selected feature names, individual predictions and residuals are **never printed or uploaded**. Only numeric aggregate model scores, SHA and structural counters appear on standard output. Even aggregate output must be handled as study-controlled when needed.
- The four attestations indicate **manual researcher statements**, not independent checks conducted by software. They cannot prove clinical validity by themselves. The returned \`PRIVATE_SOURCE_ATTESTED_EXPLORATORY_E3_LODO_NOT_EXTERNAL_VALIDATION\` status is intentionally conservative.

## Validation and evidence limits

Synthetic-only regression tests in \`tests/test_rosuvastatin_heldout.py\` cover source SHA mutation, omitted attestations, row shuffle invariance, duplicate biological donors, target leakage, non-numeric signals, extra postdose outcomes, incomplete features and the clinical-vs-constant comparator distinction. \`python -m unittest discover -s tests -v\` runs the offline tests; no original CSV is bundled.

**Not demonstrated yet:** source-concordant numerical n40 results; independent replication of published AUC/Cmax R²; external generalizability; SWEET benefit; incremental value beyond genotype or clinical covariates; SAEM/NLME, a concentration-time PopPK fit, PBPK, or a dosing recommendation.

## Next mandatory scientific steps

After receiving the five authentic author CSVs, verify source hash and exact subject crosswalk, actual plasma PK units, predose timings and basic feature annotations; compare clinical-only/clinical+metabolomics. Then freeze an independent preregistered nested-XGBoost/MOFA2/SWEET plan with inner **training-only** tuning and pathway/feature selection, plus properly defined permutation nulls and optional external cohort. Omics2PK is developed separately from MultiOmics, PMx Explain and MIPD-Engine.
