# Jarmusch13 real private MS1 → MultiOmics input → Omics2PK PK parity (19 October 2026)

## What has actually been checked

The owner personally downloaded original public GNPS FBMN job `b3966445abad4658a3cdc63c8198a666`; the original ZIP (SHA256 `2047ea2032d73e3b48e2696fb1f652a7dedaa6f8cd5d40d0926d810bdd46c026`) remains **PRIVATE** and is not part of this repository. The independent, pseudonymized, non-published input package `jarmusch13_private_MultiOmics_Omics2PK.zip` has SHA256 `95a556d7ca288da9cf41f45d58212e05bd8e4e531b05e830505ede37cef8f294`. It contains 13 genuine human **predose, day1 T0 plasma MS1** assays, 6349 original GNPS-processed MS1 ion features, a wide and long matrix, assay metadata, tentative spectral annotations, and a separate **private clinical/PGx sidecar**.

New importer `scripts/validate_private_jarmusch13_multiomics_pipeline.py` fails closed on the exact private source bytes, declared GNPS job, expected CSV member names, one biological human per sample, common tissue, day1 T0, strict patient+assay+specimen crosswalk, nonnegative finite area, exact wide-vs-long comparison, clinical/PGx matching, valid observed day1 PK n13 and day9 n12. **All 82,537 (=6,349×13) original MS1 area cells agree between wide/long imports.** Of 6,349 peaks, 6,089 are nonzero in at least one predose plasma sample.

The real-world original midazolam oral CL/F benchmark was reproduced independently **from the private MultiOmics-ready ZIP**, no GNPS full QIIME2 ZIP, with original 13 and 12 same-person PK outcomes. Numeric MAE is in original clinical metadata's **unresolved physical units**:

| Actual original endpoint | Independent humans | Clinical weight+sex | Clinical+CYP3A4/5 PGx | + global MS1 PCA2, train-fold-only | Delta omics vs PGx |
|---|---:|---:|---:|---:|---:|
| Day1 oral midazolam CL/F | 13 | 6.875329 | 7.343413 | 9.888019 | **+2.544607, worse** |
| Day9 oral midazolam CL/F | 12 | 3.984289 | 4.601193 | 4.748841 | **+0.147648, worse** |

This is **real molecular measures linked to true systemic PK**, not merely in vitro CYP evidence. The addition of *global unsupervised* MS1 PCA2 **did not improve PK prediction**. This does **not** refute all biologically constrained metabolomic covariates, RNA networks or mechanistic PopPK. Very small n, partial PGx panel, insufficient physical PK unit certainty, preprocessed GNPS feature universe, internal LOPO rather than genuine independent validation.

## What has NOT been tested (important)

The private input files have syntax compatible with the documented MultiOmics metabolomics CSV and MS peak-AUC long import. But **this does NOT prove an actual browser execution or that the actual browser-produced downstream exchange preserved sample IDs, numerical intensities or annotation identities**.

The new `--native-exchange` option is an explicit second gate that requires a **genuine MultiOmics-generated `exchange/`** with the versioned DOWNSTREAM_EXCHANGE_V1 manifest + SHA256, and compares its **actual numeric values and exact original subjects/assays** to the private input. It refuses a swapped-patient join, modified feature intensity, missing feature, mismatched SHA, incompatible privacy scope, or model-inferred molecule identities. NO substitute local hand-built export is described as a validated native producer.

The MultiOmics tool's QC/filtering and normalization may lawfully change a matrix; if that happens, the native-export parity gate will intentionally **abstain** until an explicit transformations/provenance rule and training-fold-only molecular feature processing are specified. T0 is a relative clinical study baseline, **not an ISO absolute sampling timestamp**; never invent a collection date to satisfy another PK matching contract.

## Exploratory biochemical class proxy (NOT a source-confirmed CYP marker)

The GNPS original MS2 spectral match output contains 346 putative identification rows, and 9 matched MS1 feature IDs with names matching **cholic/chenodeoxycholic/deoxycholic/glycocholic/taurocholic acid** as a predeclared illustrative *bile acid class* proxy. All nine had nonzero signals for at least 10/13 individual baseline plasma samples. Spectral matches can be ambiguous and include overlapping adducts/isomers, so this is NOT 9 independently confirmed metabolite molecules. A uniformly oriented score of source-normalized log-peak areas, feature-z standardized using train people only, was evaluated with the **same person folds** and fixed PGx Ridge reference:

| Endpoint | Baseline clinical+PGx MAE | + exploratory MS2 bile-class proxy MAE | Delta |
|---|---:|---:|---:|
| Midazolam day1, n13 | 7.343413 | 9.975204 | **+2.631791 (worse)** |
| Midazolam day9, n12 | 4.601193 | 5.119382 | **+0.518189 (worse)** |

**This was a second exploratory analysis, NOT an independent test of a preregistered confirmatory hypothesis.** Correlation does not imply direct CYP3A catalysis; an exploratory bile-acid class is not a measured CYP3A flux or liver enzyme activity. No positive biological PK claim. Future score groups need identities assigned using authentic standards or high-confidence spectral annotation, pathway/organ-specific mechanism predefinition and validation on a **separate cohort**. SWEET is an RNA network method, **not applicable** to unidentified MS1 peak areas.

## Reproduction in a private, authorized environment

```bash
python -m pip install -e '.[benchmark]'
python -m scripts.validate_private_jarmusch13_multiomics_pipeline \
  --private-input-zip /private/jarmusch13_private_MultiOmics_Omics2PK.zip
# Only AFTER downloading the *real* MultiOmics browser export:
python -m scripts.validate_private_jarmusch13_multiomics_pipeline \
  --private-input-zip /private/jarmusch13_private_MultiOmics_Omics2PK.zip \
  --native-exchange /private/multiomics_analysis/exchange
```

The command prints **aggregate counts and error metrics only**. Both zipped original source and pseudonymized participant-level input remain strictly outside GitHub/Actions (no upload, artifact, or public fixture). Public CI tests use invented subject/sample/peak values. Source provenance checks alone do not certify ethics/consent, nor validate a research claim of omics superiority.

## Next decisions

1. Run the actual uploaded CSV + metadata in MultiOmics browser, export the actual downstream V1 ZIP, and privately supply only the export for numeric parity. **Do not publish participant-level data in a public issue or PR.**
2. Develop genuine omic-normalization/metadata transform provenance and a fold-honest numerical score loader, not arbitrary DE group fold changes.
3. Preselect biologically meaningful CYP/transport/metabolic pathways and identities without inspecting PK outcomes; seek a larger external PK-plus-omics dataset (busulfan132 and Jarmusch independent cohort are more informative than repeated tuning n13).
4. Only then attempt M3 verified PK concentration-time NLME and future Omics2PD. Clinical dosing is out of scope.
