# Jarmusch14: Qualifying original GNPS feature abundance matrix, after real oral PK/PGx source validation (2026-10-18)

**Clinical PK and CYP genetics are already genuinely measured and source-qualified** on n13–14 volunteers, DOI [10.1111/cts.12785](https://doi.org/10.1111/cts.12785), author original research repository pinned SHA `95586620ead4ac5241dc1e99a98e047c1f99b17b`. Our [real genotype vs prior probe-function observed PK LOPO benchmark](JARMUSCH14_PUBLIC_REAL_CLINICAL_PK_PGX_20261018.md) is already complete, with original authentic clinical metadata SHA validated. But original participant-level **untargeted baseline MS ion-area** matrix is **NOT** yet imported, so this is currently a functional-phenotype vs genetic comparator, not Omics2PK molecular-abundance-vs-genotype superiority.

## Precise original GNPS task and source file found

The original publication [source article](https://pmc.ncbi.nlm.nih.gov/articles/PMC7485946/) provides **original feature-based molecular networking (FBMN) GNPS task** [b3966445abad4658a3cdc63c8198a666](https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=b3966445abad4658a3cdc63c8198a666), completed 2019-07-09, archived workflow DONE. The public GNPS status page explicitly lists:

- **Feature Quantification Table**: `alan_jarmusch/DrugMetabolism_Tsunoda/MS1basedMolNets_gapfilled_filtered/GNPS_quant.csv`
- **GNPS MS2 MGF**: same original folder, `GNPS.mgf`
- **Original GNPS metadata**: same original folder, `FULL Metadata.txt`

The [author's source R preprocessing](https://github.com/alan-jarmusch/UntargetedMS-DrugMetabolism/blob/master/Code/20191030_UntargetedMS_DrugMetabolism.R) reads `GNPS_quant.csv`, formats and joins peak intensities to `Materials/Metadata.txt` by ORIGINAL study MS run filename. `GNPS_quant.csv` is **not checked into the author GitHub code repository**, but GNPS status proves it was an original 2019 input at the source location. Raw spectra deposited in [MassIVE MSV000082493](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=b8ffc66bab8d435aa863a0fed8fff6fa), which is described as **Partial Public**, and other blood ion sources may be in [MSV000084008](https://www.omicsdi.org/dataset/gnps/MSV000084008).

## Real read-only public GNPS API verification executed

New `scripts/audit_jarmusch_gnps_quant_original_public_task.py` and separate GH Actions CI <https://github.com/rberrah/Omics2PK/actions/runs/39171187796> checked exact documented original public task APIs without login or permissions bypass. No participant names/IDs, metabolite areas, source filenames beyond known official path or original sample records were printed/uploaded.

| Exact original API attempt | Actual result |
|---|---|
| `ManageParameters?task=b3966445abad4658a3cdc63c8198a666` | HTTP original XML returned, **6,394 bytes**, source SHA256 `4740120b61af355c180c0386e241437f55275cb923963c59d376dd2c9b6b5bd0`; verifies job parameters, NOT ion matrix |
| GNPS original `DownloadResultFile` original task, `quantification_table_reformatted/` output | **HTTP 500** on GitHub Actions runner; no peak intensities retrieved |
| GNPS original `DownloadResult` `view=download_cytoscape_data` | **HTTP 405** on GitHub Actions runner; no spectral feature abundance file retrieved |

These HTTP responses are **access limitations for attempted documented endpoints**. They do NOT establish that the authors never deposited or that no researcher can retrieve `GNPS_quant.csv` through GNPS browser/my-files interface or a different authenticated/authorized download interface.

## Specific public data assistance the researcher could offer (NO source patient identifiers in public GitHub)

If the user has a GNPS/MassIVE account and can access the ORIGINAL public completed job, obtain **only the processed feature-by-analytical-sample peak area matrix `GNPS_quant.csv`** from the above source path, **not multi-GB raw mass-spec `.d`/`.mzXML` spectra**. Upload the processed original matrix to the private chat or provide a legitimate non-expiring original public source download URL. If direct access blocked, source-access assistance from the original GNPS submitting user/institution may be required. Do not share passwords or access tokens.

Once obtainable, validate SOURCE sample column headers against original `Materials/Metadata.txt` individual `filename` and source `ATTRIBUTE_Subject`/study day and minutes:
- **Only pre-dose T0 plasma MS** source samples (original published 13 original volunteers on day1, **12 on day9**), not post-drug exposure or study fecal/skin data, as predictors for prospective day9 PK;
- Train/reference `sex+bodyweight+CYP2C19/CYP3A4/CYP3A5 typed SNPs`, then `+predose T0 metabolomic PCA or restrained signature`; fixed held-person-group CV, primary next target actual midazolam `CL/F` or CYP2C19 omeprazole/5-OH metabolite AUC ratio, same original complete subject intersection, and no original-day9 probe PK in training of test subject;
- Functional comparator `+prior clinical day1 probe PK` evaluated as a separate ACTIVE function phenotype, **not** incorrectly called independent multi-omics;
- Source experiment n≤13, dimensionality huge: external clinical replication and cautious uncertainty obligatory. No claims of validated multi-omics improvement if adding thousands of peak intensities worsens out-of-person performance.

**Status:** `REAL_CYP_PGX_AND_CLINICAL_PK_14SUBJECT_DATA_READY; EXACT_ORIGINAL_GNPS_FEATURE_TABLE_SOURCE_PATH_LOCATED; PUBLIC_API_GNPS_XML_OK_BUT_PEAK_MATRIX_HTTP500_405; OMICS_QUANTITATIVE_MATRIX_ACTUALLY_DOWNLOADED_FALSE`.
