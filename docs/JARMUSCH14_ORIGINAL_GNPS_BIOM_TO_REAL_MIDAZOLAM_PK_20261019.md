# Real GNPS QIIME2-MS1 baseline metabolomics → observed human CYP3A PK (19 Oct 2026)

## Source supplied **privately** by project owner

**A new breakthrough in actual accessible data**, not an inference from the status page. Researcher provided a real ZIP download from original completed GNPS feature-based molecular-networking job `b3966445abad4658a3cdc63c8198a666`:

- Original user-uploaded archive `ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-b3966445-download_cytoscape_data.zip`, ~88 MB, exact SHA256 `2047ea2032d73e3b48e2696fb1f652a7dedaa6f8cd5d40d0926d810bdd46c026`. **Archive intentionally absent from public GitHub, CI artifacts and releases.**
- Original GNPS `params.xml` proves author's **input** `GNPS_quant.csv`, but original CSV itself is still not contained in this ZIP.
- Actual GNPS **output** `qiime2_output/qiime2_table.qza` contains an HDF5 BIOM `feature-table.biom`, with **6,349 unique MS1 features × 1,198 analytical sample IDs and 5,215,628 nonzero entries**. Unlike the original input, this is a GNPS/2019-QIIME2 **processed output**. Its peak areas are not verified to be numerically identical to the original `GNPS_quant.csv`, nor guaranteed independent of the original all-study feature alignment.
- Source `qiime2_output/qiime2_metadata.tsv` contains **1,440 rows** but only 1,202 distinct MS analytical sample IDs because several sample names repeat. Verified identical patient/day/time, genotype and PK source values across repeated metadata copies; 1,198 unique sample IDs match BIOM exactly. Do **not** count analytical rows/replicates as independent people.
- Filter **MS_SampleType=sample + blood_plasma + day 1 + time 0 min**. Yields **13 distinct authentic biological subjects** with **13 unique original plasma MS1 assays**. The original study explicitly describes time 0 as **before administration** of the drug probes. The same subjects have true day1 targeted midazolam oral CYP3A CL/F and original CYP3A4/CYP3A5 genotype calls, baseline sex + measured weight. **12 of these 13** have genuine day9 oral CL/F after cefprozil.
- 6,089/6,349 BIOM MS1 features have nonzero measurement in at least one qualifying original predose plasma specimen. The biom abundance values are nonnegative finite ion-signal measures, not molar tissue enzyme activity. Original physical `CL/F` table units are not yet securely verified and errors below use original metadata's numeric scale.

Original open-access methods Jarmusch et al. DOI [10.1111/cts.12785](https://doi.org/10.1111/cts.12785); original GNPS task [b3966445...](https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=b3966445abad4658a3cdc63c8198a666).

## First *ACTUAL original-source in-vivo PK* LOPO run (13/12 humans)

Executed **locally on authentic user-supplied ZIP**. Completely fixed research comparison:

1. Baseline: observed source patient body weight + sex, Ridge alpha=1 predicting log(midazolam oral CL/F).
2. Baseline + true source CYP3A4 and CYP3A5 categorical calls one-hot encoded on the training subjects only. It is **a limited PGx panel**, not complete phased CYP genetic variation.
3. Same baseline+PGx + 2 **train-only** unsupervised principal components from original day1 T0 GNPS peak intensity vectors. Data preprocessing: each person's summed peaks → relative ion intensities × 1e6 then log1p; in each fold keep MS1 features present in >= 50% of **training** people, standardize on train, fit PCA2 on train, fit fixed Ridge(alpha=1). No post-dose MS samples or same-person followup PK outcome enters features. One original held-out person per fold; no source IDs or person-row predictions printed/committed.

| Actual measured target | Independent persons | Sex+weight MAE | Sex+weight+PGx MAE | +PGx+predose-MS1 PCA2 MAE | Omic ΔMAE vs PGx |
|---|---:|---:|---:|---:|---:|
| Day1 midazolam oral CL/F, original numeric scale | 13 | **6.875329** | 7.343413 | 9.888019 | **+2.544607** (+34.65%, worse) |
| Day9 midazolam oral CL/F, original numeric scale | 12 | **3.984289** | 4.601193 | 4.748841 | **+0.147648** (+3.21%, worse) |

Exploratory *paired errors*, 1,999 bootstrap resamples WITHOUT refitting any molecular preprocessing/model, **not prospective CI**:
- Day1 incremental omic ΔMAE 95% descriptive interval **[+0.4421, +5.9959]**, only 1/13 individual errors improved over PGx baseline. 25 shuffled person-to-metabolome negative controls had median omics-augmented MAE **8.0924**, all better than actual omics MAE (empirical 25/25).
- Day9 incremental omic ΔMAE 95% interval **[-0.4273, +0.6899]**, 4/12 individual errors improved. 25 shuffled omic-person negative controls had median omics-augmented MAE **4.6313**, 13/25 at-or-below actual model error.

### Correct falsifiable interpretation

This is a **genuine paired original human omic-to-in-vivo clinical PK prediction analysis**, technically more advanced than liver CYP activity evidence. It **DOES NOT** demonstrate improvement over genuine source CYP3A4/CYP3A5 variants (on either day). A massive unsupervised untargeted 6,349-feature MS1 space projected onto just 2 global PCs in n12–13 patients is a weak, nonspecific metabolite covariate; the negative result does not establish that mechanism-targeted metabolite scores, stronger genetically stratified data, SWEET *RNA* networks, PopPK or multiomic variables cannot improve PK elsewhere.

**Important methodological limits:** cohort tiny, many MS1 ions unidentified/isobars, no biological mechanism claim, source preprocessing MZmine/GNPS across the complete source assay set *before* LOPO introduces possible unsupervised transductive feature-universe dependence, no hold-out study, no externally frozen normalization/reference, original CYP3A CL/F measurement units need verification, incomplete pharmacogenetic star-allele inference, multiple exploratory outcomes, no observed C(t)-likelihood population-PK fit. The 25 molecular-shuffle controls are diagnostic only and not enough for exact statistical significance.

This source meets E3 **clinical real-PK linkage and baseline comparison**, but does **not** satisfy evidence of omics superiority or external validation. It is metabolomics (NOT transcriptomics), so do **not** apply RNA-only SWEET to these peaks or claim measured per-patient CYP enzyme activity. Main next step: preregister endogenous drug-metabolism pathway **candidate metabolite classes** using source spectral annotation, matching published biology/organ provenance; compare with nonspecific PCA on a genuinely independent sufficiently sized study before a strong conclusion.

## Reproducibility without publishing source data

Source reader: `scripts/validate_local_jarmusch14_gnps_biom_pk.py` (source-byte SHA-gated, local only, prints aggregate only). For an **authorized local** researcher with the exact archive:

```bash
python -m pip install -e .
python -m pip install numpy==2.3.5 scikit-learn==1.7.2 scipy==1.16.3 h5py==3.15.1
python -m scripts.validate_local_jarmusch14_gnps_biom_pk \
  --archive /path/to/ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-b3966445-download_cytoscape_data.zip
```

Only artificial unidentifiable fixtures are used in public GitHub Actions. Uploaded exact raw GNPS ZIP, extracted BIOM, patient source metadata, original per-person MS values and PK source rows are **not** committed or uploaded to workflows. Research results are historical exploratory internal validation, not clinical decision support.
