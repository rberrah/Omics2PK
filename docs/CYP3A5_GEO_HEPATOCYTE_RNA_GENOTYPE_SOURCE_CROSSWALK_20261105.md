# CYP3A5 rs776746: qualified public hepatocyte RNA/genotype identity, not yet exon-QTL reconstruction

**Research source audit, 5 November 2026.** This follows the bidirectional genomic variation → MultiOmics → Omics2PK programme. Original source sample metadata were fetched and independently counted in [successful GitHub Actions #41998586213](https://github.com/rberrah/Omics2PK/actions/runs/41998586213). The original RNA-seq and genotyping SOFT metadata were read independently: **this did not create a joint genotype/exon expression table**.

## Published biological effect

Zhong et al., *PLoS Genetics* 2020, DOI [10.1371/journal.pgen.1008662](https://doi.org/10.1371/journal.pgen.1008662), primary human hepatocyte regulatory genomics: CYP3A5*3 / rs776746 is a well-established splice-disrupting allele. The publication reports a CYP3A5 **exon-level expression association (eExon-QTL)** that is not equivalently recovered as a significant gene-level eQTL. **This is a publication-level result, not an independently reproduced per-sample eExon-QTL in Omics2PK.**

The important mechanistic hypothesis is **not** “RNA is uniformly normal in every *3 carrier”; instead it is that **total gene-level quantification can obscure an isoform-/exon-specific splice effect**. Which layer predicts functional CYP3A5 activity or drug PK is still unverified.

## Actual public original source audit

| Official source | Sample records discovered | Matching identity requirement |
|---|---:|---|
| GEO [GSE124074](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124074) primary hepatocyte RNA-seq | **60** original GSM samples | Exact original source sample titles, pre-QC |
| GEO [GSE147628](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147628) Illumina genotyping | **96** original GSM samples | Exact original source sample titles, pre-QC |
| Pairing by matching original sample title | **60** overlaps | **Candidate source-metadata pairing, not validated individual genotype and quantitative transcript join** |

Because GEO does not supply a `series_matrix.txt.gz` for both series, a first source fetch correctly returned **HTTP 404**. The audit was repaired to use GEO's **actual official family SOFT source** for RNAseq/genotyping; it streams the original archive and collects **only accession and sample-title identities in memory**. No individual original donor titles or expression values are committed, displayed in CI, or exported. Original source compressed file SHA-256 hashes are logged in the successful CI run.

The audit **does NOT** establish (and deliberately marks false) each of these:
- Independent original **rs776746 genotype** calls and complete CYP3A5*3/*6/*7 germline diplotypes;
- Original quantitative **gene-level RNA, exon-level RNA, splice junction/isoform RNA** joined on the exact samples;
- Independent **CYP3A5-selective enzyme activity**, canonical/aberrant proteoform or real clinical drug disposition;
- **Tacrolimus exposure, CL/F, C0 or actual PD** in these subjects.

An exact title match is not proof that the genotype arrays and gene expression have been normalized and quality-controlled to the same donor; genetics and repeated technical assays need original source QC.

## Next quantitative experiments, with explicit source gates

**MultiOmics (biological signal/known variant):**
1. Source-retrieve and independently authenticate the actual `rs776746` genotype in the genotyping original array/call matrix, assembly, strand and allele conventions.
2. Source-retrieve and quantify original per-sample **total CYP3A5 gene RNA** versus **exon-/junction-specific CYP3A5 RNA** from the 60 samples. Ensure annotation build, transcript/probe method and exon orientation are matched, and account for different gene isoforms.
3. Only then test actual source-called genotype groups, adjusting for batch/cell model if warranted. Compare original gene-RNA and exon-RNA representations with independent donor-level folds and report whether the known allele-specific splicing effect is detected. This is **RNA mechanism validation**, not yet measured enzyme function.
4. A high total RNA with low canonical isoform warrants the hypothesis of aberrant splicing; **this does not prove proteasomal decay or catalytic inactivity in the person**.

**Omics2PK (drug PK/PD relevance):**
1. Add source-qualified same-donor full PGx and actual measured PK/PD endpoint.
2. Clinical+PGx+PopPK baseline versus total RNA, canonical splice/junction RNA, relevant protein/validated enzymatic activity/metabolites, and no additional omics.
3. Assay- and organ-aware hypothesis selection strictly *inside* training folds; heldout patients and independent study for MAE/RMSE/log-fold exposure bias, clinical calibration and decision metrics.
4. No clinical PK effect is inferable from the current GEO hepatocyte expression/genotype metadata or the publication's eExon-QTL.

## Implementation

- Original public source qualification: `scripts/audit_public_cyp3a5_geo_rna_genotype_join.py`;
- Synthetic source negative controls: `tests/test_public_cyp3a5_geo_rna_genotype_join.py`;
- Successful original GEO replay: [Actions #41998586213](https://github.com/rberrah/Omics2PK/actions/runs/41998586213).

**Scientific status: source-metadata matching SUCCESS; original variant/exon quantitative phenotype and clinical PK NOT TESTED.** This distinction must survive downstream exchange.