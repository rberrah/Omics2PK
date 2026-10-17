# Jarmusch n14 — original publisher supplementary files recovered from official Europe PMC (2026-10-17)

Primary source: Jarmusch et al., *Clinical and Translational Science*, DOI https://doi.org/10.1111/cts.12785, PMC7485946. Human research protocol includes measured quantitative targeted plasma PK (caffeine, midazolam, omeprazole, metabolites), untargeted plasma/urine/fecal MS, CYP genotypes and stool 16S before/after antibiotic perturbation.

## Actual source access achieved in GitHub CI

**Authoritative original endpoint**: https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7485946/supplementaryFiles

**Original article ZIP** SHA-256 `60ea77fbc54833cb049e2a8084befcc2d38d11335b2d12cc5d8822b2d83f0e51`, 1,920,454 bytes. [Source-archive GitHub Actions run](https://github.com/rberrah/Omics2PK/actions/runs/39023183567), status success. `scripts/audit_jarmusch14_europepmc_supplement_zip.py` inspects only archive metadata and original Office structure. No original individual patient rows or PK values logged/redistributed.

**Two publisher supplements originally rejected by direct Wiley access are NOW legitimately retrieved from Europe PMC:**

| Original article supplement | SHA-256 / structure | Scientific significance |
|---|---|---|
| `CTS-13-972-s001.docx` | `9c96a723f6bf012b8a441ea166d3dbf1389eaeb01acba2eb697ab450fa799944`, 1,644,596 bytes, one original **15-row × 11-column** table | Header columns: subject, age, sex, race, CYP2C19 genotype/phenotype, CYP3A4 genotype, CYP3A5 genotype, CYP3A phenotype, CYP2D6 genotype/phenotype; **14 actual participants**. It is a real **individual human pharmacogenomic covariate table**, not quantified patient C(t). |
| `CTS-13-972-s002.xlsx` | `4d84d0ae3a26c6d9f0f46a4804fd87e828d80d9c1a3aabf9d65ff59deafe4ace`, 12,426 bytes, **1 worksheet**, sheet name `Serene_FDR_sOTU_annot.txt`, 30 original rows including header, seven apparent columns; 116 numeric cells | Source contains a compact taxonomic/ASV annotation result, **not a complete raw 14-person × microbial features** matrix or a clinical PK concentration table. |

The full original MassIVE project is https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=b8ffc66bab8d435aa863a0fed8fff6fa (MSV000082493), with publicly available raw MS and original author's `GNPS_quant.csv` processing script https://github.com/alan-jarmusch/UntargetedMS-DrugMetabolism . **MassIVE project must still be checked at individual sample-ID and plasma timepoint level**; original source package alone does not establish a full human PK↔omics crosswalk.

## Evidence boundary

- [x] Original publication and PMC OA metadata verified.
- [x] Original publisher supplements retrieved from official Europe PMC without violating robots controls / scraping.
- [x] Original human n14 **CYP genotype phenotype source** located in full, dimension-qualified, SHA frozen.
- [x] Original microbiome taxon/ASV summary source located, dimension-qualified, SHA frozen.
- [ ] Original independent targeted LC-MS/MS subject×dose×time×concentration (plus LLOQ) **not located as numeric 14-person longitudinal matrix** in supplied Office files.
- [ ] Original pre-dose **individual** untargeted plasma metabolomics intensity matrix (per human subject and assay batch) **not validated with subject crosswalk** in publisher supplements; original MSV000082493 may hold raw files and `quant.csv`.
- [ ] NCA or model estimation from original dose, correct C(t) units and time metadata withheld until patient crosswalk passes.
- [ ] No real E3 predictive uplift or numeric PopPK model can be claimed.

**Consequence**: PGx covariates are authentic and available (locally / subject privacy), so the next true data need is targeted quantitative PK source values and actual ID mapping to untargeted plasma features. The existence of n14 rows in PGx table **cannot** create a PK concentration table. Do not match subjects by row order or interpolate from plot pixels.

Follow-up in issue #44 and primary audit `docs/JARMUSCH14_MASSIVE_TARGETED_PK_CROSSWALK_ORIGINAL_AUDIT_20261017.md`.
