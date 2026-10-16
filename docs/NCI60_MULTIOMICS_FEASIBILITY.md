# Public dual-omics feasibility gate — NCI-60 cross-source identity audit

**Purpose:** determine whether separately hosted public NCI-60 tumor
cell-line expression/metabolomics can be safely aligned to measured
doxorubicin drug response on the same **cell-line identities**.

Sources are pinned:
- mixOmics multidrug: 60 NCI-60 tumor cell lines; 48 ABC RT-qPCR genes;
  measured drug GI50, 1,429 compounds.
- Mathelab IntLIM formatted NCI60: public gene-expression and
  metabolomics matrices plus cell-line metadata; associated IntLIM
  manuscript and vignettes. This repository does not contain a
  directly usable drug outcome column in its sample metadata.

**This is not a same-time or same-sample molecular integration.**
Measurements derive from different studies/technologies and may
reflect different experimental conditions. Exact *cell-line name*
matching alone is the first necessary condition, not evidence of
temporally paired assays or independence.

The gate reports **both exact literal intersections** and the intersection
after a transparent punctuation/case-only key
`toupper(gsub("[^A-Z0-9]", "", name))`. All normalized-key
collisions are fatal. This is **not a biological alias map**.
Fuzzy string matching, unverifiable aliases, row-index joins,
imputation, modality selection, and drug-response fitting remain forbidden
until a source-verified crosswalk and study-provenance manifest exist.
A value of >=30 collision-free syntax-only matches demonstrates only
feasibility of an *exploratory* identity crosswalk, not same-time
multimodal measurement or valid outcome modeling.

Future qualified data harmonization must explicitly distinguish:
1. Within-cohort gene expression measured in different assays (NCI60
   qPCR ABC vs separate microarrays).
2. In-vitro drug GI50 scores (not clinical PopPK), with measured
   duplicates, missingness and sampling times.
3. Tumor class as a crucial confounder. Train/test by *cell line*,
   not assay or inferred cell label.
4. Feature preprocessing inside every outer training fold, with
   held-out cell lines excluded from data integration and SWEET
   network reference construction.
5. Compare transcriptomics-only, metabolomics-only, both, and SWEET
   features with/without evidence priors, *without choosing positive
   outcomes post hoc*.
