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


## Actual source audit results (2026-10-16)

Reproducible GitHub Actions evidence:
https://github.com/rberrah/Omics2PK/actions/runs/38911376280

- MixOmics drug/ABC cell lines: **60**.
- IntLIM sample metadata: **58**; RNA matrix: **57**; metabolite matrix: **58**.
- **11** exact four-way string matches.
- **48** four-way matches using solely punctuation and uppercase
  normalization with explicit no-collision testing.
- **12** source drug lines remain unresolved by this syntactic operation:
  A549, H226, H23, H322M, H460, H522, HL60, MCF7ADRr,
  MDAMB231, MDAN, RXF393, SKMEL2.

The unmatched IntLIM identifiers include A549/ATCC, HL-60(TB),
MDA-MB-231/ATCC, NCI-H226, NCI-H23, NCI-H322M, NCI-H460,
NCI-H522, NCI/ADR-RES and RXF-393. These require a
**source-verifiable per-cell-line mapping**, not an invented merge.
Particularly, historically changing cell-line identities and
reclassifications (e.g. putative breast-derived resistant lines)
must not be assumed to be synonyms.

**Decision:** 48 syntactic matching line IDs make a guarded exploratory
joint cohort technically plausible. This does not establish that the
gene-expression and metabolomic measurements were collected on the
same biological specimen, at the same time, or in equal culture
conditions as the GI50 experiment. No predictive model has been
run on the joined data. A reviewed identity crosswalk and
provenance/cancer-type stratification remain blocking gates.
