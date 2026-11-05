# Public EATRIS-Plus v4 original WBC versus plasma omic integration — no clinical PK

**Study source:** EATRIS-Plus multiomic reference cohort, official
[Zenodo v4 dated 2026-03-20](https://zenodo.org/records/19136013).
Only original **publicly downloadable** data; do not assume future
BIOMARGIN/ABIS, clinical tacrolimus, controlled WGS or other patient
cohorts. This is a **MultiOmics integration/negative-control**
study, not an Omics2PK patient-PGx outcome study.

## Physical source of the omics

The v4 authors describe:
- **White blood cells:** bulk mRNA-seq, mature miRNA-seq,
  EM-seq methylation.
- **Plasma:** protein assays, acylcarnitines,
  amino acids, very-long-chain fatty acids,
  lipidomics; plasma miRNA qRT-PCR.
- **Do not call WBC RNA hepatic CYP3A mRNA**. Plasma
  proteins/metabolites are indirect systemic physiology
  measurements rather than tissue-specific active
  hepatic enzymes.
- These original human blood/plasma profiles have no
  actual **drug PK outcomes** or unrestricted full
  germline pharmacogenetic panel. No complete CYP3A5
  diplotypes or simulated tacrolimus AUC/C0 may be
  appended to fake a source-ready study.

**Source verification:** genuine original public
`EATRIS_PLUS_mudata.h5mu` 141.7 MB, exact official
MD5 `bbf0342797adcbcd6a3c5b242ca64573`,
verified in
[MultiOmics original data-source CI SUCCESS #42036020160](https://github.com/rberrah/MultiOmics/actions/runs/42036020160).
Actual source `/mod/<assay>/obs/_index` was
read and joined across each view by the **original
sample ID**, not sorted row position. Original individual
records were never committed or printed.

## Real cohort and quality of original matched views

Original h5mu **125** global individual identifiers
and **12** named assay matrices. Examples:

| Source original assay | n with source view | Features | Physical source |
|---|---:|---:|---|
| `mRNA` | **90** | 16,386 | WBC |
| `miRNA-mature` | 89 | 2,130 | WBC |
| `EM-seq` | 117 | 100,000 | WBC |
| `proteins` | 125 | 1,393 | Plasma (source-description level) |
| `proteomics-ipp` | 123 | 6,831 | Plasma-related distinct assay (do not equate platform) |
| `acylcarnitines` | 125 | 33 | Plasma |
| `aminoacids` | 125 | 46 | Plasma |
| `lipidomics-neg` / `lipidomics-pos` | 125 / 125 | 164 / 196 | Plasma |
| `miRNA-qRT-PCR` | 86 | 179 | Plasma |

Only **55 original subjects** occur in all twelve matrices.
The useful targeted 4-matrix **WBC mRNA + plasma
acylcarnitines + plasma amino acids + plasma `proteins`**
intersection is **n90**. This is a genuine participant
link; it is **not** proof that all assays came from
one simultaneous blood draw or one plasma aliquot.
A detailed per-assay original ISA biological material
crosswalk should be checked rather than inferred from
feature names or an organ-specific gene label.

## Direct quantitative public E1/analysis result

**Outcome:** reconstruction of original measured
`proteins` plasma profile using source original
WBC mRNA, plasma targeted metabolites, both or
source donor-label permutation. **80 source proteins**
were selected using only the corresponding
outer-training fold. Evaluation MAE on
training-standardized targets averaged over each
heldout independent donor, 3 repeats × 5
**original subject-disjoint** outer folds, exact
same 90 subjects.

| Prespecified source arm | Fixed 3 PCs per modality, Ridge alpha10 | 3 PCs **total** per model, alpha chosen on inner 3-fold CV |
|---|---:|---:|
| No molecular predictors: training mean | **0.760615** | **0.760615** |
| WBC mRNA | 0.783413 | 0.760756 |
| Plasma targeted metabolome | 0.792970 | 0.761520 |
| WBC RNA + plasma metabolites | **0.819299** | 0.761429 (1 RNA + 2 metabolites) |
| WBC RNA + plasma metabolites alternate allocation | — | 0.760720 (2 RNA + 1 metabolite) |
| WBC RNA train correspondence permuted control | 0.784368 | 0.760747 |

Primary original source replays:
[fixed-source #42037048717](https://github.com/rberrah/MultiOmics/actions/runs/42037048717)
and
[matched-PC/nested-alpha #42038211111](https://github.com/rberrah/MultiOmics/actions/runs/42038211111),
both succeeded. Detailed methods, confidence intervals,
negative controls, feature count and data contracts:
[MultiOmics original-source scientific report](https://github.com/rberrah/MultiOmics/blob/feat/variant-omic-functional-evidence-sidecar-20261105/docs/multiomics/EATRIS_PLUS_V4_PUBLIC_WBC_PLASMA_EVIDENCE_20261105.md).

### What this demonstrates

1. **Cross-compartment assays are not interchangeable.**
   WBC mRNA is a direct white-cell transcript measurement,
   not a hepatocyte CYP3A activity assay. Plasma targeted
   metabolites and proteins measure systemic pools, not
   each tissue's catalytic contribution.
2. **Excess information may degrade an *estimated finite-data
   model*.** In the first fixed alpha10 model, joint
   3-layer information degraded protein reconstruction by
   ΔMAE **+0.058684** versus the no-omics training
   mean [source original-donor descriptive 95% CI
   **+0.045340,+0.072923**].
3. **The algorithm can learn to ignore additional views.**
   After matching dimensionality and choosing shrinkage
   by nested train-only CV, the joint arms approached
   no-omics (ΔMAE +0.000814 or +0.000105,
   depending on WBC/plasma PC allocation).
   RNA-containing arms often selected **alpha 100,000 or
   1,000,000** — near-zero molecular coefficients.
4. **No demonstrated positive incremental contribution**
   to this particular measured protein target on these
   available public views and fixed evaluation designs.
   Difference over plasma metabolites alone in
   the capacity-matched models: Δ−0.000091
   [−0.001086,+0.000764] or Δ−0.000800
   [−0.001862,+0.000176]; neither descriptive interval
   excludes zero. The null (permuted donor pairing)
   is similar to the true WBC predictor.
5. **Not an intrinsic information-theoretic penalty.**
   A sufficiently flexible optimal predictive function
   can ignore redundant features. These source
   examples specifically show the cost of exploiting
   them with a finite-sample estimator and
   under-/over-regularization.

**Critical design qualification:** fixed versus nested
arms differ both in **total PC budget** and
**regularization optimization**, so their contrast
cannot isolate which mechanism caused the reduction
in prediction error. No external cohort, site
transport, multiplicity correction, nested-bootstrap
refitting, genuinely measured enzyme activity, clinical
transplant outcome or patient dose utility was obtained.

## Omics2PK research relevance without pretending these data contain clinical PK

This is a *method-selection and biological provenance
control* to be retained in the intended two-way
MultiOmics→Omics2PK pipeline:
- Genotype interpretation must use **actual germline/
  transplanted-donor source** and allele/transcript
  coverage; not available here.
- A real tissue or biofluid specimen origin,
  method and pre/postdose landmark must be
  recorded before a signal can be called
  mechanistically relevant.
- MultiOmics identifies measured
  layer/variant-associated evidence and possible
  discordance, abstaining where source-specific
  functional tests are missing.
- Omics2PK must still compare each view against
  a **real clinical + complete PGx reference on
  independently observed drug AUC/CL/C0**,
  if/when such original data actually exist.
- For now: legitimate public clinical PK exists
  in Rotroff metformin n32 with **predose
  heparin plasma metabolomics**, but not
  a source-qualified full PGx comparator;
  public EATRIS gives *no human PK endpoint*.
  Neither should be relabelled as tacrolimus
  genotype-plus-omics success.

**No new datasets need to be requested to reproduce
these public-only E1 research steps.**


## Addendum — sensitivity to the measured plasma protein panel (2026-11-05)

MultiOmics [research PR #121](https://github.com/rberrah/MultiOmics/pull/121)
now contains a reproducible experimental extension comparing
the two independent source matrices \`proteins\` (1,393 original
features, 125 people) and \`proteomics-ipp\` (6,831 features,
123 people) on their **five-assay common subject intersection**.
The same WBC-RNA + targeted-plasma-metabolite predictors,
person-disjoint 3 × 5 folds, 3/9 total-PC budgets, and fixed
Ridge penalties 10 and 1e6 are used against both targets.

For Omics2PK, this is a **measurement-panel sensitivity
benchmark**: only improvement relative to each panel's
own train-mean no-omics predictor is compared. An absolute
standardized error for the 1,393-feature panel is not
commensurate with one from the 6,831-feature panel.
Furthermore, the source ISA has **not yet established
equivalence of the panel-specific collection times,
aliquots, analytes and platforms**.

**Qualification required before any positive claim:**
(i) successful original-source checksum-verified benchmark
and paired-person checks; (ii) separate material, method,
missingness and analyte-identity audit; (iii) independent
clinical PK and PGx observations if the proposed conclusion
mentions improvement over pharmacogenetics, AUC or clearance.

The **original public-source replay has now succeeded**:
[MultiOmics original EATRIS two-panel CI #42077407584](https://github.com/rberrah/MultiOmics/actions/runs/42077407584),
88 matched donors, exact original source checksum verified.
This establishes stability or instability of a particular
*interomic reconstruction objective*, not that a new omic
layer improves individual drug PK or that either measured
panel measures hepatocyte CYP3A4/5 activity. The two-target
comparison remains explicitly exploratory.


### Original matched-source results: protein target choice changes the size of the negative result

Same 88 original persons, same predictors and person-disjoint
3 x 5 folds; exactly 80 panel-specific protein targets per
outer training fold. Fixed Ridge alpha10, 3 total PCs unless
otherwise specified. Values below are heldout **standardized
plasma-protein profile MAE**, not drug PK.

| Source target | No omics MAE | WBC RNA MAE | Plasma metabolites MAE | WBC + plasma MAE | WBC + plasma 9-PC MAE |
|---|---:|---:|---:|---:|---:|
| \`proteins\` | **0.753535** | 0.776767 | 0.766483 | 0.772035 | 0.805722 |
| \`proteomics-ipp\` | **0.331838** | 0.367851 | 0.361322 | 0.367762 | 0.421688 |

Relative excess error of 3-PC joint omics versus **the
matching target's own no-omics reference** is
**+2.455%** for \`proteins\` and **+10.826%**
for \`proteomics-ipp\`. At 9 PCs these excess
errors rise to **+6.926%** and **+27.076%**.
Ridge alpha=1e6 restores prediction to each
panel's no-omics mean. The original subject-paired
difference in **relative skill** (IPP minus proteins)
for joint 3-PC alpha10 is **−8.3706 percentage
points** [descriptive bootstrap 95%:
**−11.2219, −5.8733**]. The inferential
bootstrap holds fitted CV predictions fixed;
it does not adjust for repeated analysis choices
or multiple comparisons.

**Critical:** The two protein panels quantify
*different targets*, and no exact original
analyte / assay / aliquot crosswalk has been
verified. Thus these observations mean that
**the negative finding depends on which
original measured protein matrix is predicted**,
not that one proteomic technology is inferior.
Nor do they demonstrate that transcriptomics
or metabolomics cannot improve an enzyme
activity outcome, a patient-specific AUC/CL
or genuine clinical+genetic PGx prediction.


### Original-source feature identity diagnostic: no direct shared assay feature keys

[Original revised source replay #42078275344 **SUCCESS**](https://github.com/rberrah/MultiOmics/actions/runs/42078275344)
independently verified the published H5MU checksum, passed
source-ID/target-isolation synthetic tests and audited the
source \`var\` feature keys in both plasma-related proteomic panels.
**Exact original feature-key overlap: 0.**
**Case-folded, whitespace-trimmed feature-key overlap: 0**
(with no resulting case-fold collisions).

This is an **identifier-namespace finding, not proof that the
two assays quantify zero shared biological proteins**. The
source-key strings are not a reliable cross-platform protein
crosswalk; an original analyte ontology/UniProt/gene/proteoform
and assay-method mapping is required before any matched-analyte
or measurement-technology causal comparison. No protein names
or individual donor rows were emitted.

Both replayed panels had complete heldout observations for
the 80 fold-selected evaluation targets in each fold (80
times the heldout fold size), so the very different
baseline standardized MAEs are not explained simply by
unmeasured heldout outcomes *among these selected targets*.
Target distribution, scaling, assay chemistry and analyte
identity remain possible explanations, not established causes.

**Disposition:** keep the negative result as a robust
observed *target-panel sensitivity* (n=88 original persons);
mark any claim that one protein platform is superior to
the other as **NOT ESTABLISHED**. Exact feature identity
alone cannot qualify a common target set.
