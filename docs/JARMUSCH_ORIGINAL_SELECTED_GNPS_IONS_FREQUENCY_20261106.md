# Actual original Jarmusch molecular ion selection (per-fold), not just k counts

**Source:** real original GNPS 2019 task b3966445abad4658a3cdc63c8198a666, donor-qualified n10
person-heldout nested \`scripts/benchmark_original_jarmusch_urine_t0_nested.py\`, updated
library-feature lookup \`scripts/benchmark_original_jarmusch_nested_pk_germline.py\`.
**Reproduction:** [GitHub Actions original source #42134309314, success](https://github.com/rberrah/Omics2PK/actions/runs/42134309314).

**The selection counts below refer to independently trained patient-heldout models.**
They are not independent studies, patient-level recommendations, identities confirmed
against analytical reference standards, or clinical laboratory targets.

Source feature references below are exact **original GNPS public MS1 cluster/scan IDs**,
not individual donor IDs. Library name refers to source putative MS2 match, not exact
chemical identity. Original feature frequency is computed after source-verified train-only
selection in each outer donor-LOPO fit. Repeated selections of distinct ions with the
same tentative spectral library hit do NOT imply distinct biologically confirmed
metabolites (adducts/isotopes/structural analogues may explain duplication).

## Candidates from models that numerically improve heldout original day9 CL/F

| Source specimen | Pool | Actual original GNPS feature | Selected across outer models | Source MS2 library match | Claim |
|---|---|---|---:|---|---|
| Day1 predose **plasma**, n10 | 6349 untargeted | **6171** | **4/10** | None passing source annotation quality criteria | Candidate to identify/verify, not CYP function |
| Day1 predose plasma | 6349 | **936** | 2/10 | None | Candidate to identify/verify |
| Day1 predose plasma | 6349 | 3103, 3337 | 1/10 each | None | Selection unstable |
| Day9 predose **urine**, n10 | 6349 untargeted | **369** | **6/10** | None | Candidate to identify/verify |
| Day9 predose urine | 6349 | **3018** | **4/10** | None | Candidate to identify/verify |
| Day9 predose urine | 6349 | **4321** | **3/10** | None | Candidate to identify/verify |
| Day9 predose urine | 6349 | 6135 | 2/10 | None | Candidate to identify/verify |
| Day9 predose urine | 6349 | 410, 414, 1025, 1150, 1786, 1948 | 1/10 each | None | Selection unstable |

**Clinical outcome**: exactly same observed day9 oral midazolam CL/F. Baseline
weight + sex + source CYP3A4/5 PGx **MAE 3.356567**.
J1 plasma untargeted **3.081850** (8.19% numeric reduction,
descriptive paired bootstrap ΔMAE 95% [−0.790684,+0.141669]),
J9 urine untargeted **3.227443** (3.85% reduction,
descriptive ΔMAE 95% [−0.470548,+0.240795]).
**Both intervals include zero**, n10, exploratory multiple assays, no external replication.
Not proof of a functional biomarker or tacrolimus relevance.
In source J9 urine untargeted, k3 was selected in **8/10 outer folds**
(k1 1/10, k0 1/10).

## Source-putatively annotated molecular groups selected but DID NOT improve PK

| Matrix, model | Original feature ID | Fold frequency | Original GNPS source putative library match | Model final MAE (reference 3.356567) |
|---|---|---:|---|---:|
| Day1 urine, 18 steroid/bile candidate ions | **740** | **4/10** | Chenodeoxycholic acid, MQScore .7801 | **3.794865**, worse |
| Day1 urine, 18 candidate ions | **1889** | **4/10** | Same putative chenodeoxycholic acid, MQScore .7801 | same panel |
| Day1 urine, 18 candidate ions | **2496** | 2/10 | Glycoursodeoxycholic acid, MQScore .7887 | same panel |
| Day1 urine, 18 candidate ions | 739, 742, 2497 | 1/10 each | Putative cholic acid / glycochenodeoxycholic acid | same panel |
| Day9 urine, 18 candidate ions | **453**, **739** | 1/10 each | Putative glycodeoxycholic acid, cholic acid | **3.976507**, worse |
| Day1 plasma, 18 candidate ions | **454** | 2/10 | Putative lithocholic acid, MQScore .8977 | **3.535090**, worse |
| Day1 plasma, all source MS2 annotated | **6205** | 7/10 | Putative 1,2-dioleoyl-sn-glycero-3-phosphatidylcholine | **3.991813**, worse |
| Day1 plasma, all MS2 annotated | **258** | 5/10 | Same putative phosphatidylcholine | same panel |
| Day1 plasma, all MS2 annotated | **339** | 4/10 | Putative N-alpha-acetyl-L-arginine | same panel |

**Why included in candidate pools?** The 18-ion panel was **not a verified CYP3A
enzyme-function panel**. It was only an exploratory **GNPS annotation-keyword
filter** for steroids and bile acids after MS2 match quality threshold
MQScore >=.7, shared peaks >=6, absolute mass error <=20 ppm. Being a bile
acid does not mean its measured concentration reports CYP3A5 expression or
specific CYP3A4 catalytic activity. The 207 annotated pool accepts
all quality-passing source MS2 library assignments, regardless of
biological CYP3A relevance. The 6349-MS1 pool has no source chemical
identity requirement and is best understood as *discovery*.

**Why k0 in some folds?** The baseline model always competes against
k1/k3. A positive best-nonzero minus best-zero fourfold inner mean
absolute logCLF error forces no molecular features to be added.
This is *estimator-level* non-retention, not a declaration that the
compound has no molecular activity. In the original day1 midazolam
analysis with n13, k0 wins 13/13 for both 18 and 207 source panels,
despite all 18 source candidates being detectable predose.

## Scientific next step before discussing enzyme mechanisms

To interpret GNPS source ions 369, 3018, 4321, 6171 and 936:
1. Retrieve original **precursor m/z and retention time** and original MS2 spectra, and assess whether any ions reflect blanks, batch, isotopes, adducts or duplicated compounds.
2. Search GNPS/MassBank/analytical library annotations, request actual experimental standards, use validated targeted LC-MS/MS and precise sample processing/timing. A library match alone is not identification.
3. Construct candidate biochemical substrate/product relationships and check actual enzyme/subcellular/tissue localization, relevant negative-control pathways, renal excretion and microbiome confounding.
4. Perform **CYP3A4/5 functional verification**, not only original-model SHAP or posthoc ranking. No mechanism claim until independently observed.
5. Freeze a new qualified panel. Original n10 observations have already informed selection, so cannot be reused as unbiased confirmatory validation.

**Source-sampling caution:** urine pre-dose J1/J9 reflects an excreted,
concentrated and strongly kidney/timing-sensitive matrix; while plasma
has better clinical applicability, neither directly measures liver
protein. Original Jarmusch genotype was not sequenced into full allele
haplotypes and transcript/protein not measured from liver tissue.
