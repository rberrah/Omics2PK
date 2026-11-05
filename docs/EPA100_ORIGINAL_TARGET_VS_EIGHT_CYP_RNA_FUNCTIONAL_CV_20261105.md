# Original human liver EPA: does eight-CYP RNA improve functional prediction beyond target RNA?

**Study type:** independently reproduced, exploratory, public-original same-donor analysis, 2026-11-05.
**Run:** [successful original-source GitHub Actions #42101612511](https://github.com/rberrah/Omics2PK/actions/runs/42101612511).
**Data:** original [US EPA DOI 10.23719/1519197](https://doi.org/10.23719/1519197),
published with [Liu et al. 2021, 10.1080/00498254.2020.1867929](https://doi.org/10.1080/00498254.2020.1867929).
The source XLSX file is SHA256 pinned to
\`cfc9dc334064ea61e99ca6fd39c2c8ef7c621d000a37c4da2c95f26ca85e36dd\`.
Donors are aligned by the **two publisher source sample identifiers**, never by
row order or demographics. n=100 original donor bank; per-CYP all-eight-transcript
complete-case n=93–95. Original donor IDs and quantified values stay transient in CI.

## Research question

Can preselected additional CYP qPCR transcripts improve prediction of a
**separately measured functional CYP activity** beyond the matched CYP transcript,
or does the increased molecular dimension make the model less reliable?

**Critical boundary:** qPCR and microsomal enzyme probe rates come from matched
**human liver samples**, not peripheral plasma. This is **targeted single-omic
transcript measurements**, not multi-omics and **NOT observed patient AUC, CL,
drug concentrations, a complete germline genotype test, or clinical dosing**.

## Prespecified computational implementation (but exploratory post-source)

Eight measured CYP endpoints with separate same-donor complete-case analysis:
CYP1A2, CYP2B6, CYP2C9, CYP2C19, CYP2D6, CYP2E1, CYP3A4, CYP4A11.
For every endpoint, identical 3 repetitions x 5 original-donor-disjoint heldout
CV splits, outcome original measured enzyme activity in pmol/min/mg protein.
Training-only feature standardization, outcome log1p fitting, fixed Ridge
\`alpha=1\` or \`alpha=100\`, reverse transform predictions and MAE in observed
original activity units. **No model or molecular feature selected on test labels**.
Original input qPCR is log1p-transformed before train-only standardization.

Four paired model arms:
1. age + sex only;
2. age + sex + **matching CYP transcript**;
3. age + sex + other seven nonmatching CYP transcripts;
4. age + sex + matching + seven nonmatching (all eight).

The complete-case donor intersection is held identical for all four arms of
each CYP. Repeated CV losses are averaged **within unique person** before
1999 subject-level paired residual bootstrap draws. Descriptive 95% intervals
are neither model-refit uncertainty nor corrected for eight related endpoints.
The experiment/alphas were chosen **after viewing other original EPA results**:
there is no pre-registration or independent validation.

## MD5/SHA qualified actual original-source results

**Ridge alpha=1, MAE in each enzyme assay's own units (lower is better):**

| Functional endpoint | n true donors | age/sex | + own CYP RNA | + other 7 only | + all 8 RNA | Relative change **all8 vs own CYP** |
|---|---:|---:|---:|---:|---:|---:|
| CYP1A2 | 93 | 30.673469 | 25.318876 | 29.886436 | 24.351345 | **−3.82%** |
| CYP2B6 | 95 | 59.372846 | 48.726999 | 51.950211 | 51.329090 | **+5.34%** |
| CYP2C9 | 95 | 634.415749 | 614.495808 | 627.868273 | 635.868082 | **+3.48%** |
| CYP2C19 | 93 | 69.842287 | 69.116785 | 75.211538 | 72.440369 | **+4.81%** |
| CYP2D6 | 95 | 138.411519 | 117.798583 | 147.290322 | 128.536874 | **+9.12%** |
| CYP2E1 | 95 | 970.637830 | 986.221007 | 975.745888 | 988.093599 | **+0.19%** |
| CYP3A4 | 94 | 2080.331982 | 1766.969175 | 1808.334515 | 1755.430270 | **−0.65%** |
| CYP4A11 | 94 | 618.514415 | 606.155002 | 574.631818 | 583.810054 | **−3.69%** |

For fixed alpha1, adding all seven other transcripts worsens MAE
for **5/8** activities and improves **3/8** versus the model with
only the matched transcript. This is not independent eight-cohort
replication; all endpoints come from the same human liver bank.

**CYP3A4 specific:** adding matched CYP3A4 qPCR decreases
MAE **2080.332 → 1766.969** vs age/sex (**15.06% reduction**).
Including all eight mRNAs yields **1755.430**, only a further
**0.65% relative reduction**. The exploratory same-person
descriptive 95% bootstrap interval for the MAE difference
(all8 minus own target) includes zero:
**[−168.929,+142.326] pmol/min/mg protein**.
No generalizable gain from the seven extra transcripts is established.

For fixed alpha100, relative extra-seven-RNA effects differ substantially:

| Functional endpoint | own CYP RNA MAE, alpha100 | all8 RNA MAE, alpha100 | Relative change all8 vs own CYP |
|---|---:|---:|---:|
| CYP1A2 | 27.677 | 27.280 | −1.43% |
| CYP2B6 | 53.539 | 51.200 | −4.37% |
| CYP2C9 | 606.993 | 587.678 | −3.18% |
| CYP2C19 | 70.041 | 70.884 | +1.20% |
| CYP2D6 | 127.106 | 128.783 | +1.32% |
| CYP2E1 | 968.114 | 964.256 | −0.40% |
| CYP3A4 | 1862.037 | 1794.116 | −3.65% |
| CYP4A11 | 601.606 | 574.991 | −4.42% |

Thus a substantially stronger prespecified ridge penalty
changes **which CYP endpoints** seem to benefit from
the broader panel: the observed penalty from additional
features is model-/regulation-sensitive, not a universal
information-theoretic penalty. Alpha100 is **not selected
in a nested model-choice protocol** here, so do not cherry-pick
the better alpha per CYP and advertise the resulting accuracy.

A positive nominal paired interval for worsened CYP2D6
(all8 vs targeted RNA, alpha1; ΔMAE +10.738,
descriptive 95% [0.920,21.096]) is **not** adjusted
for the eight correlated outcomes, both fixed alphas,
and post-source question selection. Do not call this
an individually confirmatory significant effect.

## What can be claimed in an application for data or funding

**Qualified empirical evidence:** a same-donor *human tissue function*
bridge can be tested reproducibly using source-verified data;
CYP3A4 target transcript shows a substantial held-donor
predictive improvement beyond age/sex in this specific
in-vitro activity assay. Broadening the predictor panel to
eight CYP transcripts does not necessarily add benefit;
its performance depends on model shrinkage.

**Necessary clinical hypothesis still unproven:**
genotype + pre-dose circulating molecular signals improve
future tacrolimus AUC/CL prediction beyond the best
clinical + genuine complete CYP3A5/3A4 genotype/PopPK
baseline. That requires same-person patient PK and valid
tissue, assay/time, mutation and specimen provenance.

**Minimum proposed prospective / banked-sample pilot:**
authentic individual clinical+PGx and prior-dose/time-stamped
tacrolimus PK, banked **pre-dose plasma** for endogenous
targeted metabolomics, recorded molecular/renal/hepatic
covariates (inflammation, co-medications, hematocrit),
unrelated 1st-dose and repeat prediction validation.
Optional separate PBMC/EV/proteomics remains
**mechanistically unvalidated liver proxy**. Split by patient
and clinical center/time if possible, nested omic selection
and penalty tuning and external replication. Primary outcome
held-out observed independently measured AUC/CL/F,
not model-predicted AUC reused as a tautological target.

The source-controlled code is
\`scripts/benchmark_original_epa_cyp_8rna_saturation.py\`,
tests \`tests/test_epa_cyp_8rna_saturation.py\`
and workflow \`.github/workflows/epa-original-eight-rna-target-saturation.yml\`.
