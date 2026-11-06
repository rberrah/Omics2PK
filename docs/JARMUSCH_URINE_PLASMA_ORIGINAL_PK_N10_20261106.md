# Jarmusch original public paired urine / plasma predose and measured day9 midazolam PK
_2026-11-06 · source-authenticated original n10, exploratory posthoc · no tacrolimus evidence_

## Source, sampling and source-consistency
- Jarmusch et al., CTS 2020 https://doi.org/10.1111/cts.12785
- Public GNPS official task \`b3966445abad4658a3cdc63c8198a666\`, canonical SHA256 of metadata, BIOM and task params checked. Source provenance \`scripts/audit_original_jarmusch_compartments.py\`.
- Real original **14 urine day1 time0 / 14 urine day9 time0** (urine before day-specific oral midazolam), 13 day1 plasma time0, 12 day9 plasma time0. Of 11 people with two plasma predose collections, **10** have **both urine predose samples**, both plasma predose samples, original clinical+genotype and genuine day9 oral midazolam CL/F. All tested matrices rejoined on the **same original biological person and assay source IDs**, not by row position.
- Original urine collection after midazolam at 2h/4h/6h/8h exists but is **not** analyzed as a predictor for current-visit CL/F because midazolam and its downstream drug metabolites would leak the endpoint. Baseline day1 urine to day9 CL/F is an 8-day prognostic prediction before cefprozil; day9 urine to day9 CL/F is contemporary predose signal after cefprozil.
- Urine untargeted MS1 is **relative peak area** (per-sample rowsum normalization); urinary creatinine, osmolality, diuresis / timed urine excretion rates and exact confirmed steroids not documented in this model. Spot urinary cortisol ratio without concentration standards, circadian time, renal adjustment and independent chemistry should not be equated to an endogenous CYP3A function test.
- Genuine source CYP3A4/5 genotype is limited categorical PGx, not full CYP star allele panel. Baseline additionally has **body weight and sex** (NOT age).

**Reproducible real source CI:** [#42133661009, successful](https://github.com/rberrah/Omics2PK/actions/runs/42133661009). Source code \`scripts/benchmark_original_jarmusch_urine_t0_nested.py\` and dedicated CI \`.github/workflows/jarmusch-original-urine-vs-plasma-day9-pk.yml\`.

## Same ten patients, same observed day9 CL/F, same nested folds

MAE original CL/F numeric source units; original physical unit not independently verified. Lower better. Outer leave-one-original-person-out; inner 4fold train-only k=0/1/3 and Ridge alpha1/20/100. Nested selection criterion log CL/F error, report outer original numeric CL/F MAE. **No independent validation, no confidence adjustment for multiple models**.

| Predictor arm | J1 plasma | J1 urine | J9 plasma | J9 urine |
|---|---:|---:|---:|---:|
| Weight+sex+PGx | 3.356567 | 3.356567 | 3.356567 | 3.356567 |
| +18 GNPS **putative** steroid/bile ions | 3.535090 | 3.794865 | 3.611280 | 3.976507 |
| +207 GNPS high-quality **putative** MS2 ions | 3.991813 | 3.395145 | 3.450288 | 3.642981 |
| +6349 untargeted MS1 ions | **3.081850** | 3.385119 | 3.705224 | **3.227443** |

### Exactly how many molecular features were actually accepted?

Numbers below describe *independently trained outer fold models*, not personalized feature allocation to the heldout patient.

| Matrix and molecular candidates | k=0 | k=1 | k=3 | Difference original CL/F MAE vs weight+sex+PGx | Paired patient residual bootstrap 95%, descriptive |
|---|---:|---:|---:|---:|---|
| J1 urine +18 source-putative steroid/bile | 5/10 | 1/10 | 4/10 | +0.438298 | [+0.052917,+0.866734] |
| J1 urine +207 MS2 annotated | 9/10 | 1/10 | 0 | +0.038578 | [0,+0.115734] |
| J1 urine +6349 untargeted | 8/10 | 1/10 | 1/10 | +0.028552 | [0,+0.076424] |
| J9 urine +18 putative steroid/bile | 8/10 | 2/10 | 0 | +0.619940 | [0,+1.480244] |
| J9 urine +207 MS2 annotated | 5/10 | 3/10 | 2/10 | +0.286414 | [−0.024040,+0.841742] |
| J9 urine +6349 untargeted | 1/10 | 1/10 | **8/10** | **−0.129124** (3.85% lower MAE) | [−0.470548,+0.240795] |
| J1 plasma +6349 untargeted | 6/10 | 2/10 | 2/10 | −0.274717 (8.19% lower MAE) | [−0.790684,+0.141669] |

**Result:** urine J9 untargeted exhibits a small *numerical* incremental prediction gain and often k3; it is **not robust** at n10 and depends on a high-dimensional candidate set with uncertain molecular identity. Its descriptive bootstrap interval includes zero. Urine J1 and the hypothesis-like putative steroid/bile GNPS panel do not improve held-out PK. Do not confuse "putatively annotated urinary ions" with orthogonally validated **6β-hydroxycortisol/cortisol** or other true endogenous CYP3A4 functional ratios.

## Why a mechanism-first assay may still be worth testing

The original source tested **untargeted MS1** and GNPS MS2 spectral matches, not clinically validated targeted analytes. No actual source-specific quantified urinary cortisol/6β-hydroxycortisol pair with retention time/internal standards/urinary creatinine/collection interval is established in this benchmark. A negative arbitrary molecular panel therefore **does not falsify a specific CYP3A functional substrate-product hypothesis**. Candidate selection must be stage-aware:
- CYP3A5*3 is splice-loss (not missense) and motivates liver/intestinal isoform RNA and enzyme/protein function, not a blanket decision to skip transcriptomics.
- Urinary 6β-hydroxycortisol/cortisol is principally CYP3A4 (not CYP3A5-exclusive), impacted by renal excretion, dilution and circadian rhythm; a spot MS1 ion is not a validated phenotype. Literature: https://pmc.ncbi.nlm.nih.gov/articles/PMC13688978/, https://pubmed.ncbi.nlm.nih.gov/10805065/.
- For tacrolimus, the candidate mechanistic network includes hepatic and enterocyte CYP3A4/5 -> 13-O-demethyltacrolimus, while transported/gut/renal/hematocrit functions affect blood exposure and are not proxies for CYP3A5-specific catalysis: https://pubmed.ncbi.nlm.nih.gov/15951320/.

**Next genuinely new question:** hypothesis-led targeted LC-MS/MS urine metabolite ratios at fixed morning predose times, plasma endogenous ratio, matched genotype and independently observed systemic PK (ideally including target drug tacrolimus); prospective/frozen panel then independent validation. This source is already multiply inspected; further optimizations on n10 should be clearly labelled exploratory.
