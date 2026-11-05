# Latham PHH real public-source reproduction (2026-11-04)

## Qualified original material

- Latham et al., *Drug Metabolism and Disposition* 2025; DOI **10.1016/j.dmd.2025.100196**, PMC **PMC12990092**, printed **Table 4**.
- Source retrieved **live from official PubMed Central HTML** in GitHub Actions run [41805204594](https://github.com/rberrah/Omics2PK/actions/runs/41805204594) (Europe PMC JATS attempted first); official fulltext source SHA-256 from the first successful replay: `14c1256ab782007a156a22f7a4814427e6952fcdad07cf4f989c6494757911a9`. Source HTML hashes may change with page rendering; scripts always report the current acquired-byte hash.
- The source table contains 14 separately identified primary human hepatocyte donors with CYP2C8/CYP3A4/CYP3A5 genotypes, measured CYP2C8 and **shared CYP3A4/5** probe activities, targeted CYP2C8/CYP3A4/CYP3A5 protein concentrations, apparent imatinib depletion clearance and its N-desmethyl formation. **One donor lacks all three CYP3A protein fields.** Four-arm protein comparisons require the same **13** complete biological donors.
- Original CYP3A5 phenotype classes: **NM 1, IM 3, PM 10** (CPIC-derived qualitative groups). *3/*6 and *3/*7 are represented in the original rows. This does NOT validate calling complete allele coverage in another cohort.
- Original Table 4 clearance header: **µL/million cells/min**, consistent with the paper's PHH clearance equation; a prose paragraph separately reports µL/mg protein/min for sex-specific means. Report model values on the **original Table 4 numeric scale** and preserve the unit discrepancy, rather than inventing a conversion.

## Four deliberately fixed original E2 comparators

Implementation: `src/omics2pk/latham_phh_source_gate.py` and `scripts/benchmark_latham_phh_original.py`. Five checks before fitting: correct DOI/table structure, source genotype-class balance, one and only one missing CYP3A protein panel, original donor identity and complete-case common subset, and CYP3A4+CYP3A5 protein total consistency. The original table is held in memory only; neither sample keys nor original rows, predictions or losses are committed or printed.

**Models:** one fixed `Ridge(alpha=1)` on the log apparent *in-vitro* imatinib clearance; one original donor held out at a time, and training-only feature mean/SD. The genotype baseline encodes **CYP2C8 non-*1/*1 carrier flag plus CYP3A5*1 copy count**. This encoding is a simple descriptive comparator, **not a comprehensive CYP2C8 functional allele model**. Do not confuse CYP2C8 enzymatic probe and CYP3A probe with CYP3A5-specific activity.

| Same 13 PHH donors, same heldout folds | Replayed source-table MAE |
|---|---:|
| CYP2C8 + CYP3A5 PGx | 2.906576 |
| PGx + CYP2C8, CYP3A4, CYP3A5 targeted protein abundances | 3.003704 |
| PGx + CYP2C8 and shared CYP3A probe activity | **2.762944** |
| PGx + proteins + both probe activities | 3.262812 |

These are **genuine source-read real human-donor results**, first independently verified in a successful GitHub workflow on 2026-11-04. They are a single tiny exploratory study, not a new external cohort, not clinical PK (E3), not predicted tacrolimus AUC, not a validated effect of curated omic priors, and not a demonstration that proteomics is redundant. No nested model selection or clinically meaningful interval is supported by n13.

## 2026-11-05: independent endpoint-sensitive error metrics and donor-paired uncertainty

**Live original-source replay SUCCESS**:
[GitHub Actions #42015617333](https://github.com/rberrah/Omics2PK/actions/runs/42015617333).
The source reader acquired the official PMC HTML and reproduced the exact
original **13 protein-complete biological donor** rows. The original absolute
source-unit MAEs remain unchanged. Additional *held-donor-out* predictions
were retained **in memory only**, with no donor IDs/predictions released.

| Original same-n13 model arm | Linear MAE ↓ | Linear RMSE ↓ | Mean absolute log(pred/obs) ↓ | Log RMSE ↓ | Within ±20% |
|---|---:|---:|---:|---:|---:|
| PGx only | 2.906576 | 3.845941 | 0.432147 | 0.548539 | **5/13** |
| PGx + targeted CYP2C8/3A4/3A5 proteins | 3.003704 | 4.184066 | **0.408474** | **0.535347** | **5/13** |
| PGx + CYP2C8 and shared CYP3A probe activities | **2.762944** | **3.778989** | 0.415360 | 0.551937 | **5/13** |
| PGx + both | 3.262812 | 4.597902 | 0.445258 | 0.588263 | **5/13** |

The ratio-based P20 criterion (0.8–1.2) is an **accuracy summary of
in-vitro clearance prediction**, not a dose attainment probability or
clinical therapeutic outcome. All four also have exactly **10/13
(76.9%)** observations within a factor of two.

**Paired original donor-residual bootstrap**, 999 fixed-seed draws,
95% percentile interval; same independent donors and original
heldout predictions, **without refitting or correcting post-hoc model
selection**:

| Prespecified arm minus PGx alone | ΔMAE (negative better), 95% interval | Δ mean absolute log error, 95% interval |
|---|---|---|
| + proteins | +0.097128 [−0.970203, +1.602390] | −0.023673 [−0.149653, +0.144656] |
| + probe activities | −0.143632 [−0.953050, +0.709198] | −0.016787 [−0.129581, +0.106302] |
| + both | +0.356235 [−0.764823, +1.862689] | +0.013110 [−0.122345, +0.191671] |

**Scientific conclusion:** neither molecular augmentation shows a
reliably established incremental benefit in n13; the apparent
“best omic” **changes with the chosen loss function**.
The metric should be prespecified according to the intended
drug-exposure and dosing consequence, not chosen post hoc to favor a
modality. Linear MAE weights high-clearance absolute deviations more
than log-scale error, whereas log errors emphasize multiplicative
prediction ratios and relative error. P20 cannot discriminate these
arms here. Neither assay measures CYP3A5-selective metabolism;
**none is human patient plasma PK or validation of a dosing strategy**.

Source and code:
`src/omics2pk/latham_phh_source_gate.py`,
`tests/test_latham_phh_source_gate.py`. Published main text has a
unit discrepancy; do not silently change the source numeric scale.

## Critical biological counterexample: molecular signal is not synonymous with function

**Of the ten CYP3A5 CPIC nonexpresser-class donors**, nine have complete source quantitative CYP3A5-targeted peptide measurements, and **all nine have a positive measured CYP3A5-targeted protein signal**. All **10/10** nonexpresser-class donors have positive **shared CYP3A4/5 midazolam-probe activity**. This is not a paradox to 'correct' by assigning high CYP3A5 function:

- The proteomics assay targets a peptide, not independently assayed complete active CYP3A5 catalytic capacity.
- Shared CYP3A probe activity can originate largely from CYP3A4.
- CYP3A5*3/*6/*7 remain CPIC no-function categories; a low or detectable peptide signal does not imply a universal numeric residual activity.

These counts provide a **real-data stress test** directly supporting the proposed distinction between measured RNA/peptide/metabolite signal and true specific enzyme function. They do **not** establish the molecular source of the peptide or the magnitude of in-vivo CYP3A5 clearance.

## How to reproduce and what is deliberately not claimed

```bash
python -m pip install -e .
python -m pip install 'beautifulsoup4==4.12.3'
python -m unittest discover -s tests -p 'test_latham_phh_source_gate.py' -v
python -m scripts.benchmark_latham_phh_original
```

The source reader first tries official Europe PMC JATS, then official PMC HTML; both are restricted to PMC12990092 and source DOI. SHA-256 of the acquired fulltext is logged, not the source donor rows. The execution returns an error/abstention if the current source is unavailable, modified or structurally unqualified. Synthetic source-parser unit tests do not substitute for a successful source replay.

**Next experimental gate:** new matched donor sample with full *3/*6/*7 and CYP3A4 coverage plus **independent CYP3A5-specific activity**, pre-landmark single/multiomics and then genuine drug PK. The PHH data can calibrate scientific caution but cannot train and validate that full target.
