# Experimental SWEET-inspired inductive networks — methodological gate

This method is a **separate exploratory extension**. It is **not the original
SWEET algorithm**, has not been scientifically validated against patient PK/PD,
and should not be called a CYP/ABCB1 activity measurement.

## What is implemented

- `fit_frozen_reference`: correlations and query-weight calibration on
  **training subjects only**; gene order, tissue and reference IDs frozen.
- `predict_network`: a new subject modifies only their own graph, without
  updating the frozen training reference.
- `crossfit_training_features(folds=0)`: leave one training subject out
  when building that subject's network. `folds>=2` holds out entire
  deterministic, subject-ID-based groups; fold/reference size is recorded.
- `build_outer_split_features`: explicitly separates training and test
  subjects; uses cross-fitted features for training and the full **outer
  training reference** for each unseen test subject. This function has no
  endpoint parameter and will not accept duplicate or overlapping subject IDs.
- The outputs are per-gene unsigned connection strengths, unvalidated
  descriptions of the network, not mechanistic covariates.

Two constructions are retained.

**v0 (historical Omics2PK prototype):**
`A = corr_genes(X_train)`, `B = corr_genes([X_train, x])`, then
`G0(x) = A + w0(x) * (B-A)`.

**v1 (recommended candidate for new experiments, query-anchored):**
`A_x = corr_genes([X_train, x])` and
`B_x = corr_genes([X_train, x, x])`, then
`G1(x) = A_x + w1(x) * (B_x-A_x)`.

v1 mirrors the perturbation structure of published SWEET more closely: the
target sample exists once in its query-specific anchor and is duplicated once
in the perturbed matrix. The **training reference remains frozen**: no other
held-out subject, PK endpoint or outcome is used. The query's own expression is
allowed to construct its own network, just as a deployed model may use the new
patient's measured RNA.

v0 remains the **API compatibility default** so historical Omics2PK scripts
do not silently change numerical results. New qualification workflows must
request v1 explicitly with `construction="query_anchored_duplicate_v1"`.

### Sample weighting and difference from official SWEET

The published method computes sample weights using the **full expression
universe**, independently from the smaller subset of genes selected for
individual network edges. This prototype now accepts separate
`weight_expression` / `weight_genes` arrays, aligned to the **same
training subjects**, and `weight_query` for each unseen subject.
If these arguments are omitted, the candidate uses the selected gene panel
for both weights and edges. This difference is explicit in the outputs.

For a new, truly unseen patient, weight normalization still uses only the
training subject-similarity distribution and clips values outside its range;
this is **not** the transductive published sample-weight procedure. The
reference sizes also differ by one between LOO training and held-out tests.

**Measured public-data diagnostic (updated 2026-10-22):**
the original SWEET example has 80 expression samples and 100 genes.
Using all 100 genes for sample weighting and a fixed first-25-gene edge panel,
the first 60 samples form the frozen training cohort and the last 20 are
held out. Training features are cross-fitted; every test subject uses the same
60-person reference fingerprint.

Against the published SWEET **transductive** 80-person formula, relative
Frobenius differences across the 20 held-out subjects are:

- historical v0: mean **0.27257**, median **0.27380**;
- query-anchored v1: mean **0.26729**, median **0.27115**;
- v1 is closer for **19/20** held-out subjects;
- first held-out subject: **0.27065 → 0.25586**.

This is a modest mathematical improvement, **not equivalence**. The official
comparison includes all 80 subjects in its reference and is therefore
transductive; a deployable inductive method cannot be expected to reproduce it
exactly without reintroducing held-out cohort information. No outcome labels,
PK/PD concentrations or biological claims enter this diagnostic.

This remains **SWEET-inspired**: not a clinical predictor and not a
validated replacement for official SWEET. Do not present correlation
network strengths as mechanistic enzyme activities.


### Real human-liver expression stress test — GSE9588

A second diagnostic uses the authentic processed **GSE9588** human-liver
microarray matrix: 427 GEO sample columns, first 75 as a fixed reference and
352 as unseen query columns. The same 24 deterministic, training-variable but
**unannotated** array probes are used for v0 and v1.

Observed on the official GEO matrix (2026-10-22):

- v0: **238/352 unique graph hashes**, 114 collisions;
- v1: **352/352 unique graph hashes**, 0 collisions;
- both: **117 clamped weights** and **115 zero weights**;
- clamp direction: **115 low**, **2 high**;
- mean query similarity z-score versus the 75-reference distribution:
  **−1.757**, range **−5.476 to +2.382**; **69/352** queries have |z| > 3;
- unclipped weight position: median **0.156**, range **−0.677 to 1.052**.

Thus v1 fixes a specific software-personalization failure of v0: when the
network-update weight is zero, v1 still retains a query-specific anchor rather
than collapsing every affected subject to the identical training-only graph.
It does **not** solve weight calibration or source distribution shift. The
strongly negative similarity shift also makes the historical choice of the
**first 75 GEO columns** a methodological sensitivity factor: GEO source order
was not designed as a representative training sample. A separate outcome-free,
order-independent reference-selection sensitivity is therefore required before
interpreting the clamp rate as an intrinsic property of v1. There is still no
verified GSM→independent-donor crosswalk, no probe→gene mechanistic annotation
in this smoke test and no matched CYP activity/PK endpoint.
Therefore this result is not evidence of predictive PK utility.

### Sample-weight universe sensitivity on the same GSE9588 split

The historical GSE9588 smoke used the same 24 probes both for network edges and
for sample-weight similarity. This is narrower than the published SWEET design,
where sample weights are derived from a much broader expression universe.
A source-identical sensitivity therefore kept the 24 edge probes fixed but
expanded sample weighting to **512 finite probes selected for non-zero variance
inside the same 75-reference training set only**.

On the same 75-reference / 352-query split, query-anchored v1 changed from:

- clamped weights: **117 → 24**;
- low clamps: **115 → 22**; high clamps remain **2 → 2**;
- zero weights: **115 → 22**;
- queries with |similarity z| > 3: **69 → 17**;
- mean similarity z: **−1.757 → −1.301**;
- median unclipped weight position: **0.156 → 0.333**;
- minimum unclipped position: **−0.677 → −0.230**;
- v1 graph hashes remain **352/352 unique** in both settings.

For legacy v0, broader weighting also improves uniqueness from **238/352 to
331/352** graphs because fewer queries collapse to a zero-update reference.

This strongly indicates that the narrow 24-probe sample-weight universe was a
major source of the apparent calibration failure. It does **not** prove that
512 probes are sufficient or equivalent to a genome-wide SWEET input.

A second outcome-free sensitivity replaced the source-order `first75`
reference with 75 GSM IDs selected deterministically by SHA-256 rank. Each
reference strategy independently selected its 512 train-variable weight probes;
neither expression values nor outcomes influenced which GSMs entered the
reference. With v1:

- clamps: **24 → 11**;
- low clamps: **22 → 5**; high clamps: **2 → 6**;
- zero weights: **22 → 5**;
- |similarity z| > 3: **17 → 1**;
- mean similarity z: **−1.301 → −0.200**;
- median unclipped position: **0.333 → 0.545**;
- v1 graph hashes remain **352/352 unique**.

This demonstrates that GEO source order was another major contributor to the
apparent distribution shift. It does **not** establish SHA-256 selection as an
optimal biological reference; it is deliberately value/outcome-independent and
serves only as an order-bias diagnostic.

A nested SHA-256 reference-size sensitivity then used n=50, 75, 100 and 150,
with 512 train-variable sample-weight probes selected independently inside each
reference:

| Reference n | Query n | v1 clamps | Clamp fraction | Low / high | Zero weights | |z|>3 | Mean z |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 50 | 377 | 17 | 4.51% | 10 / 7 | 10 | 0 | −0.097 |
| 75 | 352 | 11 | 3.13% | 5 / 6 | 5 | 1 | −0.200 |
| 100 | 327 | 12 | 3.67% | 6 / 6 | 6 | 2 | −0.161 |
| 150 | 277 | **0** | **0%** | 0 / 0 | 0 | 0 | −0.201 |

All four v1 runs retain **100% unique graph hashes**. The result is not
monotonic at small/medium n, so it does not justify a universal threshold.
The original zero-clamp n=150 result was therefore tested across four
additional, predeclared salted SHA-256 reference rankings. Reference membership
used GSM identifiers only; expression values and outcomes never selected the
150 reference samples.

Across the **five** n=150 references (277 held-out queries each):

| Reference | v1 clamps | Clamp fraction | Low / high | Zero weights | |z|>3 | Mean z |
|---|---:|---:|---:|---:|---:|---:|
| historical sha256 baseline | 0 | 0% | 0 / 0 | 0 | 0 | −0.201 |
| repeat 1 | 3 | 1.083% | 3 / 0 | 3 | 0 | −0.060 |
| repeat 2 | 1 | 0.361% | 0 / 1 | 0 | 0 | +0.151 |
| repeat 3 | 5 | 1.805% | 3 / 2 | 3 | 1 | −0.007 |
| repeat 4 | 0 | 0% | 0 / 0 | 0 | 1 | −0.137 |

Thus **2/5** deterministic n=150 references have zero clamps and the worst
observed clamp fraction is **1.805%**. n=150 is materially better calibrated
than the historical first-75 reference, but the single zero-clamp result was
not universally reproducible. This does **not** justify a universal
"minimum n=150" rule. Reference composition still matters.

The next gate is therefore no longer to search for a convenient reference
size. It is to qualify upstream train-only preprocessing/normalisation and then
test stability across biologically independent donor/site cohorts when such
source identities and outcomes are available.

## Gates and known limitations

1. **Public outcome-free sanity check:** official SWEET 80-sample RNA
   example, sample weights from all 100 genes and edges from a
   predeclared first-25-gene panel; first 60 samples as reference
   training, last 20 as held-out. Leave-one-out among training subjects,
   training-frozen reference for every held-out subject, and a test for
   invariance of one held-out result to perturbing another held-out profile.
   This **cannot** establish an improvement in PK/PD prediction.
2. **No outer-fold normalization leak:** this function receives a
   preprocessed matrix. Upstream feature selection, batch correction,
   scaling, imputation and gene mapping must be fit **inside each outer
   training fold** and frozen for test subjects. Not currently implemented.
3. **Training/test distribution mismatch and weight calibration:** for v1, a
   cross-fitted training query has an anchor reconstructed from its n-1
   reference + itself, whereas an outer test query uses the full n-person
   training reference + itself. Real GSE9588 removes graph collisions but still
   yields 117 clamped / 115 zero weights, so calibration remains an explicit
   unresolved gate.
4. **Biological scope:** experimental 128-gene maximum; heterogeneous
   RNA cell types and tissue availability may invalidate transfer.
5. **External scientific comparison:** compare original transductive SWEET,
   v0/v1 frozen-reference scores, raw expression, simple PCA/pathway scores, null
   and permuted networks on a **measured molecular or pharmacological target**.
   Avoid post-hoc claims of causality or CYP functional activity.
6. **Graph-to-PopPK:** do not build M3/M4 until an appraised, patient-paired
   PK and omics dataset exists and a real NLME model can be fitted.

## Local execution

```bash
python -m pip install -e '.[sweet-inductive]'
python -m unittest discover -s tests -p test_crossfit.py -v
python -m unittest discover -s tests -p test_inductive_sweet.py -v
python scripts/pilot_inductive_public.py --upstream /path/to/SWEET \
  --report inductive_pilot.json
```

Public CI downloads a pinned upstream SWEET example **only**; individual
study RNA and any PK/clinical data remain on authorized local machines.
