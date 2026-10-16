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

For reference expression `X_train` (genes × subjects),
`A = corr_genes(X_train)`. For one unseen subject's expression `x`,
`B = corr_genes([X_train, x])`. This preliminary variant uses
`G(x) = A + w(x) * (B-A)`, where the weight is a bounded
function of query-to-training expression-profile correlations.

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

**Measured public-data diagnostic (CI, 2026-10-16):**
the original SWEET example has 80 expression samples and 100 genes.
Using all 100 genes for sample weighting and the first 25 genes for
network edges: training 60, held-out 20; training LOO references 59 and
held-out references 60. All test networks share one training-frozen
fingerprint; 1 of 20 test weights was clipped outside the reference range.
For one test subject, the relative Frobenius difference between this
inductive network and the published SWEET **transductive formula**
applied to all 80 subjects was **0.27065** (27.1%), illustrating
non-equivalence under different reference-cohort definitions.
No outcome labels or PK/PD concentrations are present in this diagnostic.

This remains **SWEET-inspired**: not a clinical predictor and not a
validated replacement for official SWEET. Do not present correlation
network strengths as mechanistic enzyme activities.

## Gates and known limitations

1. **Public outcome-free sanity check:** official SWEET 80-sample RNA
   example, sample weights from all 100 genes and edges from a\n   predeclared first-25-gene panel; first 60 samples as reference
   training, last 20 as held-out. Leave-one-out among training subjects,
   training-frozen reference for every held-out subject, and a test for
   invariance of one held-out result to perturbing another held-out profile.
   This **cannot** establish an improvement in PK/PD prediction.
2. **No outer-fold normalization leak:** this function receives a
   preprocessed matrix. Upstream feature selection, batch correction,
   scaling, imputation and gene mapping must be fit **inside each outer
   training fold** and frozen for test subjects. Not currently implemented.
3. **Training/test distribution mismatch:** LOO reference size n-1,
   test reference size n; K-fold references are smaller still. Quantify
   impact on feature scale and calibration.
4. **Biological scope:** experimental 128-gene maximum; heterogeneous
   RNA cell types and tissue availability may invalidate transfer.
5. **External scientific comparison:** compare original transductive SWEET,
   frozen-reference scores, raw expression, simple PCA/pathway scores, null
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
