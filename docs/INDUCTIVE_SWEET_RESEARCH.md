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

### Crucial difference from published SWEET

Published SWEET computes sample weights using the **entire gene-expression
input**, while selected network edges may use a targeted gene subset. The
current inductive prototype uses **the same selected panel both for weights
and edges**. That is a material methodological difference requiring
independent calibration before comparison or naming equivalence. The weight
normalization also differs for an unseen patient. This variant must **never**
be presented as numerically equivalent to official SWEET.

## Gates and known limitations

1. **Public outcome-free sanity check:** official SWEET 80-sample RNA
   example, predeclared first 25 genes; first 60 samples as reference
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
