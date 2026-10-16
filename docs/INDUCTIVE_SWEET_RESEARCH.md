# Experimental inductive gene network: separate candidate, not official SWEET

Status: **research hypothesis only, no measured-PK evaluation**.
Upstream SWEET is transductive; its archived edge outputs may be reproduced
directly, but that does **not** establish how an unseen patient should be scored.

## Hypothesis

Let \`X_ref\` be a gene × subject expression matrix for **training individuals
only**. Define the reference gene correlation \`A = corr_gene(X_ref)\`.
For the unseen subject's aligned expression \`x\`, recompute
\`B(x) = corr_gene([X_ref, x])\`. Use the separately trained sample-correlation
distribution to assign a bounded multiplier \`w(x)\`, and form

\`G(x) = A + w(x) × (B(x) - A)\`.

This is **SWEET-inspired, not equivalent to published SWEET**. Specifically,
the proposed new-patient weight uses correlations of \`x\` with reference
subjects and calibrates the range from training-only subject–subject
similarities. It is clipped to a predeclared range to avoid unbounded
extrapolation. The reference and scale are never updated using other held-out
patients or outcomes.

## Research controls already implemented

- Matched gene order, tissue and sample identities are required.
- Query subject must be new; no shared subject can appear in the reference.
- Frozen read-only reference arrays, stable SHA-256 fingerprint, deterministic
  outputs and tests verifying predictions are invariant to other queries.
- Reject constant genes, query expression with zero variance, NaN/Inf,
  insufficient reference samples, degenerate weights, shape/identifier errors.
- Hard limit of 128 genes for this preliminary test: panels must be chosen
  using training data or truly prior evidence, never held-out PK outcomes.
- Returns network and unsigned connection strength descriptors only; those
  are **not** estimated CYP3A enzymatic activities.

## Unresolved before M1bis/M3/M4 evaluation

1. **Training network mismatch:** training subjects must obtain features
   under a comparable leave-one-out or cross-fitted reference protocol;
   there is no validated approach here yet. Using all training subjects in
   the reference for their own feature construction would bias the comparison.
2. Compare with upstream transductive SWEET on its public expression example
   and quantify the effect of reference size, batch, tissue, scale and outliers.
3. Subject-level nested CV must refit expression preprocessing and reference
   for every outer fold; no target-informed feature/gene selection outside folds.
4. Network weight clipping may distort effects for OOD patients; quantify.
5. Only then test the *measured* PK/PD endpoint, vs raw expression, PCA/MOFA2,
   and no-network controls, with and without Pillar A. Do not publish this
   formulation as validated SWEET or clinical utility.

## Local testing

\`\`\`sh
python -m pip install -e '.[sweet-inductive]'
python -m unittest discover -s tests -p test_inductive_sweet.py -v
\`\`\`

No patient data are downloaded or sent over the network.
