# SWEET official numeric reproduction — qualification gate

The first gate runs the original upstream scripts *unmodified*, at commit
\`cf064eb5f5b996c40bba3cce507384f6c0fdecb5\`, on upstream's own publicly
available example. It compares **all published sample weights**, **all edges in
the first published individual network**, its **z-score file**, and the global
mean/std. The archived z-score files are now **forensically explained** rather than merely
flagged as inconsistent. In the upstream git history, the older Step 3 source
used the expression `raw - mean/std`; current upstream code correctly uses
`(raw - mean)/std`. The archived example z-score files were generated with the
older precedence-bug formula and were not regenerated after the code correction.

On both archived patients, all **4,753 edges per patient** match
`raw_edge_score - mean/std` to floating-point precision. With the published
`mean=0.021918999549565974` and `std=0.1922225916111927`, the legacy offset is
approximately `0.114029258298116`. Thus the first archived value
`-0.08907350579171248` is exactly explained by
`0.024955752506403597 - 0.114029258298116`.

The current pinned script still reproduces all sample weights, raw edge scores,
mean/std and the **correct** z-score formula. Archived z-score files are retained
only as historical outputs; they must not be used as validated standardized
SWEET features.
It verifies the outputs for the source snapshot, not reproducibility
of an independent biological discovery. Tolerance is absolute/relative 1e-8.

Run on an authorized **local** machine with Python 3.11, NumPy 2.3.5 and the
separately cloned pinned upstream repository:

    python scripts/verify_sweet_upstream.py --upstream /path/to/SWEET --report report.json
    # Optional strict gate: intentionally errors because the archived files are stale legacy outputs
    python scripts/verify_sweet_upstream.py --upstream /path/to/SWEET --strict-archived-zscores

The GitHub Actions job executes the same **public-only software regression**
as a small automated check. It downloads no controlled data, exports no
individual molecular profiles, and retains only aggregate comparison metrics.

## Important research limitations

- Upstream sample weights and correlations are transductive: they include all
  samples in the input cohort. Splitting downstream XGBoost/PK subjects
  **after** this computation is insufficient evidence of an inductive model.
- The archived example is not an external test set. Numerical reproduction
  cannot establish correct biological interpretation, predictive performance
  or calibration.
- Upstream converts certain missing values to zero; this project defaults to
  refusing unhandled missing expression values.
- Upstream \`-z\` accepts \`type=bool\`, for which passing \`False\` is evaluated
  as true by Python; the runner explicitly passes \`True\`.
- The official edge networks are undirected co-expression. They must remain
  distinct from typed, directed Pillar-A pharmacological relationships.
- Full genome-wide pairwise outputs scale quadratically; this numeric
  comparison intentionally reuses only the published small example.

### Qualification decision

Issue #4 is no longer a normalization-meaning blocker. For Omics2PK, the allowed
SWEET numerical reference is the **current pinned source algorithm** and its
recomputed raw edges / mean / std / parenthesized z-scores. The stale archived
`*_zscore.txt` files are explicitly rejected as reference targets.

This resolves numerical provenance only. The next methodological gate remains
construction of a new-patient network from a **training-frozen** reference
without using held-out expression cohorts or PK/PD endpoints.
