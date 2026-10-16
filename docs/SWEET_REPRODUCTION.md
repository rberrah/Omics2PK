# SWEET official numeric reproduction — qualification gate

The first gate runs the original upstream scripts *unmodified*, at commit
\`cf064eb5f5b996c40bba3cce507384f6c0fdecb5\`, on upstream's own publicly
available example. It compares **all published sample weights**, **all edges in
the first published individual network**, its **z-score file**, and the global
mean/std. It verifies the outputs for the source snapshot, not reproducibility
of an independent biological discovery. Tolerance is absolute/relative 1e-8.

Run on an authorized **local** machine with Python 3.11, NumPy 2.3.5 and the
separately cloned pinned upstream repository:

    python scripts/verify_sweet_upstream.py --upstream /path/to/SWEET --report report.json

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

The next separate methodological gate must demonstrate the construction of
a new-patient network from a **training-frozen** reference without using
the held-out individual's PK/PD endpoint.
