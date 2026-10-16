# SWEET outputs → tabular M1bis (research POC)

Upstream: [SysMednet/SWEET](https://github.com/SysMednet/SWEET), commit \`cf064eb5f5b996c40bba3cce507384f6c0fdecb5\` (MIT). This repository **does not vendor SWEET source**.

## Offline workflow

1. Export normalized/transformed RNA from MultiOmics (native \`exchange/\` contract).
2. \`omics2pk prepare-sweet exchange/ sweet-input/\` — creates \`expression.txt\`, \`patient.txt\`, \`gene.txt\`. Strictly input preparation, not network inference.
3. Run the original SWEET scripts **locally** in a separately pinned checkout using their own \`-f\`, \`-w\`, \`-p\`, \`-g\` arguments. These have scientific caveats below.
4. \`omics2pk summarize-sweet /path/to/edge-outputs sweet-input/patient.txt sweet-features/\`

Produced outputs: \`network_summary.csv\`, \`gene_absolute_strength.csv\`, \`gene_signed_strength.csv\`. Values are **topological descriptive features**, not biological activity estimates.

**Important:** this reader expects a complete undirected graph (all n(n-1)/2 unique pairs) on the same gene universe in each patient; stops if an edge is missing, duplicated, infinite/NaN, or if sample IDs are unsafe. A 200,000-edge input cap applies to this **proof of concept**. Targeted candidate gene panels should be fixed before fitting PK, ideally on biological grounds or training-only folds. The cap is a software safety guardrail, not evidence that genome-scale SWEET has been validated.

## Scientific blockers for meaningful M1bis / M3 / M4

- **Reference cohort leakage:** the original SWEET sample weights and aggregated correlations are computed using the entire input cohort (transductive). Merely splitting an XGBoost dataset afterwards is not enough. A new unseen-patient reference procedure and its calibration must be demonstrated.
- **Malformed / missing data:** upstream edge script converts missing-like strings to zero, and a constant gene may create non-finite Pearson correlations. Omics2PK rejects NA inputs by default and rejects NaN/Inf edges.
- **Boolean CLI:** upstream \`-z\` uses \`type=bool\` in argparse; avoid passing \`-z False\`, which Python treats as true.
- **Potential size:** a 20,000-gene network contains approximately 200 million edges **per sample**. Avoid indiscriminate full-network fits.
- **No mechanistic semantics:** gene connectivity is statistical association, not CYP activity, metabolism, a genotype, or a directed causal graph.
- **No clinical validation:** synthetic tests only. Reproducing the original published reference outputs numerically and assessing a measured PK/PD dataset are separate future gates.

This module is deliberately useful for initial M1bis comparisons but should not be presented as proof that SWEET improves PopPK.
