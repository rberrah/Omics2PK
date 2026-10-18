# Authorized Human Liver Cohort Synapse import: source gate for Omics2PK (18 October 2026)

## What is actually available

Original HLC project [Synapse syn4499](https://www.synapse.org/Synapse:syn4499) has the following **real verified public metadata** (not publicly downloadable original content):

- [syn4629 — curatedPhenotype](https://www.synapse.org/Synapse:syn4629): candidate original biological phenotypes, including potentially CYP enzyme-substrate activities; **content and assay columns NOT verified**.
- [syn88644 — curatedExpressionLiver](https://www.synapse.org/Synapse:syn88644): original liver expression/crosswalk *candidate*.
- [syn89614 — curatedGenotype](https://www.synapse.org/Synapse:syn89614): original original genotype evidence *candidate*.

Metadata is returned for all files via the official Sage Synapse API without authentication, but **anonymous source-content downloads returned HTTP 403 on 18 October 2026**. That could mean account requirement, file permission or data-specific conditions; neither the current rights nor the exact file contents can be determined from HTTP 403 alone. The GEO [GSE9588](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9588) expression matrix is separately and truly **public**, processed GPL4372 microarrays with 39,302 probe rows × 427 GEO GSM columns; no original activity or protein present in GEO processed expression.

**CYP substrate activity values are NOT presently available to our programmatic analysis, and an original GSM→Synapse donor codebook has NOT been verified.** Do not announce a 398-donor clinical/metabolic benchmark until those prerequisites pass.

## User action on synapse.org

1. Create an individual account at [accounts.synapse.org/register1](https://accounts.synapse.org/register1), verify email and configure required two-factor authentication.
2. Open the three Synapse links above when signed in. Look for **Download** or **Request Access**. A basic registered account may be enough to download shared files; other files have controlled access rules, requiring an agreement, certified user account and occasionally validated profile/institutional ethics approval. Follow the actual file-specific permissions; do not assume all original records are public.
3. If access is refused, use the on-page **Request Access** workflow; cite original academic purpose: compare donor-matched liver CYP expression (and genotype/protein if available) with ORIGINAL observed measured CYP-substrate enzyme activity under donor-grouped crossvalidation. Do not claim clinical CE-marked dosing or promise public data redistribution.

Official current Synapse docs: [Accounts](https://docs.synapse.org/synapse-docs/managing-your-account), [Account tiers](https://docs.synapse.org/synapse-docs/synapse-user-account-types), [Conditions for Use](https://docs.synapse.org/synapse-docs/sharing-settings-permissions-and-conditions-for-use), [Python authentication](https://python-docs.synapse.org/en/stable/tutorials/authentication/).

## Local access and original-file audit, now implemented

Never send credentials/PATs through chat or place them in code, notebooks, a GitHub Actions secret, PR or file upload. Authorized **controlled** source data generally must not be uploaded to the public repository or shared with persons without equivalent authorization.

After account registration and actual approved file permissions, locally:

```bash
# In a checked-out Omics2PK repository:
python -m pip install -e .
python -m pip install synapseclient
# Set up your own PRIVATE Synapse auth once via synapse's website Settings
# and either a locally stored .synapseConfig or environment variable
# SYNAPSE_AUTH_TOKEN, never committed to source control.
python -m scripts.download_authorized_synapse_hlc
# For genotype source too (optional; only if authorized):
python -m scripts.download_authorized_synapse_hlc --include-genotype
```

The downloader uses your locally authorized Synapse SDK account, asks Synapse itself for **syn4629 and syn88644** and optionally syn89614, and stores original files **only in `private/hlc_synapse`** which is already gitignored. It does not circumvent any entitlement checks. If denied, it returns role+source-ID+error-class only, never a private response text or token.

The new local-only source auditor `src/omics2pk/hlc_authorized_local_source.py` provides source SHA256, byte length, source format class, and **aggregate number of rows, columns and presence flags for CYP/activity/RNA/protein/donor-key-like headers** where source text is CSV/TSV/gzipped. It never prints original individual identifiers, assay values, file paths or any token. ZIP member count and serialized binary detection are supported; proprietary historical Sage/RData/RDS formats are **not** silently interpreted as missing values. If source is native R or historical SageBioCurated binary, perform a **rights-respecting local conversion** and re-audit; do not invent row joins.

Scripts/tests:
- `scripts/download_authorized_synapse_hlc.py`
- `src/omics2pk/hlc_authorized_local_source.py`
- `tests/test_hlc_authorized_local_source.py` (fabricated-only, includes fake 403/PAT non-disclosure check)
- `.github/workflows/hlc-authorized-local-source-contract.yml` (contract tests only; it cannot download controlled real participant data and does not try)

## Precise original biological source gate — before any real enzyme-activity model

From original `syn4629`, determine original sample/patient IDs, endpoints (e.g. CYP3A4 midazolam hydroxylation; CYP2C9 diclofenac hydroxylation; CYP2D6 dextromethorphan) and **actual units** and missing/unquantifiable values. From syn88644 or original codebook, document a source-proven one-to-one crosswalk of GEO `GSM` ↔ actual biological donor ID. Inspect 427 GEO GSM labels for replicate specimens; neither matching sample counts nor sorted row order proves the biological join.

Only then:
- Same independent original humans, same endpoint, and same held-donor-out folds for age+sex±PGx reference versus RNA-only and adjusted RNA covariates.
- If genotype source available and source-qualified, add genuine PGx reference instead of inventing genotype.
- An original independently observed same-human CYP protein panel is STILL required to answer the **RNA vs RNA+proteins** question. `syn4629` has not been shown to contain any CYP protein abundance; functional CYP probe activities are **phenotypes, not proteomics**.
- Use training-only normalization, probe annotation GPL4372, handling assay LLOQ, source temporal/batch/study provenance, cross-fitted SWEET reference and careful uncertainty. Current SWEET-inspired graph smoke on 24 **unannotated** GEO probes produced 238 distinct graphs / 352 queries, including similarity weight clamping. This is SOFTWARE qualification with duplicated projected graphs, NOT an interpretable enzyme-network feature and must remain nonclinical until corrected.
- No source-level participant measurement or identity should go to public CI, GitHub or a paper appendix if Synapse rights do not allow redistribution.

**Scientific status:** `SYNAPSE_PROJECT_METADATA_VERIFIED; ANONYMOUS_SOURCE_CONTENT_HTTP403; PRIVATE_LOCAL_IMPORT_READY; ACTUAL_DATA_ACCESS_AND_SAME_DONOR_CYP_SOURCE_QUALIFICATION_PENDING`.
