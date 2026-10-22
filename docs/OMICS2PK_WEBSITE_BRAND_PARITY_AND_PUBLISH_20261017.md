# Omics2PK research website — PMx Explain / MultiOmics visual parity

**Implementation:** `site/build_site.py` generates a static research portal independently of the Python pharmacometric engine. The shared top-level navigation is **Apprendre / Analyze / Écosystème / Ressources**. Under **Apprendre**, **Objectif** opens the scientific “why” page and **Démo guidée** opens the step-by-step methodology “how” page. **Analyze** remains explicitly in development. The public UI cannot run SWEET, PK modeling, personal data uploads or patient dosing.

## Design and links

- Four accessible drop-down menus: **Learn**, **Analyze**, **Ecosystem**, **Ressources** (English Resources).
- PMx editorial paper (#fbfaf8), ink, terracotta (#a4441f), light and dark palette, Source Serif 4 / Inter / JetBrains Mono font stacks, fixed 59px header, 1100px mobile breakpoint and purpose-based native `details` dropdowns.
- FR/EN toggle, dark/light theme, mobile hamburger menu, outside click, Escape focus behavior, reduced-motion and original site footer.
- Omics2PK entry is present in its own Ecosystem page alongside PMx Explain and MultiOmics. Other dedicated repos must link here only after the URL has passed public deployment verification.
- Evidence is organized pedagogically under **Apprendre → Objectif**: mechanistic plausibility (EPA100), functional phenotype versus available PGx (Jarmusch), then the authentic metformin n32 predose metabolomics → measured clinical AUC/CLr benchmark.
- **Analyze is intentionally an in-development page**: no data upload, no public model execution and no dosing workflow are exposed yet.

## Offline preview and CI

Run from the root of Omics2PK:

    python site/build_site.py --out site/dist
    python -m unittest discover -s tests -p 'test_omics2pk_pages.py' -v

The pull request CI builds and uploads the artifact `omics2pk-static-pages-preview`, which is a preview artifact, **not proof of live deployment**. Canonical routes are the homepage plus `demo/`, `guide/`, `analyze/`, `ecosystem/`, `disciplines/` and `resources/`. The former `learn/` route is kept only as a compatibility alias for the guided course and is excluded from the sitemap.

## First deployment of the standalone repository

1. In GitHub `rberrah/Omics2PK` open **Settings → Pages → Build and deployment → Source: GitHub Actions**.
2. Open **Actions → Publish Omics2PK website to GitHub Pages (manual)** and select **Run workflow**, branch **main**. The build validates all pages and uses official `actions/configure-pages@v5`, `actions/upload-pages-artifact@v4` and `actions/deploy-pages@v4`.
3. If Pages is enabled, the expected project path is **https://rberrah.github.io/Omics2PK/**. Verify live HTTP 200 for the homepage and canonical `demo/`, `guide/`, `analyze/`, `ecosystem/` and `resources/` routes before announcing availability; the presence of a workflow in GitHub is not verification of a deployed website.
4. On live verification, add its exact URL to **Ecosystem** in PMx Explain and MultiOmics. Do not change the codebases owned by their separate maintainers unintentionally.

## Scientific status

- SMX 26 original hepatocyte donors have aggregate source-verified exploratory LODO results; the same 15 donor comparison shows activities versus proteins and genotype but NO statistically supported uplift. Numeric original report `docs/SMX_NAT_26_ORIGINAL_MEASURED_LODO_RESULTS_20261017.md`. Original publisher donor rows are not republished.
- Metformin n32 is now a genuine original-source human **Omics→clinical PK** benchmark with predose metabolomics and measured AUClast/CLr. Its prespecified incremental comparison does **not** confirm added value of omics over the clinical reference. See `docs/E3_REAL_METFORMIN33_PREDOSE_384_OMICS_AUC_CLR_HELDOUT_20261018.md`.
- This does not constitute external validation, Omics→PopPK validation, patient-specific dosing validation or a clinical recommendation.
