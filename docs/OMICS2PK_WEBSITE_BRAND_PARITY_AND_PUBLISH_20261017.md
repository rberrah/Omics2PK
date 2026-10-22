# Omics2PK research website — PMx Explain / MultiOmics visual parity

**Implementation:** `site/build_site.py` generates a static research portal independently of the Python pharmacometric engine. The public information architecture now separates **scientific demo / why**, **guided methodology course / how**, and **Analyze / in development**. The public UI cannot run SWEET, PK modeling, personal data uploads or patient dosing.

## Design and links

- Four accessible drop-down menus: **Learn**, **Analyze**, **Ecosystem**, **Ressources** (English Resources).
- PMx editorial paper (#fbfaf8), ink, terracotta (#a4441f), light and dark palette, Source Serif 4 / Inter / JetBrains Mono font stacks, fixed 59px header, 1100px mobile breakpoint and purpose-based native `details` dropdowns.
- FR/EN toggle, dark/light theme, mobile hamburger menu, outside click, Escape focus behavior, reduced-motion and original site footer.
- Omics2PK entry is present in its own Ecosystem page alongside PMx Explain and MultiOmics. Other dedicated repos must link here only after the URL has passed public deployment verification.
- Evidence cards explicitly present prior **real** original donor-level study analyses in aggregate without donor values: EPA100 CYP, TKI21, SMX26. Clinical PK is still unverified.
- Analyze is a static research-result browser and does not pretend to be a hosted omics analysis or dosing service.

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
- Metformin n17 same-patient 255 original concentrations and pre-dose/peak plasma metabolomics is a **priority E3 lead** only. The official PMC OA article metadata is accessible under CC BY but supplies no verified original person×omics intensity matrix to the audit. See `docs/METFORMIN17_CLINICAL_OMICS_PK_OPEN_SOURCE_AUDIT_20261017.md`.
- No clinical omics→PopPK patient-specific prediction or dose recommendation is approved.
