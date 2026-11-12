# Reproducible report

`omics2pk_results.Rmd` renders the executed objectives — H1 on ARIANA and H5 on
Müller 2026 — into manuscript tables and figures.

**Every number in the report is read from JSON**, never retyped, so the report
cannot drift from the computation that produced it. The JSON in `data/` is
committed; the rendered HTML is not, because it is 1.4 MB and regenerable.

## Rendering

```r
# from the repository root
source("reports/render.R")
```

or directly:

```bash
Rscript -e "rmarkdown::render('reports/omics2pk_results.Rmd')"
```

If pandoc is not on `PATH`, point R at the copy RStudio ships. On Windows that
is usually:

```r
Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")
```

`render.R` does this automatically when it finds it.

Requires R with `rmarkdown`, `knitr`, `jsonlite`, `ggplot2` and `dplyr`.
Verified with R 4.4.1, rmarkdown 2.32, knitr 1.52.

## Regenerating the JSON

The report's inputs come from the analysis scripts, which verify the source
checksum and exit non-zero on mismatch. Neither source workbook is committed —
acquire them first (see `registry/PUBLIC_DATA_EVIDENCE_REGISTRY.json` for the
URLs and checksums), then:

```bash
ARIANA="path/to/Solanidine_metoprolol_ratios_PGx.xlsx"
MULLER="path/to/Solanidine standardized potato meal.xlsx"

python scripts/run_h1_ariana_incremental_omic_value.py --source "$ARIANA" \
    > reports/data/h1_ariana_trough.json

python scripts/run_h1_ariana_incremental_omic_value.py --source "$ARIANA" \
    --target "c(metoprolol)weighted [ng/ml]" \
    > reports/data/h1_ariana_dose_weighted.json

python scripts/run_h5_solanidine_marker_robustness.py \
    --muller "$MULLER" --ariana "$ARIANA" \
    > reports/data/h5_muller_robustness.json
```

The JSON carries aggregate statistics only. A test asserts no numeric array of
length ≥ 10 appears in it, so a per-patient vector cannot be committed by
accident.

## What the report is allowed to claim

Both objectives are capped at **exploratory**: ARIANA has been analysed
repeatedly by this project, and Müller administers no drug. The report prints
the prohibited conclusions alongside the results rather than only the results,
and reports `ABSTAIN` where a quantity could not be estimated instead of
filling it in.
