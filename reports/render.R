# Render the Omics2PK results report.
#
# Usage, from the repository root:
#     source("reports/render.R")
# or:
#     Rscript reports/render.R

local({
  here <- if (file.exists("reports/omics2pk_results.Rmd")) "reports" else "."
  rmd <- file.path(here, "omics2pk_results.Rmd")
  if (!file.exists(rmd)) {
    stop("run this from the repository root, or from reports/", call. = FALSE)
  }

  # pandoc: use whatever is on PATH, else the copy RStudio ships.
  if (!nzchar(Sys.getenv("RSTUDIO_PANDOC")) && !nzchar(Sys.which("pandoc"))) {
    candidates <- c(
      "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools",
      "C:/Program Files/RStudio/bin/quarto/bin/tools",
      "/usr/lib/rstudio/bin/quarto/bin/tools",
      "/Applications/RStudio.app/Contents/Resources/app/bin/quarto/bin/tools"
    )
    found <- candidates[file.exists(file.path(
      candidates, if (.Platform$OS.type == "windows") "pandoc.exe" else "pandoc"))]
    if (length(found)) {
      Sys.setenv(RSTUDIO_PANDOC = found[1])
      message("using bundled pandoc: ", found[1])
    } else {
      stop("pandoc not found. Install pandoc or set RSTUDIO_PANDOC.", call. = FALSE)
    }
  }

  needed <- c("rmarkdown", "knitr", "jsonlite", "ggplot2", "dplyr")
  absent <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) {
    stop("missing R packages: ", paste(absent, collapse = ", "), call. = FALSE)
  }

  # The report refuses to render on stale or missing inputs rather than
  # silently reporting whatever happens to be in data/.
  inputs <- file.path(here, "data", c(
    "h1_ariana_trough.json", "h1_ariana_dose_weighted.json",
    "h5_muller_robustness.json"))
  missing <- inputs[!file.exists(inputs)]
  if (length(missing)) {
    stop("missing report inputs:\n  ", paste(missing, collapse = "\n  "),
         "\nSee reports/README.md for how to regenerate them.", call. = FALSE)
  }

  out <- rmarkdown::render(rmd, output_file = "omics2pk_results.html", quiet = FALSE)
  message("\nwrote ", out)
  invisible(out)
})
