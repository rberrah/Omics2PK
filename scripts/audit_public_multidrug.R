#!/usr/bin/env Rscript
# Public, real measured pharmacology; NO clinical patient or PK data.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Expected upstream mixOmics checkout directory")
repo <- normalizePath(args[[1L]], mustWork = TRUE)
file <- file.path(repo, "data", "multidrug.rda")
if (!file.exists(file)) stop("Pinned original multidrug.rda missing")
load(file)
if (!exists("multidrug")) stop("Unexpected official dataset R object")
x <- multidrug
stopifnot(is.matrix(x$ABC.trans), is.matrix(x$compound),
          nrow(x$ABC.trans) == 60L, ncol(x$ABC.trans) == 48L,
          nrow(x$compound) == 60L, ncol(x$compound) == 1429L,
          length(x$comp.name) == ncol(x$compound),
          length(x$cell.line$Sample) == 60L,
          !anyDuplicated(x$cell.line$Sample))
wanted <- grep("doxorubicin|adriamycin|paclitaxel|taxol|vinblastine|cisplatin|etoposide", x$comp.name,
               ignore.case=TRUE)
cat(sprintf("measured cells=%d, ABC transporters=%d, compounds=%d\n",
            nrow(x$ABC.trans), ncol(x$ABC.trans), ncol(x$compound)))
if (length(wanted)) {
  for (i in wanted) {
    y <- x$compound[,i]
    cat(sprintf("candidate[%d]=%s; finite=%d; sd=%0.5f\n",
                i, x$comp.name[i], sum(is.finite(y)),
                sd(y[is.finite(y)])))
  }
} else {
  cat("No named candidate matches found; inspect pre-registered compound IDs.\n")
  cat("First 24 names:", paste(head(x$comp.name, 24), collapse = " | "), "\n")
}
