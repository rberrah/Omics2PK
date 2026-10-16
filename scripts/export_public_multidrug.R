#!/usr/bin/env Rscript
# Export one pre-registered real NCI-60 pharmacological endpoint, no patients.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=2L) stop("Usage: export_public_multidrug.R <mixOmics checkout> <out_dir>")
repo <- normalizePath(args[[1]],mustWork=TRUE)
output <- args[[2]]
if (!dir.exists(output)) dir.create(output, recursive=TRUE)
load(file.path(repo,"data","multidrug.rda"))
x <- multidrug
stopifnot(nrow(x$ABC.trans)==60L, ncol(x$ABC.trans)==48L,
          nrow(x$compound)==60L, ncol(x$compound)==1429L,
          length(x$cell.line$Sample)==60L)
# This decision is FROZEN before the first predictive evaluation.
drug <- "Doxorubicin"
matches <- which(x$comp.name == drug)
if (length(matches)!=1L) stop("Pre-registered Doxorubicin not uniquely found")
index <- matches[[1]]
y <- as.numeric(x$compound[, index])
if (sum(is.finite(y)) != 60L) stop("Expected 60 observed measured GI50 scores")
if (!all(is.finite(x$ABC.trans) | is.na(x$ABC.trans)))
  stop("Unexpected infinite expression measurement")
ids <- as.character(x$cell.line$Sample)
if (anyDuplicated(ids)) stop("Duplicate biological cell line")
# Pre-registered dataset paired rows by row index in the upstream object.
if (nrow(x$ABC.trans)!=length(ids) || nrow(x$compound)!=length(ids))
  stop("Expression/outcome sample count mismatch")
g <- colnames(x$ABC.trans)
if (length(g)!=48L || anyDuplicated(g) || any(is.na(g)))
  stop("Unexpected transporter names")
raw <- data.frame(subject_id=ids, cancer_class=as.character(x$cell.line$Class),
                  x$ABC.trans, check.names=FALSE)
write.csv(raw,file.path(output,"expression.csv"),row.names=FALSE,na="NA",quote=TRUE)
write.csv(data.frame(subject_id=ids, y=y),
          file.path(output,"measured_doxorubicin_response.csv"),
          row.names=FALSE,quote=TRUE)
cat(sprintf("FROZEN target: %s, compound index: %d, n=%d, p=%d, outcome_sd=%0.5f\n",
            drug,index,length(ids),length(g),sd(y)))
