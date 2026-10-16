#!/usr/bin/env Rscript
# Cross-study, research-only NCI-60 lineage comparison. NOT same-specimen omics.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3L) stop("Usage: export_public_nci60_twoomics.R <mixOmics> <IntLIM NCI60> <output>")
mixdir <- normalizePath(args[[1]], mustWork=TRUE)
intdir <- normalizePath(args[[2]], mustWork=TRUE)
out <- args[[3]]
dir.create(out, recursive=TRUE, showWarnings=FALSE)

load(file.path(mixdir, "data", "multidrug.rda"))
source_ids <- as.character(multidrug$cell.line$Sample)
source_class <- as.character(multidrug$cell.line$Class)
if (length(source_ids)!=60L || anyDuplicated(source_ids)) stop("Invalid source drug IDs")
drug <- which(multidrug$comp.name=="Doxorubicin")
if (length(drug)!=1L) stop("Pre-registered drug is not unique")
y <- as.numeric(multidrug$compound[,drug])
if (length(y)!=60L || any(!is.finite(y))) stop("Measured response incomplete")
meta <- read.csv(file.path(intdir,"pData.csv"), check.names=FALSE, stringsAsFactors=FALSE)
if (!all(c("cell_line","cancertype") %in% names(meta))) stop("Missing source metadata")
rna <- read.csv(file.path(intdir,"geneData.csv"), check.names=FALSE, stringsAsFactors=FALSE)
met <- read.csv(file.path(intdir,"metabData.csv"), check.names=FALSE, stringsAsFactors=FALSE)
if (ncol(rna)<50L || ncol(met)<50L) stop("Missing source data columns")

# This canonicalization handles typography only: NO biological aliases.
key <- function(x) gsub("[^A-Z0-9]", "", toupper(x))
safe_keys <- function(x, label) {
  k <- key(x)
  if (any(!nzchar(k)) || anyDuplicated(k)) stop(paste("Ambiguous",label,"keys"))
  k
}
drug_keys <- safe_keys(source_ids,"drug")
meta_keys <- safe_keys(meta$cell_line,"metadata")
rna_keys <- safe_keys(names(rna)[-1],"RNA")
met_keys <- safe_keys(names(met)[-1],"metabolites")
matched_keys <- drug_keys[drug_keys %in% meta_keys & drug_keys %in% rna_keys &
                          drug_keys %in% met_keys]
if (length(matched_keys)!=48L) stop(sprintf("Drifted cross-study intersection: %d, expected 48",length(matched_keys)))
ds <- match(matched_keys,drug_keys)
ms <- match(matched_keys,meta_keys)
rs <- match(matched_keys,rna_keys)+1L
ts <- match(matched_keys,met_keys)+1L

# Independently collected lineage labels must not disagree. No data-driven reassignment.
tumor_key <- function(x) {
  z <- key(x)
  z[z %in% c("LUNG","NONSMALLCELLLUNG","NONSMALLCELLLUNGCANCER")] <- "NSCLC"
  z[z %in% c("LEUK","LEUKAEMIA")] <- "LEUKEMIA"
  z[z=="MELAN"] <- "MELANOMA"
  z[z=="OVAR"] <- "OVARIAN"
  z
}
class_mix <- tumor_key(source_class[ds])
class_intlim <- tumor_key(as.character(meta$cancertype[ms]))
if (any(!nzchar(class_mix)) || any(!nzchar(class_intlim)))
  stop("Unknown tumor lineage category")
discordant <- which(class_mix != class_intlim)
if (length(discordant)!=1L || source_ids[ds][discordant]!="MDAMB435")
  stop(paste("Unexpected cross-study lineage disagreement:",
             paste(source_ids[ds][discordant],collapse=", ")))
# Prespecified SOURCE-CONCORDANT cohort: exclude the historic MDA-MB-435
# source=BREAST versus IntLIM=Melanoma, never silently relabel either.
cat(sprintf("EXCLUDED %s: mixOmics=%s vs IntLIM=%s; no biological alias inference\n",
            source_ids[ds][discordant],source_class[ds][discordant],
            meta$cancertype[ms][discordant]))
keep <- setdiff(seq_along(ds),discordant)
ds <- ds[keep]; ms <- ms[keep]; rs <- rs[keep]; ts <- ts[keep]
class_mix <- class_mix[keep]
if (length(ds)!=47L) stop("Expected 47 source-lineage-concordant lines")

# Separate public, curated log-scale source measurements; do not log-transform again.
write_modality <- function(tbl, col_ids, file) {
  feats <- as.character(tbl[[1]])
  if (anyNA(feats) || any(!nzchar(feats)) || anyDuplicated(feats))
    stop(paste("Duplicate or missing feature IDs:",file))
  values <- as.matrix(tbl[,col_ids,drop=FALSE])
  suppressWarnings(storage.mode(values) <- "double")
  if (any(is.infinite(values)) || all(is.na(values))) stop(paste("Invalid",file))
  samples_by_features <- t(values)
  colnames(samples_by_features) <- feats
  out_tbl <- data.frame(subject_id=source_ids[ds], samples_by_features,check.names=FALSE)
  write.csv(out_tbl,file.path(out,file),row.names=FALSE,quote=TRUE,na="NA")
  c(features=length(feats), missing=sum(is.na(values)))
}
rna_info <- write_modality(rna,rs,"rna.csv")
met_info <- write_modality(met,ts,"metabolite.csv")
write.csv(data.frame(subject_id=source_ids[ds],cancer_class=class_mix,
                     y=y[ds]),file.path(out,"response.csv"),row.names=FALSE)
# Local temporary provenance, not uploaded as workflow artifact.
write.csv(data.frame(drug_cell_line=source_ids[ds],
                     intlim_metadata=meta$cell_line[ms],
                     intlim_rna=names(rna)[rs],
                     intlim_metabolite=names(met)[ts],
                     tumor_class_source=source_class[ds],
                     tumor_class_intlim=meta$cancertype[ms],
                     match_method="punctuation_case_only"),
          file.path(out,"provenance_crosswalk.csv"),row.names=FALSE)
cat(sprintf("EXPLORATORY CROSS-STUDY CELL LINES=%d, RNA=%d features, metabolites=%d features; missing RNA=%d, metabolites=%d\n",
            length(ds),rna_info[["features"]],met_info[["features"]],
            rna_info[["missing"]],met_info[["missing"]]))
cat("NOT SAME BIOLOGICAL SPECIMENS; NO CLINICAL PK ENDPOINT; 12 SYNTAX-UNMATCHED + 1 HISTOLOGY-DISCORDANT LINE EXCLUDED\n")
