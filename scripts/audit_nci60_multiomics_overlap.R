#!/usr/bin/env Rscript
# Inspect real NCI-60 study overlap WITHOUT fuzzy joining, prediction or raw export.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop("Usage: audit_nci60_multiomics_overlap.R <mixOmics checkout> <IntLIM formatted NCI-60 checkout>")
source_mix <- normalizePath(args[[1]], mustWork=TRUE)
source_multi <- normalizePath(args[[2]], mustWork=TRUE)
load(file.path(source_mix,"data","multidrug.rda"))
source_names <- as.character(multidrug$cell.line$Sample)
meta <- read.csv(file.path(source_multi,"pData.csv"), check.names=FALSE, stringsAsFactors=FALSE)
if (!"cell_line" %in% names(meta)) stop("No cell_line identity column in pData")
met_names <- as.character(meta$cell_line)
# Import first header of huge RNA matrix, not full ~18 MB until linkage is confirmed.
rna_names <- colnames(read.csv(file.path(source_multi,"geneData.csv"),
                           check.names=FALSE,nrows=1,stringsAsFactors=FALSE))[-1]
protein_names <- colnames(read.csv(file.path(source_multi,"metabData.csv"),
                               check.names=FALSE,nrows=1,stringsAsFactors=FALSE))[-1]
if (anyDuplicated(source_names) || anyDuplicated(met_names) ||
    anyDuplicated(rna_names) || anyDuplicated(protein_names))
   stop("Duplicate sample names in one source")
if (length(source_names)!=60) stop("mixOmics source has unexpected subjects")
exact_four <- Reduce(intersect,list(source_names,met_names,rna_names,protein_names))
exact_rna <- intersect(source_names,rna_names)
exact_metab <- intersect(source_names,protein_names)
cat(sprintf("Source NCI-60 drug+ABC=%d, second sample metadata=%d, RNA=%d, metabolites=%d\n",
            length(source_names),length(met_names),length(rna_names),length(protein_names)))
cat(sprintf("EXACT overlap drug+RNA=%d, drug+metabolites=%d, all four=%d\n",
            length(exact_rna),length(exact_metab),length(exact_four)))
cat("EXACT matching cell-line labels:",paste(sort(exact_four),collapse=" | "),"\n")
cat("Drug set lacking exact full overlap:",paste(sort(setdiff(source_names,exact_four)),collapse=" | "),"\n")
cat("Other multiomics IDs not in source:",paste(sort(setdiff(union(rna_names,protein_names),source_names)),collapse=" | "),"\n")
if (length(exact_four) < 30L) {
 cat("NOT_QUALIFIED: insufficient exact overlap for multiomic modeling; formal source-by-source alias crosswalk is required\n")
 quit(status=2L)
}
cat("PASS exact cell-line identity overlap audit; NO measured multiomic prediction and no evidence of contemporaneous sample collection\n")
