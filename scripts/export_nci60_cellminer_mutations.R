#!/usr/bin/env Rscript
# Independently sourced, published NCI-60 predicted-deleterious *somatic* gene
# flags (rcellminerData), joined EXACTLY to the already source-qualified n47
# original cell-line keys. These are not CYP3A5 star alleles/germline PGx.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop("Usage: export_nci60_cellminer_mutations.R <source_response.csv> <out_mutation.csv>")
if (!requireNamespace("rcellminerData", quietly=TRUE) ||
    !requireNamespace("Biobase", quietly=TRUE))
  stop("ABSTAIN: original Bioconductor mutation source missing")
source_version <- as.character(utils::packageVersion("rcellminerData"))
if (source_version != "2.34.0") stop(paste("ABSTAIN: unpinned CellMiner mutation package version",source_version))
utils::data("molData", package="rcellminerData", envir=environment())
if (!exists("molData", inherits=FALSE)) stop("ABSTAIN: source object not loaded")
# rcellminerData stores ExpressionSets in the MolData eSetList S4 slot.
# Access the documented slot directly: [[ on this MolData class is unavailable
# in the pinned release, even after loading rcellminer.
esets <- methods::slot(molData, "eSetList")
if (!is.list(esets) || !("mut" %in% names(esets)))
  stop("ABSTAIN: original CellMiner mutation ExpressionSet absent")
mut <- Biobase::exprs(esets[["mut"]])
if (!is.matrix(mut) || ncol(mut) != 60L || nrow(mut) < 1000L ||
    is.null(rownames(mut)) || is.null(colnames(mut)) ||
    anyDuplicated(rownames(mut)) || anyDuplicated(colnames(mut)))
  stop("ABSTAIN: unexpected original CellMiner mutation source matrix shape")
# Source v2.34.0 can expose nonbinary positive gene-level mutation
# intensities despite the released "binary" prose. Independent upstream
# rcellminer::getBinaryMutationData() defines presence by >0 at its final
# step (source R/getBinaryMutationData.R, exact public implementation).
# Here we only represent a *nonzero mutation signal*, never allele dosage.
observed <- mut[!is.na(mut)]
if (!length(observed) || any(!is.finite(observed)) || any(observed<0))
  stop("ABSTAIN: original source mutation values nonfinite or negative")
input_was_strictly_binary <- all(observed %in% c(0,1))
if (!input_was_strictly_binary)
  cat(sprintf("SOURCE_MUTATION_ENCODING nonbinary_positive_values=%d original_max=%s; apply published >0 presence rule; NO_ALLELE_DOSAGE_INFERENCE\n",
              sum(!(observed %in% c(0,1))),format(max(observed),digits=7)))
mut[!is.na(mut)] <- as.numeric(mut[!is.na(mut)] > 0)
if (any(!is.na(mut) & !(mut %in% c(0,1))))
  stop("ABSTAIN: source mutation presence conversion invalid")
# Reuse *exact* source-concordant 47 cell-line keys from existing audit.
responses <- utils::read.csv(args[[1]], check.names=FALSE, stringsAsFactors=FALSE)
if (!identical(names(responses), c("subject_id","cancer_class","y")) ||
    nrow(responses) != 47L || anyDuplicated(responses$subject_id))
  stop("ABSTAIN: original 47 distinct lineage-qualified drug-response lines required")
canon <- function(x) gsub("[^A-Z0-9]","",toupper(as.character(x)))
# Official CellMiner labels are tissue-prefixed (e.g. BR:MCF7, LC:NCI-H23).
# Tissue is independently cross-checked against IntLIM/Guo cancer_class.
# Only the official LC:A549/ATCC -> A549 spelling variant is a declared
# source-specific alias (see NCI CellMiner RNA-seq official cell-line help).
raw_labels <- colnames(mut)
prefix <- sub(":.*$","",raw_labels)
known_prefix <- c(BR="BREAST",CNS="CNS",CO="COLON",LE="LEUKEMIA",
                  ME="MELANOMA",LC="NSCLC",OV="OVARIAN",PR="PROSTATE",RE="RENAL")
if (any(!grepl("^[A-Z]+:",raw_labels)) || any(!(prefix %in% names(known_prefix))))
  stop("ABSTAIN: unknown original CellMiner source tissue prefix")
cell_name <- sub("^[A-Z]+:","",raw_labels)
cell_name[cell_name=="A549/ATCC" & prefix=="LC"] <- "A549"
rk <- canon(responses$subject_id); mk <- canon(cell_name)
if (any(!nzchar(rk)) || anyDuplicated(rk) || any(!nzchar(mk)) || anyDuplicated(mk))
  stop("ABSTAIN: ambiguous original mutation-to-drug source keys")
idx <- match(rk,mk)
if (anyNA(idx) || anyDuplicated(idx))
  stop(sprintf("ABSTAIN: original CellMiner source-name join %d/47, no free aliases allowed",sum(!is.na(idx))))
if (any(unname(known_prefix[prefix[idx]]) != responses$cancer_class))
  stop("ABSTAIN: cross-source mutation/cancer-lineage disagreement")
# Binary gene-level mutation presence is already algorithmically qualified by
# CellMiner. Original missing assay rows remain NA, NEVER treated as wildtype.
prespecified <- c("TP53","TOP2A","TOP2B","ABCB1","ABCC1","ABCG2","CBR1","CBR3")
genes <- ifelse(startsWith(rownames(mut),"mut"),substring(rownames(mut),4L),rownames(mut))
if (anyDuplicated(genes)) stop("ABSTAIN: ambiguous duplicate CellMiner mutation gene symbols")
available <- prespecified %in% genes
if (!("TP53" %in% genes)) stop("ABSTAIN: original TP53 mutation assay unavailable")
# Absence of a gene in mutation source is not evidence of wild type;
# omit absent genes, report missing canonical gene labels as metadata.
present <- prespecified[available]
flags <- t(mut[match(present,genes),idx,drop=FALSE])
if (nrow(flags) != 47L || ncol(flags) != length(present))
  stop("ABSTAIN: source-transposed mutation subset drift")
colnames(flags) <- paste0("mut_",present)
# No missing values in selected flags; we prefer abstaining to imputation.
if (anyNA(flags)) stop("ABSTAIN: original selected mutation flags contain missing calls")
# Mutation burden is informative but NOT a clinical germline genotype;
# refuse if mutation platform has incomplete calls, to avoid NA as wildtype.
all_calls <- mut[,idx,drop=FALSE]
if (anyNA(all_calls)) stop("ABSTAIN: missing genome-wide mutation calls; cannot compute trustworthy burden")
burden <- colSums(all_calls == 1)
out <- data.frame(subject_id=responses$subject_id, flags, mutation_burden=burden,
                  check.names=FALSE)
utils::write.csv(out,args[[2]],row.names=FALSE,quote=TRUE)
cat(sprintf("CELL_MINER_MUTATIONS source_version=%s original_n60=%d matched_n47=%d measured_mutation_genes=%d prior_gene_flags_available=%d TP53_mutated=%d missed_prior_gene_symbols=%s; SOMATIC_NOT_GERMLINE_PGX\n",
            source_version,ncol(mut),length(idx),nrow(mut),length(present),
            sum(flags[,"mut_TP53"]==1),paste(setdiff(prespecified,present),collapse=";")))
