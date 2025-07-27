# ------------------------------------------------------------------------------
# Script Name: gen_10X_reduce.r
#
# Description:
#   This script generates multiple random subsampled versions of 10X Genomics
#   single-cell RNA-seq design and readcount files. For each specified mean sample size,
#   it creates 20 random subsets, writes the reduced design.tsv and readcounts.xls
#   files to new directories, and prints progress messages for monitoring execution.
#
# Usage:
#   Rscript --vanilla gen_10X_reduce.r
#
# Input:
#   - Full design file:   data/design/10X_FACS/all/design.tsv
#   - Full readcount file: data/10X_FACS/cellranger/readcounts.xls
#
# Output:
#   - For each mean sample size and iteration, reduced design and readcount files are written to:
#     data/design/10X_FACS_reduce_[mean_sample]_[i]/all/design.tsv
#     data/10X_FACS_reduce_[mean_sample]_[i]/cellranger/readcounts.xls
#
# Author: eric.audemard@umontreal.ca
# Date: [Date]
# ------------------------------------------------------------------------------

library(data.table)
library(here)

##VARIABLE
script_dir = here()

dir_design = file.path(script_dir, "data", "design", "10X_FACS", "all")
file_design = file.path(dir_design, "design.tsv")
dir_readcount = file.path(script_dir, "data", "10X_FACS", "cellranger")
file_readcount = file.path(dir_readcount, "readcounts.xls")

cat("Reading design file:", file_design, "\n")
design = fread(file_design, header = TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")
cat("Reading readcount file:", file_readcount, "\n")
readcount = fread(file_readcount, header = TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")

for (mean_sample in c(3000, 5000, 8000, 10000)) {
  # Set seed for reproducibility with paper results
  set.seed(42)
  cat("Processing mean sample size:", mean_sample, "\n")
  for (i in c(1:20)) {
    cat("  Iteration:", i, "\n")
    random_num <- round(rnorm(1, mean=mean_sample, sd=1000))
    cat("    Number of samples selected:", random_num, "\n")
    selected_samples = sample(1:length(design$sample), random_num, replace=F)
    selected_samples = sort(selected_samples)
    design_reduce = design[selected_samples,]

    dir_reduce = file.path(script_dir, "data", "design", paste("10X_FACS_reduce", mean_sample, i, sep="_"), "all")
    file_reduce = file.path(dir_reduce, "design.tsv")
    cat("    Writing reduced design to:", file_reduce, "\n")
    dir.create(dir_reduce, recursive = TRUE, showWarnings = FALSE)
    write.table(design_reduce, file_reduce, quote=FALSE, row.names=FALSE, sep="\t")
    
    readcount_reduce <- readcount[, .SD, .SDcols = colnames(readcount) %in% c("ID", design_reduce$sample)]
    dir_reduce = file.path(script_dir, "data", paste("10X_FACS_reduce", mean_sample, i, sep="_"), "cellranger")
    file_reduce = file.path(dir_reduce, "readcounts.xls")
    cat("    Writing reduced readcount to:", file_reduce, "\n")
    dir.create(dir_reduce, recursive = TRUE, showWarnings = FALSE)
    write.table(readcount_reduce, file_reduce, quote=FALSE, row.names=FALSE, sep="\t")
    rm(design_reduce)
    rm(readcount_reduce)

    gc()
  }
  cat("Finished mean sample size:", mean_sample, "\n")
}
cat("All reductions completed.\n")