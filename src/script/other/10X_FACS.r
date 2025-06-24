# ------------------------------------------------------------------------------
# Script Name: 10X_FACS.r
#
# Description:
#   This script loads, concatenates, and processes 10X Genomics single-cell RNA-seq
#   data from multiple cell types. It merges count matrices for all specified cell types,
#   filters out genes with zero counts across all cells, and writes the combined matrix
#   to a file. It also generates a design file mapping each cell to its cell type.
#
# Usage:
#   Rscript --vanilla 10X_FACS.r
#
# Input:
#   - 10X Genomics output directories for each cell type, organized as:
#     data/10X_FACS/{cell_type}_filtered_matrices_mex/hg19/
#     Each directory must contain barcodes.tsv, genes.tsv, and matrix.mtx files.
#
# Output:
#   - Combined read count matrix: data/10X_FACS/cellranger/readcounts.xls
#   - Design file: data/design/10X_FACS/all/design.tsv
#
# Author: eric.audemard@umontreal.ca
# ------------------------------------------------------------------------------

library(Matrix)
library(here)
library(purrr)
library(dplyr)
library(data.table)

##VARIABLE
script_dir = here()
set.seed(42)

types = c(
  "cd14", "naive_cytotoxic", "cytotoxic_t",
  "cd56_nk", "memory_t", "naive_t",
  "regulatory_t", "cd4_t", "cd34", "b_cells"
)


## Load and concat 10X matrix of all cells types 

load_mat_10X <- function(matrix_dir, type) {
  barcode.path <- file.path(matrix_dir, "barcodes.tsv")
  features.path <- file.path(matrix_dir, "genes.tsv")
  matrix.path <- file.path(matrix_dir, "matrix.mtx")
  mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat) = paste(type, barcode.names$V1, sep="_")
  rownames(mat) = feature.names$V1
  #print("################")
  #print(length(colnames(mat)))
  #print(length(rownames(mat)))

  df = as.data.frame(as.matrix(mat))
  df = cbind(ID = feature.names$V1, df)

  return(df)
}

read_all <- function(type) {
  matrix_dir = file.path(
    script_dir, "data", "10X_FACS",
    paste0(type, "_filtered_matrices_mex"),
    "hg19"
  )
  print(matrix_dir)
  return(load_mat_10X(matrix_dir, type))
}

all_mat = lapply(types, function(x) read_all(x))
all_in = all_mat %>% reduce(full_join, by = "ID")
all_in_filtred = all_in[apply(all_in[,-1], 1, function(x) !all(x==0)),]

dir_out = file.path(
  script_dir, "data", "10X_FACS",
  "cellranger"
)
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
write.table(all_in_filtred, file.path(dir_out, "readcounts.xls"), quote=FALSE, row.names=FALSE, sep="\t")


## Load and concat 10X design matrix

load_design_10X <- function(matrix_dir, type) {
  barcode.path <- file.path(matrix_dir, "barcodes.tsv")
  barcode.names = read.delim(barcode.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)

  df = data.frame(sample=paste(type, barcode.names$V1, sep="_"), subgroup=type, stringsAsFactors=F )

  return(df)
}

read_all_design <- function(type) {
  matrix_dir = file.path(
    script_dir, "data", "10X_FACS",
    paste0(type, "_filtered_matrices_mex"),
    "hg19"
  )
  print(matrix_dir)
  return(load_design_10X(matrix_dir, type))
}

all_design = lapply(types, function(x) read_all_design(x))
all_design = rbindlist(all_design)
head(all_design)
dir_out = file.path(
  script_dir, "data", "design", "10X_FACS",
  "all"
)
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
write.table(all_design, file.path(dir_out, "design.tsv"), quote=FALSE, row.names=FALSE, sep="\t")





#ids_reduce = sample(1:length(all_design$sample),  10000)

#dir_out = file.path(
#  script_dir, "data", "design", "10X_FACS_reduce",
#  "all"
#)
#all_design_red = all_design[ids_reduce, ]
#dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
#write.table(all_design_red, file.path(dir_out, "design.tsv"), quote=FALSE, row.names=FALSE, sep="\t")


#dir_out = file.path(
#  script_dir, "data", "10X_FACS_reduce",
#  "cellranger"
#)
#all_in_filtred_red = all_in_filtred[, c(1,which(colnames(all_in_filtred) %in% all_design_red$sample))]
#dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
#write.table(all_in_filtred_red, file.path(dir_out, "readcounts.xls"), quote=FALSE, row.names=FALSE, sep="\t")
