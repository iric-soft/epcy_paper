library(Matrix)
library(here)
library(purrr)
library(dplyr)
library(data.table)

##VARIABLE
script_dir = here()


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





types = c(
  "cd14", "naive_cytotoxic", "cytotoxic_t",
  "cd56_nk", "memory_t", "naive_t",
  "regulatory_t", "cd4_t", "cd34", "b_cells")

all_mat = lapply(types, function(x) read_all(x))
all_in = all_mat %>% reduce(full_join, by = "ID")
all_in_filtred = all_in[apply(all_in[,-1], 1, function(x) !all(x==0)),]

dir_out = file.path(
  script_dir, "data", "10X_FACS",
  "cellranger"
)
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
write.table(all_in_filtred, file.path(dir_out, "readcounts.xls"), quote=FALSE, row.names=FALSE, sep="\t")

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
