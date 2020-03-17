##NEED
# module add R/3.5.1
##RUN
# Rscript --vanilla ./src/script/other/target_aml_matrix_readcounts.r

##LIB
library(GenomicDataCommons)
library(SummarizedExperiment)
library(data.table)
library(here)
library(xlsx)

##VARIABLE
script_dir = here()

##FUNCTION
create_matrix <- function(design, tcga_se, dir_out, file_out) {
  matrix = assay(tcga_se)
  clinical = colData(tcga_se)
  matrix = matrix[,match(colnames(matrix), clinical$file_id)]
  matrix = matrix[which(rowSums(matrix) != 0),]

  #remove samples sequenced severals times
  ids_duplicate = clinical$submitter_id.main[duplicated(clinical$submitter_id.main)]
  ids2del = which(!clinical$submitter_id.main %in% ids_duplicate)
  clinical = clinical[ids2del,]
  matrix = matrix[,ids2del]
  ids_duplicate = gsub("-", "_", ids_duplicate)
  design = design[which(!design$sample %in% ids_duplicate), ]

  colnames(matrix) = clinical$submitter_id.main
  matrix = cbind(data.frame(ID=rownames(matrix)), matrix)
  colnames(matrix) = gsub("-", "_", colnames(matrix))

  #emove sample not in matrix (witout htseq count)
  design = design[which(design$sample %in% colnames(matrix)), ]

  col_selected = which(colnames(matrix) %in% design$sample)
  matrix = matrix[,c(1, col_selected)]

  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
  write.table(matrix, file_out, quote=FALSE, row.names=FALSE, sep="\t")

  return(design)
}

##MAIN

dir_other = file.path(script_dir, "data", "other")
file_clinical = file.path(dir_other, "TARGET_AML_ClinicalData_Discovery_20181213.xlsx")
clinical = read.xlsx(file_clinical, sheetIndex = 1,header = TRUE, stringsAsFactors=FALSE)

## Sample=TARGET USI (column 1), subgroup=Primary Cytogenetic Code (column 39)
design = data.frame(sample=clinical[,1], subgroup=clinical[,39], stringsAsFactors=FALSE)
design$sample = gsub("-", "_", design$sample)
design$subgroup[which(design$subgroup == "inv(16)")] = "inv16"
design$subgroup[which(design$subgroup == "t(8;21)")] = "t8_21"

dir_out = file.path(script_dir, "data", "TARGET_AML", "htseq")

file_out = file.path(dir_out, "readcounts.xls")
tcga_se = gdc_rnaseq('TARGET-AML', 'HTSeq - Counts')
design = create_matrix(design, tcga_se, dir_out, file_out)


dir_out = file.path(script_dir, "data", "design", "TARGET_AML", "all")
file_out = file.path(dir_out, "design.tsv")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
write.table(design, file_out, quote=FALSE, row.names=FALSE, sep="\t")
