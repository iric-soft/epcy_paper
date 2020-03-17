##NEED
# module add R/3.5.1
##RUN
# Rscript --vanilla ./src/script/other/tcga_aml_matrix_readcounts.r

##LIB
library(GenomicDataCommons)
library(SummarizedExperiment)
library(data.table)
library(here)

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
file_clinical = file.path(dir_other, "nationwidechildrens.org_clinical_patient_laml.txt")
clinical = fread(file_clinical, header = TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")
setDF(clinical)

inv16 = which(clinical$cytogenetic_abnormality_type %in% "inv(16)")
nk = which(clinical$cytogenetic_abnormality_type %in% "Normal")
t8_21 = which(clinical$cytogenetic_abnormality_type %in% "t(8;21)")
t15_17 = which(clinical$cytogenetic_abnormality_type %in% "t(15;17) and variants")
#take care for Complex
ck = which(clinical$cytogenetic_abnormality_type %like% "Complex")

design = data.frame(sample=clinical$bcr_patient_barcode, subgroup=clinical$cytogenetic_abnormality_type)
design$sample = gsub("-", "_", design$sample)

dir_out = file.path(script_dir, "data", "TCGA_LAML", "htseq")

file_out = file.path(dir_out, "readcounts.xls")
tcga_se = gdc_rnaseq('TCGA-LAML', 'HTSeq - Counts')
design = create_matrix(design, tcga_se, dir_out, file_out)


dir_out = file.path(script_dir, "data", "design", "TCGA_LAML", "all")
file_out = file.path(dir_out, "design.tsv")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
write.table(design, file_out, quote=FALSE, row.names=FALSE, sep="\t")

#TODO gdc_rnaseq('TCGA-BRCA', 'HTSeq - FPKM') doesn't work
#file_out = file.path(dir_out, "fpkm.xls")
#tcga_se = gdc_rnaseq('TCGA-BRCA', 'HTSeq - FPKM')
#create_matrix(tcga_se, dir_out, file_out)
