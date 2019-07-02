##NEED
# module add R/3.5.1
##RUN
# Rscript --vanilla ./src/script/other/brca_matrix_readcounts.r

##LIB
library(GenomicDataCommons)
library(data.table)
library(here)

##VARIABLE
script_dir = here()

##FUNCTION
create_matrix <- function(tcga_se, dir_out, file_out) {
  matrix = assay(tcga_se)
  clinical = colData(tcga_se)
  matrix = matrix[,match(colnames(matrix), clinical$file_id)]
  colnames(matrix) = clinical$submitter_id.main
  matrix = cbind(data.frame(ID=rownames(matrix)), all_count)

  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
  write.table(all_count, file_out, quote=FALSE, row.names=FALSE, sep="\t")
}

##MAIN
dir_out = file.path(script_dir, "data", "TCGA_BRCA", "htseq")

file_out = file.path(dir_out, "readcounts.xls")
tcga_se = gdc_rnaseq('TCGA-BRCA', 'HTSeq - Counts')
create_matrix(tcga_se, dir_out, file_out)

#TODO gdc_rnaseq('TCGA-BRCA', 'HTSeq - FPKM') doesn't work
#file_out = file.path(dir_out, "fpkm.xls")
#tcga_se = gdc_rnaseq('TCGA-BRCA', 'HTSeq - FPKM')
#create_matrix(tcga_se, dir_out, file_out)

dir_other = file.path(script_dir, "data", "other")
file_clinical = file.path(dir_other, "nationwidechildrens.org_clinical_patient_brca.txt")
clinical = fread(file_clinical, header = TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")
setDF(clinical)

tested = c("Positive", "Negative")
her2 = which(clinical$her2_status_by_ihc %in% tested)
er = which(clinical$pr_status_by_ihc %in% tested)
pr = which(clinical$er_status_by_ihc %in% tested)

her2_er = her2[which(her2 %in% er)]
her2_er_pr =  her2_er[which(her2_er %in% pr)]

triple_neg = c("Negative")
her2_tn = which(clinical$her2_status_by_ihc %in% triple_neg)
er_tn = which(clinical$pr_status_by_ihc %in% triple_neg)
pr_tn = which(clinical$er_status_by_ihc %in% triple_neg)

her2_er_tn = her2_tn[which(her2_tn %in% er_tn)]
her2_er_pr_tn = her2_er_tn[which(her2_er_tn %in% pr_tn)]

not_tn = her2_er_pr[which(!her2_er_pr %in% her2_er_pr_tn)]

design = data.frame(sample=clinical$bcr_patient_barcode[her2_er_pr_tn], subgroup="Query")
design = rbind(design, data.frame(sample=clinical$bcr_patient_barcode[not_tn], subgroup="Ref"))

dir_out = file.path(script_dir, "data", "design", "116_BRCA_tn_vs_606")
file_out = file.path(dir_out, "design.tsv")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
write.table(design, file_out, quote=FALSE, row.names=FALSE, sep="\t")
