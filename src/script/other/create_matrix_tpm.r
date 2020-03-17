
#LIB
library(data.table)
library(here)

#VARIABLE
script_dir = here()

count2tpm <- function(all_count, genes_length) {
  genes_selected = genes_length[which(genes_length$ID %in% all_count$ID),]
  all_count_common = all_count[which(all_count$ID %in% genes_length$ID),]

  genes_selected = genes_selected[match(all_count_common$ID, genes_selected$ID),]

  count_only = all_count_common[,-which(colnames(all_count_common) == "ID")]
  count_norm = count_only / genes_selected$length

  tpm = apply(count_norm, 2,
    function(x) {
      total_mass = sum(x)
      tpm_value = (x / total_mass) * 1e6
    }
  )

  all_tpm = cbind(data.frame(ID=all_count_common$ID, stringsAsFactors=F), tpm)

  return(all_tpm)
}

#MAIN
dir_data_other = file.path(script_dir, "data", "other")
file_genes_length = file.path(dir_data_other, "GRCh38_84_genes_length.tsv")
genes_length = fread(file_genes_length, header=TRUE, stringsAsFactors=FALSE,
                     sep="\t", quote = "")
setDF(genes_length)

dir_data_leucegene= file.path(script_dir, "data", "leucegene", "STAR")
file_count = file.path(dir_data_leucegene, "readcounts.xls")
all_count = fread(file_count, header=TRUE, stringsAsFactors=FALSE,
                     sep="\t", quote = "")
setDF(all_count)

all_tpm = count2tpm(all_count, genes_length)

file_out = file.path(dir_data_leucegene, "tpm.xls")
write.table(all_tpm, file_out, quote=FALSE, row.names=FALSE, sep="\t")



dir_data_tcga_laml= file.path(script_dir, "data", "TCGA_LAML", "htseq")
file_count = file.path(dir_data_tcga_laml, "readcounts.xls")
all_count = fread(file_count, header=TRUE, stringsAsFactors=FALSE,
                     sep="\t", quote = "")
setDF(all_count)

all_tpm = count2tpm(all_count, genes_length)

file_out = file.path(dir_data_tcga_laml, "tpm.xls")
write.table(all_tpm, file_out, quote=FALSE, row.names=FALSE, sep="\t")



dir_data_tcga_laml= file.path(script_dir, "data", "TARGET_AML", "htseq")
file_count = file.path(dir_data_tcga_laml, "readcounts.xls")
all_count = fread(file_count, header=TRUE, stringsAsFactors=FALSE,
                     sep="\t", quote = "")
setDF(all_count)

all_tpm = count2tpm(all_count, genes_length)

file_out = file.path(dir_data_tcga_laml, "tpm.xls")
write.table(all_tpm, file_out, quote=FALSE, row.names=FALSE, sep="\t")




dir_data_tcga_brca= file.path(script_dir, "data", "TCGA_BRCA", "htseq")
file_count = file.path(dir_data_tcga_brca, "readcounts.xls")
all_count = fread(file_count, header=TRUE, stringsAsFactors=FALSE,
                     sep="\t", quote = "")
setDF(all_count)

all_tpm = count2tpm(all_count, genes_length)

file_out = file.path(dir_data_tcga_brca, "tpm.xls")
write.table(all_tpm, file_out, quote=FALSE, row.names=FALSE, sep="\t")
