library(here)

#VARIABLE
script_dir = here()

# for GRCh38_84 host="http://mar2016.archive.ensembl.org"
# check here https://useast.ensembl.org/Help/ArchiveList to find host for an other version
mart = biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="http://mar2016.archive.ensembl.org")
genes = biomaRt::getBM(attributes = c("ensembl_gene_id", "version", "external_gene_name", "gene_biotype"), mart=mart)
genes$ensembl_gene_id = paste(genes$ensembl_gene_id, genes$version, sep=".")
genes = genes[,-2]

dir_out = file.path(script_dir, "data", "other")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
file_out = file.path(dir_out, "GRCh38_84_genes_biotype.tsv")
write.table(genes, file_out, quote=FALSE, row.names=FALSE, sep="\t")
