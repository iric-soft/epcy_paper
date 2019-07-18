
#To build data/other
Rscript --vanilla ./src/script/other/biotype.r
Rscript --vanilla ./src/script/other/create_gene_length.r

# read count to tpm
Rscript --vanilla ./src/script/other/create_matrix_tpm.r

#To build data/TCGA_BRCA
Rscript --vanilla ./src/script/other/brca_matrix_readcounts.r
