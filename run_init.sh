
#To build data/other
Rscript --vanilla ./src/script/other/biotype.r
Rscript --vanilla ./src/script/other/create_gene_length.r

#To build data/TCGA_BRCA
Rscript --vanilla ./src/script/other/brca_matrix_readcounts.r
