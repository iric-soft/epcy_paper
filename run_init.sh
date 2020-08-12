
## To build data/other
#Rscript --vanilla ./src/script/other/biotype.r
#Rscript --vanilla ./src/script/other/create_gene_length.r

## Read count to tpm
#Rscript --vanilla ./src/script/other/create_matrix_tpm.r

## To build data/TCGA_BRCA
#Rscript --vanilla ./src/script/other/brca_matrix_readcounts.r

## To build data/TCGA_LAML
#Rscript --vanilla ./src/script/other/tcga_aml_matrix_readcounts.r

## To build data/TARGET_AML
#Rscript --vanilla ./src/script/other/target_aml_matrix_readcounts.r

## Create STAG2 dataset
#python ./src/script/other/STAG2.py

## To create all design for all cohorts
#Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/TARGET_AML/ subgroup
#Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/10X/ subgroup
#Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene3/ subgroup
#Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/STAG2/ subgroup
#Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/STAG2_ko/ subgroup
#Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/STAG2_ko_all/ subgroup
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/STAG2_granulo/ subgroup
#Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/10X_FACS/ subgroup

## Create leucegene_random
#Rscript --vanilla ./src/script/other/gen_leucegene_random.r
#Rscript --vanilla ./src/script/other/gen_10X_random.R
#Rscript --vanilla ./src/script/other/gen_STAG2_random.r

## Create subsampling
Rscript --vanilla ./src/script/other/gen_rep_leucegene_subsampling.R
#Rscript --vanilla ./src/script/other/gen_leucegene_subsampling.R
#Rscript --vanilla ./src/script/other/gen_10X_subsampling.R
#Rscript --vanilla ./src/script/other/gen_STAG2_subsampling.R


##

## MNNIST fashion
#git clone git@github.com:zalandoresearch/fashion-mnist.git
#python ./src/script/other/mnist_fashion.py
#Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/mnist_fashion_train/ subgroup
##Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/mnist_fashion_test/ subgroup
