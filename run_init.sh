
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

## To create all design for all cohorts
#Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene/ subgroup
#Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_v2/ subgroup
#Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_jf/ JF
#Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/TARGET_AML/ subgroup

## Create leucegene_random
#Rscript --vanilla ./src/script/other/gen_leucegene_random.r

## Create leucegene_ss
Rscript --vanilla ./src/script/other/gen_leucegene_subsampling.r

## Build leucegene_test and leucegene_train
#Rscript --vanilla ./src/script/other/gen_leucegene_tt.r
#Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_train/ subgroup

## MNNIST fashion
#git clone git@github.com:zalandoresearch/fashion-mnist.git
#python ./src/script/other/mnist_fashion.py
#Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/mnist_fashion_train/ subgroup
##Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/mnist_fashion_test/ subgroup

## Build leucegene_test and leucegene_learning_curve
#Rscript --vanilla ./src/script/other/gen_leucegene_learning_curve.r
#Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_learning_curve/comp subgroup


## Build leucegene_test and leucegene_comp_tt
#Rscript --vanilla ./src/script/other/gen_leucegene_comp_tt.r
#Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt/comp subgroup FALSE
#Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt/test_train subgroup FALSE
