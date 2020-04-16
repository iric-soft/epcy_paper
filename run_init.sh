
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


## Build leucegene_test and leucegene_comp_tt
Rscript --vanilla ./src/script/other/gen_leucegene_comp_tt.r
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/1/comp subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/1/test_train subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/2/comp subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/2/test_train subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/3/comp subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/3/test_train subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/4/comp subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/4/test_train subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/5/comp subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/5/test_train subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/6/comp subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/6/test_train subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/7/comp subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/7/test_train subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/8/comp subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/8/test_train subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/9/comp subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/9/test_train subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/10/comp subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/10/test_train subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/11/comp subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/11/test_train subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/12/comp subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/12/test_train subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/13/comp subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/13/test_train subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/14/comp subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/14/test_train subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/15/comp subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/15/test_train subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/16/comp subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/16/test_train subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/17/comp subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/17/test_train subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/18/comp subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/18/test_train subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/19/comp subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/19/test_train subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/20/comp subgroup FALSE
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene_comp_tt_v2/20/test_train subgroup FALSE
