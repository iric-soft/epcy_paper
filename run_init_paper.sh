
## To create all design for all cohorts
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene3/ subgroup
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/10X_FACS/ subgroup
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/10X_FACS_reduce/ subgroup

## Create leucegene_random
Rscript --vanilla ./src/script/other/gen_leucegene_random.r
Rscript --vanilla ./src/script/other/gen_10X_random.R

## Create subsampling
Rscript --vanilla ./src/script/other/gen_rep_leucegene_subsampling.R
