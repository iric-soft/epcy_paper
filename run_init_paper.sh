
## To create all design for all cohorts
echo "Generate subgroup for leucegene 3" 
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene3/ subgroup
echo "Generate subgroup for 10X FACS"
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/10X_FACS/ subgroup
echo "Generate subgroup for 10X FACS reduce"
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/10X_FACS_reduce/ subgroup

## Create leucegene_random
echo "Generate random for leucegene 3"
Rscript --vanilla ./src/script/other/gen_leucegene_random.r
echo "Generate random for 10X FACS reduce"
Rscript --vanilla ./src/script/other/gen_10X_random.r

## Create subsampling
echo "Generate subsampling for leucegene3"
Rscript --vanilla ./src/script/other/gen_rep_leucegene_subsampling.r
