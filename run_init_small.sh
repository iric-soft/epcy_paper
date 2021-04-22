## Create leucegene_random
echo "Generate random for small_leucegene"
Rscript --vanilla ./src/script/other/gen_small_leucegene_random.r

## Create subsampling
echo "Generate subsampling for small_leucegene"
Rscript --vanilla ./src/script/other/gen_rep_small_leucegene_subsampling.r
