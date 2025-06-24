
## To create all design for all cohorts
echo "Generate subgroup for leucegene 3"
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/leucegene3/ subgroup
echo "Generate subgroup for 10X FACS"
Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/10X_FACS/ subgroup
echo "Generate subgroup for 10X FACS reduce"
Rscript --vanilla ./src/script/other/gen_10X_reduce.r
for num_sample in  3000 5000 8000 10000; do
  for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do
    Rscript --vanilla ./src/script/other/gen_design_from_all.r ./data/design/10X_FACS_reduce_${num_sample}_$i/ subgroup
  done
done

## Create subsampling
echo "Generate subsampling for leucegene3"
Rscript --vanilla ./src/script/other/gen_rep_leucegene_subsampling.r
