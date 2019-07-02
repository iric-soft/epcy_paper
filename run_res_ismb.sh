num_fold="loo" # "loo" "10" "100" ...
#designs="132_FLT3-ITD"
designs="28_inv16_vs_28 28_inv16 33_MLL 132_FLT3-ITD 139_NPM1_mut"



cd src/pyres

python3 -m pyres eval_cv -p ../../data \
  -m ../../data/leucegene/STAR/readcounts.xls \
  --loo --use_LR \
  --biotype protein_coding --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
  --design ${designs} \
  --outdir ../../data/res/leucegene/ \
  --top_values 1 3 5 10 50 100 200 \
  --cpm
  #--n_datasets 3
  #--shuffle_seeds 1234

python3 -m pyres eval_cv -p ../../data \
  -m ../../data/leucegene/STAR/readcounts.xls \
  --loo --use_randF \
  --biotype protein_coding --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
  --design ${designs} \
  --outdir ../../data/res/leucegene/ \
  --top_values 1 3 5 10 50 100 200 \
  --cpm
  #--n_datasets 3
  #--shuffle_seeds 1234
