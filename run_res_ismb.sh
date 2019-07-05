num_fold="loo" # "loo" "10" "100" ...
#designs="132_FLT3-ITD"
designs="28_inv16_vs_28 28_inv16 33_MLL 132_FLT3-ITD 139_NPM1_mut"

cd src/pyres

###########################################################
# Using logistic regression
python3 -m pyres eval_cv -p ../../data/design/leucegene \
  -m ../../data/leucegene/STAR/readcounts.xls \
  --loo --use_LR \
  --biotype protein_coding --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
  --design ${designs} \
  --outdir ../../data/res/leucegene/ \
  --top_values 1 3 5 10 50 100 200 \
  --cpm
#  #--n_datasets 3
#  #--shuffle_seeds 1234

###########################################################
# Using random forest
#python3 -m pyres eval_cv -p ../../data/design/leucegene \
#  -m ../../data/leucegene/STAR/readcounts.xls \
#  --loo --use_randF \
#  --biotype protein_coding --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
#  --design ${designs} \
#  --outdir ../../data/res/leucegene/ \
#  --top_values 1 3 5 10 50 100 200 \
#  --cpm
#  #--n_datasets 3
#  #--shuffle_seeds 1234

#for method in epcy deseq2 edger limma
#do
#  python3 -m pyres heatmap_cv -p ../../data/design/leucegene \
#    --loo --biotype protein_coding \
#    --design ${designs} \
#    --outdir ../../data/res/leucegene/heatmap_cv \
#    --top 10 --method ${method} \
#    --bf ../../data/other/GRCh38_84_genes_biotype.tsv
#done

for gene in ENSG00000272767.1 ENSG00000232431.3
do
  for tool in density log_reg
  do
    python3 -m pyres ${tool} -p ../../data/design/leucegene \
      -m ../../data/leucegene/STAR/readcounts.xls \
      --cpm --design 33_MLL \
      --outdir ../../data/res/leucegene/ \
      --id ${gene}
  done
done

for gene in ENSG00000133392.16
do
  for tool in density log_reg
  do
    python3 -m pyres ${tool} -p ../../data/design/leucegene \
      -m ../../data/leucegene/STAR/readcounts.xls \
      --cpm --design 28_inv16 \
      --outdir ../../data/res/leucegene/ \
      --id ${gene}
  done
done
