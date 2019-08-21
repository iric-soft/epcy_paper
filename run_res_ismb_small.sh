designs="28_inv16_vs_28"

cd src/pyres

###########################################################
# comparaison of predicted performance
python3 -m pyres eval_cv -p ../../data/design/leucegene/ \
  --cpm -m ../../data/leucegene/STAR/readcounts.xls \
  --loo --use_LR \
  --biotype protein_coding --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
  --design ${designs} \
  --outdir ../../data/res/leucegene \
  --top_values 1 3 5 10 50 100 200 500

###########################################################
# top 50 heatmap
python3 -m pyres clust_exp -p ../../data/design/leucegene \
  -m ../../data/leucegene/STAR/tpm.xls \
  --biotype protein_coding --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
  --design ${designs} \
  --outdir ../../data/res/leucegene \
  --top_values 50

###########################################################
# top 3 density
for design in $designs
do
  python3 -m pyres density -p ../../data/design/leucegene \
    --cpm -m ../../data/leucegene/STAR/readcounts.xls \
    --design ${design} \
    --biotype protein_coding --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
    --outdir ../../data/res/leucegene \
    --top 3
done
