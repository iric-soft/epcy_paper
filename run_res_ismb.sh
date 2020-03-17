
designs="28_inv16_vs_28 28_inv16 33_MLL 132_FLT3-ITD 139_NPM1_mut"

cd src/pyres

###########################################################
# comparaison of predicted performance
python3 -m pyres eval_cv -p ../../data/design/leucegene// \
  --cpm -m ../../data/leucegene/STAR/readcounts.xls \
  --loo --use_LR \
  --biotype protein_coding --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
  --design ${designs} \
  --outdir ../../data/res/leucegene \
  --top_values 1 3 5 10 50 100 200 500

###########################################################
# top 50 heatmap
#python3 -m pyres clust_exp -p ../../data/design/leucegene \
#  -m ../../data/leucegene/STAR/tpm.xls \
#  --biotype protein_coding --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
#  --design ${designs} \
#  --outdir ../../data/res/leucegene \
#  --top_values 50

###########################################################
# top 3 density
#for design in $designs
#do
#  python3 -m pyres density -p ../../data/design/leucegene \
#    --cpm -m ../../data/leucegene/STAR/readcounts.xls \
#    --design ${design} \
#    --biotype protein_coding --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
#    --outdir ../../data/res/leucegene \
#    --top 3
#done

###########################################################
# comparaison of predicted performance on TCGA BRCA
#num_fold="200"
#designs="104_triple_neg"
#python3 -m pyres eval_cv -p ../../data/design/TCGA_BRCA \
#  -m ../../data/TCGA_BRCA/htseq/readcounts.xls \
#  -q htseq -f ${num_fold} --use_LR \
#  --biotype protein_coding --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
#  --design ${designs} \
#  --outdir ../../data/res/TCGA_BRCA/ \
#  --top_values 1 3 5 10 50 100 200 500 \
#  --cpm
