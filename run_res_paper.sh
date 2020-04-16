
cd src/pyres

###########################################################
# comparaison of predicted performance
#designs="28_inv16_vs_28 28_inv16 33_MLL 132_FLT3-ITD 139_NPM1_mut 5_NUP98NSD1 9_EVI1 13_Mono5 13_Tri8 18_t8_21 30_t15_17 62_CK 62_Inter"
#methods="epcy deseq2 deseq2_pvalue edger edger_pvalue limma limma_pvalue"
#methods="epcy deseq2_pvalue edger_pvalue limma_pvalue"
#methods="epcy deseq2 edger limma"
#python3 -m pyres eval_auc -p ../../data/design/leucegene/ \
#  --biotype protein_coding --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
#  --design ${designs} \
#  --method ${methods} \
#  --outdir ../../data/res/leucegene \
#  --top_values 1 3 5 10 50 100 200 500

###########################################################
# comparaison of performance on random design
designs="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"

python3 -m pyres eval_random -p ../../data/design/leucegene_random/ \
  --biotype protein_coding --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
  --design ${designs} \
  --outdir ../../data/res/leucegene \
