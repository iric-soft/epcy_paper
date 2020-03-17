designs="132_FLT3-ITD 139_NPM1_mut 28_inv16_vs_28 28_inv16 33_MLL 62_CK 9_EVI1 62_Inter 13_Mono5 126_Normal 30_t15_17 18_t8_21 13_Tri8 5_NUP98NSD1"

cd src/pyres
###########################################################
###########################################################
#### *** Leucegne ******

# cytogenetics
# python3 -m pyres clust_exp_umap -p ../../data/design/leucegene \
#  -m ../../data/leucegene/STAR/tpm.xls \
#  --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
#  --design ${designs} \
#  --outdir ../../data/res/leucegene \
#  --top_values 3 5 10 50 100 200 500
#  #--biotype protein_coding

#designs="28_inv16_vs_28 28_inv16 33_MLL 62_CK 9_EVI1 62_Inter 13_Mono5 126_Normal 30_t15_17 18_t8_21 13_Tri8 5_NUP98NSD1"
#python3 -m pyres clust_all_umap -p ../../data/design/leucegene \
#  -m ../../data/leucegene/STAR/tpm.xls \
#  --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
#  --design ${designs} \
#  --outdir ../../data/res/leucegene \
#  --ext pdf \
#  --mcc 0 \
#  --lfc 0
#  #--biotype protein_coding


###########################################################
# umap cluster
designs_cluster="cluster0 cluster1 cluster2 cluster3 cluster4 cluster5 cluster6 cluster7 cluster8 cluster9 cluster10 cluster11 cluster12"

#python3 -m pyres clust_all_umap -p ../../data/design/leucegene \
#  -m ../../data/leucegene/STAR/tpm.xls \
#  --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
#  --design ${designs_cluster} \
#  --outdir ../../data/res/leucegene_cluster \
#  --ext pdf \
#  --mcc 0 \
#  --lfc 0
#  #--biotype protein_coding

#for design in $designs_cluster
#do
#  python3 -m pyres density -p ../../data/design/leucegene \
#    --cpm -m ../../data/leucegene/STAR/readcounts.xls \
#    --design ${design} \
#    --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
#    --outdir ../../data/res/leucegene_cluster \
#    --top 5 \
#    --mcc 0 \
#    --lfc 0
#    #--biotype protein_coding
#
#done

###########################################################
###########################################################
#### *** Leucegne train/test ******

#for design in $designs_cluster
#do
#  python3 -m pyres density -p ../../data/design/leucegene \
#    --cpm -m ../../data/leucegene/STAR/readcounts.xls \
#    --design ${design} \
#    --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
#    --outdir ../../data/res/leucegene_cluster \
#    --top 5 \
#    --mcc 0 \
#    --lfc 0
#    #--biotype protein_coding
#
#done

python3 -m pyres eval_tt -p ../../data/design/ \
  -m ../../data/leucegene_v2/STAR/tpm.xls \
  --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
  --outdir ../../data/res/eval_tt \
  --top_values 1 2 3 5 10 20 50 100 200 500 1000 2000 \
  --pvalues 0.05 \
  --mcc 0 \
  --lfc 0
  # --biotype protein_coding

###########################################################
###########################################################
#### *** Leucegne v2 ******

#designs="31_inv16 48_MLL 92_CK 13_EVI1 87_Inter 23_Mono5 275_NK 30_t15_17 18_t8_21 18_Tri8 6_NUP98NSD1"
#python3 -m pyres clust_all_umap -p ../../data/design/leucegene_v2 \
#  -m ../../data/leucegene_v2/STAR/tpm.xls \
#  --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
#  --design ${designs} \
#  --outdir ../../data/res/leucegene_v2 \
#  --ext pdf \
#  --mcc 0 \
#  --lfc 0
#  #--biotype protein_coding

#python3 -m pyres clust_exp_umap -p ../../data/design/leucegene_v2 \
#  -m ../../data/leucegene_v2/STAR/tpm.xls \
#  --biotype protein_coding --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
#  --design ${designs} \
#  --outdir ../../data/res/leucegene_v2 \
#  --top_values 3 5 10 50 100 200 500

#for design in $designs
#do
#  python3 -m pyres density -p ../../data/design/leucegene_v2 \
#    --cpm -m ../../data/leucegene/STAR/readcounts.xls \
#    --design ${design} \
#    --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
#    --outdir ../../data/res/leucegene_v2 \
#    --top 5 \
#    --mcc 0 \
#    --lfc 0
#    #--biotype "protein_coding antisense" \
#done


###########################################################
###########################################################
#### *** Leucegne jf ******

designs="31_inv16 53_MLL 90_CK 13_EVI1 80_Inter 22_Mono5 267_NK 30_t15_17 18_t8_21 18_Tri8 19_NUP98"
#python3 -m pyres clust_all_umap -p ../../data/design/leucegene_jf \
#  -m ../../data/leucegene_v2/STAR/tpm.xls \
#  --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
#  --design ${designs} \
#  --outdir ../../data/res/leucegene_jf \
#  --ext pdf \
#  --mcc 0 \
#  --lfc 0
#  #--biotype protein_coding



###########################################################
# umap cluster
designs_cluster="cluster0 cluster1 cluster2 cluster3 cluster4 cluster5 cluster6 cluster7 cluster8 cluster9 cluster10 cluster11 cluster12 cluster13 cluster14"

#python3 -m pyres clust_all_umap -p ../../data/design/leucegene_jf \
#  -m ../../data/leucegene_v2/STAR/tpm.xls \
#  --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
#  --design ${designs_cluster} \
#  --outdir ../../data/res/leucegene_jf_cluster \
#  --ext pdf \
#  --mcc 0 \
#  --lfc 0
#  #--biotype protein_coding

###########################################################
# comparaison of predicted performance on TCGA BRCA
num_fold="200"
designs="104_triple_neg"
