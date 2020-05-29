

#epcy profile_rna -d ./data/design/leucegene3/30_t15_17/design.tsv \
#                 -m ./data/leucegene3/STAR_RSEM/readcounts.xls \
#                 --log --cpm \
#                 --ids ENSG00000255248.9 ENSG00000129682.16 ENSG00000259353.1 ENSG00000122824.11 ENSG00000229508.2 ENSG00000262831.1 ENSG00000183570.16 ENSG00000008853.16 ENSG00000113389.16 \
#                 -o  ./data/res/leucegene3/t15_17/cross

#epcy profile_rna -d ./data/design/leucegene3/30_t15_17/design.tsv \
#                -m ./data/leucegene3/STAR_RSEM/readcounts.xls \
#                --log --cpm \
#                --ids ENSG00000168004.9 ENSG00000162493.16 ENSG00000089820.15 \
#                -o  ./data/res/leucegene3/t15_17/top



#epcy profile_rna -d ./data/design/TCGA_LAML/all/design.tsv \
#                -m ./data/TCGA_LAML/htseq/readcounts.xls \
#                --log --cpm --subgroup sub_reduce --query t15_17 \
#                --ids ENSG00000130707.16 ENSG00000172543.6 ENSG00000279536.1 ENSG00000106991.12 ENSG00000120093.10 ENSG00000245685.6 ENSG00000270182.1 \
#                -o  ./data/res/TCGA_LAML/t15_17/




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

python3 -m pyres diff_pred \
  --outdir ../../data/res/ \


###########################################################
# comparaison of performance on random design
designs_random="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"
designs="30_inv16 30_t15_17"
method="deseq2 edger limma epcy"

#python3 -m pyres eval_random \
#  -p ../../data/design/leucegene3/ -r ../../data/design/leucegene3_random/ \
#  --biotype protein_coding --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
#  --methods ${method} -q "STAR_RSEM" \
#  --design ${designs} \
#  --design_random ${designs_random} \
#  --outdir ../../data/res/leucegene3 \

p_ss="0 0.01 0.02 0.04 0.06 0.08 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9"
rep="1 2 3 4 5 6 7 8 9 10"
designs="30_t15_17"
method="deseq2 edger voom epcy"
ids="ENSG00000117266.15 ENSG00000230749.5 ENSG00000168004.9 ENSG00000227268.4 ENSG00000008394.13 ENSG00000085514.16"
python3 -m pyres eval_ss \
  -p ../../data/design/leucegene3_ss/ \
  --methods ${method} -q "STAR_RSEM" \
  --design ${designs} \
  --outdir ../../data/res/leucegene3 \
  --p_ss ${p_ss} \
  --reps ${rep} \
  --ids ${ids} \

designs_random="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"
designs="1079_1 1046_2 216_11"
method="trend epcy"

#python3 -m pyres eval_random \
#  -p ../../data/design/10X/ -r ../../data/design/10X_random/ \
#  --methods ${method} -q "cellranger" \
#  --design ${designs} \
#  --design_random ${designs_random} \
#  --outdir ../../data/res/10X \
