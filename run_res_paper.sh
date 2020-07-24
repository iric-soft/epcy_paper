

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


ids="ENSG00000117266.15 ENSG00000162493.16 ENSG00000230749.5 ENSG00000168004.9 ENSG00000227268.4 ENSG00000008394.13 ENSG00000085514.16"
#epcy profile_rna \
#  -d ./data/design/leucegene3/30_t15_17/design.tsv \
#  -m ./data/leucegene3/STAR_RSEM/readcounts.xls \
#  --log --cpm \
#  -o ./data/res/leucegene3/eval_ss/30_t15_17/ \
#  --ids ${ids}

#ids="ENSG00000169429 ENSG00000008394 ENSG00000163682 ENSG00000115828 ENSG00000087086"
ids="ENSG00000086730 ENSG00000090382 ENSG00000204287 ENSG00000103313 ENSG00000233927 ENSG00000163694 ENSG00000161642"
#epcy profile_rna \
#  -d ./data/design/10X/1079_1/design.tsv \
#  -m ./data/10X/cellranger/readcounts.xls \
#  --log --cpm --strip \
#  -o ./data/res/10x/1079_1/ \
#  --ids ${ids}

ids="HLA-DR ENSG00000197956 ENSG00000245532 ENSG00000136826 ENSG00000165092 CD28 ENSG00000143546 ENSG00000204472"
#epcy profile_rna \
#  -d ./data/design/10X/1079_1/design.tsv \
#  -m ./data/10X/cellranger/readcounts.xls \
#  --log --cpm --strip \
#  -o ./data/res/10x/1079_1/profile_v2/ \
#  --ids ${ids}

ids="ENSG00000143546 ENSG00000163220 CD14 ENSG00000162444 ENSG00000160255 CD27 ENSG00000196924 ENSG00000078596 ENSG00000152518 CD4"
#epcy profile_rna \
#  -d ./data/design/10X/1079_1/design.tsv \
#  -m ./data/10X/cellranger/readcounts.xls \
#  --log --cpm --strip \
#  -o ./data/res/10X/eval_ss/1079_1/ \
#  --ids ${ids}

#epcy profile_rna \
#   -d ./data/design/10X/1046_2/design.tsv \
#   -m ./data/10X/cellranger/readcounts.xls \
#   --log --cpm --strip \
#   -o ./data/res/10X/eval_ss/565_3/ \
#   --ids CD4

#epcy profile_rna \
#  -d ./data/design/10X/1046_2/design.tsv \
#  -m ./data/10X/cellranger/readcounts.xls \
#  --log --cpm --strip \
#  -o ./data/res/10X/eval_ss/565_3/ \
#  --ids CD4

#epcy profile_rna \
#  -d ./data/design/10X/413_4/design.tsv \
#  -m ./data/10X/cellranger/readcounts.xls \
#  --log --cpm --strip \
#  -o ./data/res/10X/eval_ss/413_4/ \
#  --ids CD4

#epcy profile_rna \
#  -d ./data/design/10X/387_5/design.tsv \
#  -m ./data/10X/cellranger/readcounts.xls \
#  --log --cpm --strip \
#  -o ./data/res/10X/eval_ss/387_5/ \
#  --ids CD4

#epcy profile_rna \
#  -d ./data/design/10X/361_6/design.tsv \
#  -m ./data/10X/cellranger/readcounts.xls \
#  --log --cpm --strip \
#  -o ./data/res/10X/eval_ss/361_6/ \
#  --ids CD4

#epcy profile_rna \
#  -d ./data/design/10X/316_7/design.tsv \
#  -m ./data/10X/cellranger/readcounts.xls \
#  --log --cpm --strip \
#  -o ./data/res/10X/eval_ss/316_7/ \
#  --ids CD4

#epcy profile_rna \
#  -d ./data/design/10X/311_8/design.tsv \
#  -m ./data/10X/cellranger/readcounts.xls \
#  --log --cpm --strip \
#  -o ./data/res/10X/eval_ss/311_8/ \
#  --ids CD4

#epcy profile_rna \
#  -d ./data/design/10X/307_9/design.tsv \
#  -m ./data/10X/cellranger/readcounts.xls \
#  --log --cpm --strip \
#  -o ./data/res/10X/eval_ss/307_9/ \
#  --ids CD4

#epcy profile_rna \
#  -d ./data/design/10X/246_10/design.tsv \
#  -m ./data/10X/cellranger/readcounts.xls \
#  --log --cpm --strip \
#  -o ./data/res/10X/eval_ss/246_10/ \
#  --ids CD4

#epcy profile_rna \
#  -d ./data/design/10X/216_11/design.tsv \
#  -m ./data/10X/cellranger/readcounts.xls \
#  --log --cpm --strip \
#  -o ./data/res/10X/eval_ss/216_11/ \
#  --ids CD4


ids="ENSG00000204472 ENSG00000204482 ENSG00000135047 ENSG00000211751 ENSG00000187514 ENSG00000109971"
#epcy profile_rna \
#  -d ./data/design/10X/216_11/design.tsv \
#  -m ./data/10X/cellranger/readcounts.xls \
#  --log --cpm --strip \
#  -o ./data/res/10X/eval_ss/216_11/ \
#  --ids ${ids}

ids="ENSG00000149294 ENSG00000115523 ENSG00000105374 ENSG00000100453 ENSG00000011600"
epcy profile_rna \
  -d ./data/design/10X_FACS/8385_cd56_nk/design.tsv \
  -m ./data/10X_FACS/cellranger/readcounts.xls \
  --log --cpm --strip \
  -o ./data/res/10X_FACS/profile/8385_cd56_nk/ \
  --ids ${ids}

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
p_fdr="100 200 400 1000"

#python3 -m pyres eval_random \
#  -p ../../data/design/leucegene3/ -r ../../data/design/leucegene3_random/ \
#  --biotype protein_coding --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
#  --methods ${method} -q "STAR_RSEM" \
#  --pfdr ${p_fdr} \
#  --design ${designs} \
#  --design_random ${designs_random} \
#  --outdir ../../data/res/leucegene3 \

p_ss="0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9"
#p_ss="0 0.2 0.9"
rep="1 2 3 4 5 6 7 8 9 10"
designs="30_t15_17"
method="deseq2 edger voom epcy" # epcy_bagging"
ids="ENSG00000117266.15 ENSG00000162493.16 ENSG00000230749.5 ENSG00000168004.9 ENSG00000227268.4 ENSG00000008394.13 ENSG00000085514.16"
#python3 -m pyres eval_ss \
#  -p ../../data/design/leucegene3_ss/ \
#  --methods ${method} -q "STAR_RSEM" \
#  --design ${designs} \
#  --outdir ../../data/res/leucegene3 \
#  --p_ss ${p_ss} \
#  --reps ${rep} \
#  --ids ${ids} \


designs_random="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"
designs="1079_1 216_11"
method="trend epcy mast"
p_fdr="0.001 0.01 0.05 0.1"

#python3 -m pyres eval_random \
#  -p ../../data/design/10X/ -r ../../data/design/10X_random/ \
#  --methods ${method} -q "cellranger" \
#  --pfdr ${p_fdr} \
#  --design ${designs} \
#  --design_random ${designs_random} \
#  --outdir ../../data/res/10X \


p_ss="0 0.1 0.2 0.3 0.4 0.5 0.6 0.7" # 0.8 0.9"
rep="1 2 3 4 5 6 7 8 9 10"
method="trend epcy mast"
design="1079_1"
ids="ENSG00000143546 ENSG00000163220 CD14 ENSG00000162444 ENSG00000160255 CD27 ENSG00000196924 ENSG00000078596 ENSG00000152518 CD4"
#python3 -m pyres eval_ss \
#  -p ../../data/design/10X_ss/ \
#  --methods ${method} -q "cellranger" \
#  --design ${design} \
#  --outdir ../../data/res/10X \
#  --p_ss ${p_ss} \
#  --reps ${rep} \
#  --ids ${ids} \


design="216_11"
ids="ENSG00000204472 ENSG00000204482 ENSG00000135047 ENSG00000211751 ENSG00000187514 ENSG00000109971"
#python3 -m pyres eval_ss \
#  -p ../../data/design/10X_ss/ \
#  --methods ${method} -q "cellranger" \
#  --design ${design} \
#  --outdir ../../data/res/10X \
#  --p_ss ${p_ss} \
#  --reps ${rep} \
#  --ids ${ids} \



designs_random="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"
designs="4458_ko"
method="trend epcy"

#python3 -m pyres eval_random \
#  -p ../../data/design/STAG2/ -r ../../data/design/STAG2_random/ \
#  --methods ${method} -q "cellranger" \
#  --design ${designs} \
#  --design_random ${designs_random} \
#  --outdir ../../data/res/STAG2 \


p_ss="0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9"
rep="1 2 3 4 5 6 7 8 9 10"
method="trend epcy"
design="4458_ko"
#####ids=""
#python3 -m pyres eval_ss \
#  -p ../../data/design/STAG2_ss/ \
#  --methods ${method} -q "cellranger" \
#  --design ${design} \
#  --outdir ../../data/res/STAG2 \
#  --p_ss ${p_ss} \
#  --reps ${rep} \
#  --ids ${ids} \
