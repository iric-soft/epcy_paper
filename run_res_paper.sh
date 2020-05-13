
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


epcy profile_rna -d ./data/design/leucegene/30_t15_17/design.tsv \
                 -m ./data/leucegene/STAR/readcounts.xls \
                 --log --cpm \
                 --ids ENSG00000130707.17 ENSG00000172543.7 ENSG00000279536.1 ENSG00000106991.13 ENSG00000120093.11 ENSG00000245685.6 ENSG00000270182.1 \
                 -o  ./data/res/leucegene/t15_17/

epcy profile_rna -d ./data/design/TCGA_LAML/all/design.tsv \
                -m ./data/TCGA_LAML/htseq/readcounts.xls \
                --log --cpm --subgroup sub_reduce --query t15_17 \
                --ids ENSG00000130707.16 ENSG00000172543.6 ENSG00000279536.1 ENSG00000106991.12 ENSG00000120093.10 ENSG00000245685.6 ENSG00000270182.1 \
                -o  ./data/res/TCGA_LAML/t15_17/



###########################################################
# comparaison of performance on random design
designs_random="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"
designs="28_inv16 30_t15_17 62_Inter"

#python3 -m pyres eval_random \
#  -p ../../data/design/leucegene/ -r ../../data/design/leucegene_random/ \
#  --biotype protein_coding --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
#  --design ${designs} \
#  --design_random ${designs_random} \
#  --outdir ../../data/res/leucegene \


designs_random="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"
designs="1079_1 1046_2 216_11"
method="trend epcy"

  python3 -m pyres eval_random \
    -p ../../data/design/10X/ -r ../../data/design/10X_random/ \
    --methods ${method} -q "cellranger" \
    --design ${designs} \
    --design_random ${designs_random} \
    --outdir ../../data/res/10X \
