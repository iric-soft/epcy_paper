#
# NEED:
# - epcy installed
# - arround 1 hours to run all this script


cd src/pyres
##########################################################################
# Figure 1: Synthetic examples to illustrate the diff. between DEG and PG

#python3 -m pyres diff_pred --outdir ../../data/res/ \

########################################################################
# Figure 3: Trends of significance values as a fonction of # of samples

p_ss="0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9"
rep="1 2 3 4 5 6 7 8 9 10"
designs="30_t15_17"
method="deseq2 edger voom epcy"
ids="ENSG00000117266.15 ENSG00000162493.16 ENSG00000230749.5 ENSG00000168004.9 ENSG00000227268.4 ENSG00000008394.13 ENSG00000085514.16"

#python3 -m pyres eval_ss \
#  -p ../../data/design/leucegene3_rep_ss/ \
#  --methods ${method} -q "STAR_RSEM" \
#  --design ${designs} \
#  --outdir ../../data/res/leucegene3_rep \
#  --p_ss ${p_ss} \
#  --reps ${rep} \
#  --ids ${ids} \


##############################################################################
# Figure 4: comparaison of performance on random design and select thresholds
# Figure 5: scatter plot which compare EPCY vs DEG
designs_random="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"
designs="30_t15_17 30_inv16"
method="deseq2 edger voom epcy"
e_fpr="0.00001 0.0001 0.001 0.01"

#python3 -m pyres eval_random \
#  -p ../../data/design/leucegene3/ -r ../../data/design/leucegene3_random/ \
#  --biotype protein_coding --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
#  --methods ${method} -q "STAR_RSEM" \
#  --efpr ${e_fpr} --top 20 \
#  --design ${designs} \
#  --design_random ${designs_random} \
#  --outdir ../../data/res/leucegene3 \
#  --ngenes 60564

cd ../..
#Create Fig. 5B), and selecte some genes to highlight
python3 ./src/script/other/t15_17.py
python3 ./src/script/other/inv16.py

cd src/pyres
##############################################################################
# Figure 6: scatter plot which compare EPCY vs DEG on single cell
designs_random="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"
designs="287_cd14 858_cd56_nk 948_cd34 1016_cytotoxic_t 1058_b_cells 1077_naive_t 1093_memory_t 1108_regulatory_t 1217_cd4_t 1338_naive_cytotoxic"
method="mast trend epcy"
e_fpr="0.0001 0.001 0.01 0.05"

#python3 -m pyres eval_random \
#  -p ../../data/design/10X_FACS_reduce/ -r ../../data/design/10X_FACS_reduce_random/ \
#  --methods ${method} -q "cellranger" \
#  --efpr ${e_fpr} --top 50 \
#  --design ${designs} \
#  --design_random ${designs_random} \
#  --outdir ../../data/res/10X_FACS_reduce \
#  --ngenes 21952

cd ../..
#Create UMAP on some selected genes and Fig. 6B)
#python3 ./src/script/other/10X_FACS_red_res_umap.py


########################################################################
# Plotting gene expression used in severals figures and in supplemtary
ids="ENSG00000141485.16 ENSG00000078399.18 ENSG00000286179.1 ENSG00000162493.16 ENSG00000256951.1 ENSG00000227268.4 ENSG00000089820.15 ENSG00000163701.19 ENSG00000214548.18 ENSG00000143995.20 ENSG00000232046.7 ENSG00000279536.1 ENSG00000076706.17 ENSG00000255248.9 ENSG00000102287.19 ENSG00000078399.18 ENSG00000143995.20 ENSG00000266217.2"
epcy profile_rna \
  -d ./data/design/leucegene3/30_t15_17/design.tsv \
  -m ./data/leucegene3/STAR_RSEM/readcounts.xls \
  --log --cpm  --no_density \
  -o ./data/res/leucegene3/profile/30_t15_17/ \
  --ids ${ids}

ids="ENSG00000171388.12 ENSG00000147488.11 ENSG00000133392.18 ENSG00000227502.3 ENSG00000188153.13 ENSG00000228836.8"
#epcy profile_rna \
#  -d ./data/design/leucegene3/30_inv16/design.tsv \
#  -m ./data/leucegene3/STAR_RSEM/readcounts.xls \
#  --log --cpm  --no_density \
#  -o ./data/res/leucegene3/profile/30_inv16/ \
#  --ids ${ids}

ids="ENSG00000147403 ENSG00000174059 ENSG00000167526 ENSG00000115268 ENSG00000008988 ENSG00000196126 ENSG00000019582 ENSG00000231389 ENSG00000204287 ENSG00000113389 ENSG00000233968 ENSG00000237819 ENSG00000113389 ENSG00000163106 ENSG00000204287 ENSG00000124766 ENSG00000125691 ENSG00000171858 ENSG00000165092 ENSG00000174099 ENSG00000156508 ENSG00000204472 ENSG00000223609 ENSG00000163751"
epcy profile_rna \
  -d ./data/design/10X_FACS_reduce/948_cd34/design.tsv \
  -m ./data/10X_FACS_reduce/cellranger/readcounts.xls \
  --log --cpmed --no_density --violin \
  -o ./data/res/10X_FACS_reduce/948_cd34/violin \
  --ids ${ids}
