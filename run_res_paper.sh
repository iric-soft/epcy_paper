#
# NEED:
# - epcy installed
# - arround 1 hours to run all this script


cd src/pyres
########################################################################
# Figure 2: Trends of significance values as a fonction of # of samples

p_ss="0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9"
rep="1 2 3 4 5 6 7 8 9 10"
designs="30_t15_17"
method="deseq2 edger voom epcy"
ids="ENSG00000162493.16 ENSG00000227268.4"

python3 -m pyres eval_ss \
  -p ../../data/design/leucegene3_rep_ss/ \
  --methods ${method} -q "STAR_RSEM" \
  --design ${designs} \
  --outdir ../../data/res/leucegene3_rep \
  --p_ss ${p_ss} \
  --reps ${rep} \
  --ids ${ids} \


##############################################################################
# Figure 3: comparaison of DEG vs PG, using bulk RNA

designs="30_t15_17" #"30_t15_17 30_inv16"
method="deseq2 edger voom epcy"
quantiles="0.9999 0.9995 0.999 0.995 0.99"

python3 -m pyres eval_bulk \
  -p ../../data/design/leucegene3/ \
  -m ../../data/leucegene3/STAR_RSEM/readcounts.xls --cpm \
  --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
  --methods ${method} -q "STAR_RSEM" \
  --quantiles ${quantiles} \
  --design ${designs} \
  --outdir ../../data/res/leucegene3/ \
  --ngenes 60564

cd ../..

cd src/pyres
##############################################################################
# Figure 4 and 5: scatter plot which compare EPCY vs DEG on single cell

cell_numbers="3000 5000 8000 10000"
rep="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"
method="wilcox mast trend epcy"
celltypes="cd14 cd56_nk cd34 cytotoxic_t b_cells naive_t memory_t regulatory_t cd4 naive_cytotoxic"

python3 -m pyres eval_sc \
  -p ../../data/design \
  --methods ${method} -q "cellranger" \
  --celltypes ${celltypes} \
  --outdir ../../data/res/ \
  --cellNumber ${cell_numbers} \
  --reps ${rep} \

cd ../..

########################################################################
# Plotting gene expression used in severals figures and in supplemtary
ids="ENSG00000270947.1 ENSG00000141485.16 ENSG00000078399.18 ENSG00000286179.1 ENSG00000162493.16 ENSG00000256951.1 ENSG00000227268.4 ENSG00000089820.15 ENSG00000163701.19 ENSG00000214548.18 ENSG00000143995.20 ENSG00000232046.7 ENSG00000279536.1 ENSG00000076706.17 ENSG00000255248.9 ENSG00000102287.19 ENSG00000078399.18 ENSG00000143995.20 ENSG00000266217.2"
epcy profile_rna \
  -d ./data/design/leucegene3/30_t15_17/design.tsv --condition subgroup \
  -m ./data/leucegene3/STAR_RSEM/readcounts.xls \
  --log --cpm  --no_density \
  -o ./data/res/leucegene3/profile/30_t15_17/ \
  --ids ${ids}
