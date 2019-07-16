#designs="33_MLL"
designs="28_inv16_vs_28 28_inv16 33_MLL 132_FLT3-ITD 139_NPM1_mut"

cd src/pyres

###########################################################
# Using logistic regression
#python3 -m pyres eval_cv -p ../../data/design/leucegene/ \
#  --cpm -m ../../data/leucegene/STAR/readcounts.xls \
#  --loo --use_LR \
#  --biotype protein_coding --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
#  --design ${designs} \
#  --outdir ../../data/res/leucegene \
#  --top_values 1 3 5 10 50 100 200 500

#python3 -m pyres eval_cv -p ../../data/design/leucegene/ \
#  --cpm -m ../../data/leucegene/STAR/readcounts.xls \
#  --loo --use_LR\
#  --biotype protein_coding --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
#  --design ${designs} \
#  --outdir ../../data/res/leucegene/pvalue \
#  --top_values 1 3 5 10 50 100 200 500 \
#  --pvalues 0.05 0.01 0.00001 0.0000000001 0.000000000000001 0.00000000000000000001

#python3 -m pyres clust_exp -p ../../data/design/leucegene \
#  -m ../../data/leucegene/STAR/tpm.xls \
#  --biotype protein_coding --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
#  --design ${designs} \
#  --outdir ../../data/res/leucegene \
#  --top_values 50

#python3 -m pyres clust_exp -p ../../data/design/leucegene \
#  -m ../../data/leucegene/STAR/readcounts.xls \
#  --biotype protein_coding --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
#  --design ${designs} \
#  --outdir ../../data/res/leucegene/ \
#  --top_values 3 5 10 50 100 200 500 \
#  --cpm

for design in $designs
do
  python3 -m pyres density -p ../../data/design/leucegene \
    --cpm -m ../../data/leucegene/STAR/readcounts.xls \
    --design ${design} \
    --biotype protein_coding --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
    --outdir ../../data/res/leucegene \
    --top 3
done

#for top in 10 50 200
#do
#  for method in epcy deseq2 edger limma
#  do
#    python3 -m pyres heatmap_cv -p ../../data/design/leucegene \
#      --loo --biotype protein_coding \
#      --design ${designs} \
#      --outdir ../../data/res/leucegene/heatmap_cv \
#      --top ${top} --method ${method} \
#      --bf ../../data/other/GRCh38_84_genes_biotype.tsv
#  done
#done

#for gene in ENSG00000272767.1 ENSG00000232431.3 ENSG00000277435.1 ENSG00000188626.6 ENSG00000271616.1
#do
#  for tool in log_reg
#  do
#    python3 -m pyres ${tool} -p ../../data/design/leucegene \
#      -m ../../data/leucegene/STAR/readcounts.xls \
#      --cpm --design 33_MLL \
#      --outdir ../../data/res/leucegene/ \
#      --id ${gene}
#  done
#done

#for gene in ENSG00000133392.16
#do
#  for tool in density log_reg
#  do
#    python3 -m pyres ${tool} -p ../../data/design/leucegene \
#      -m ../../data/leucegene/STAR/readcounts.xls \
#      --cpm --design 28_inv16 \
#      --outdir ../../data/res/leucegene/ \
#      --id ${gene}
#  done
#done


#num_fold="200"
#designs="104_triple_neg"
#for pred in --use_LR #--use_randF
#do
#  python3 -m pyres eval_cv -p ../../data/design/TCGA_BRCA \
#    -m ../../data/TCGA_BRCA/htseq/readcounts.xls \
#    -q htseq -f ${num_fold} ${pred} \
#    --biotype protein_coding --bf ../../data/other/GRCh38_84_genes_biotype.tsv \
#    --design ${designs} \
#    --outdir ../../data/res/TCGA_BRCA/ \
#    --top_values 1 3 5 10 50 100 200 500 \
#    --pvalue 0.01 \
#    --cpm
#done
