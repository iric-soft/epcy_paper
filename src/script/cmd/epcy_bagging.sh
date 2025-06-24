#!/bin/bash
#epcy.sh ${path_design} ${path_tpm} ${path_output} ${num_proc}

path_design=$1
path_tpm=$2
path_output=$3
num_proc=$4

echo "epcy pred_rna --cpm --log --npv --ppv --tpr --tnr --fnr --fpr --for --ts --acc --f1 --condition subgroup -b 20 -t ${num_proc} -d ${path_design}/design.tsv -m ${path_tpm}/readcounts.xls -o ${path_output}/ --randomseed 42 && epcy qc -p ${path_output}/predictive_capability.tsv -o ${path_output}/"
