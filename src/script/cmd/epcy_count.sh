#!/bin/bash
#epcy.sh ${path_design} ${path_tpm} ${path_output} ${num_proc}

path_design=$1
path_tpm=$2
path_output=$3
num_proc=$4

echo "epcy pred --cpm -l 0 -e 0 -t ${num_proc} -d ${path_design}/design.tsv -m ${path_tpm}/readcounts.xls -o ${path_output}/"
