#!/bin/bash
#LDE.sh ${path_software} ${path_design} ${path_counts} ${path_output} ${num_proc}

path_software=$1
path_design=$2
path_counts=$3
path_output=$4
num_proc=$5

echo "Rscript --vanilla ${path_software} ${path_design} ${path_counts} ${path_output} ${num_proc}"
