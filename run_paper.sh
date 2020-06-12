

type_exec="torque" # "slurm" "bash" "torque"
num_proc="4"
working_dir="/u/eaudemard/project/epcy_paper/"
#working_dir="./"

######################################
# RUN DEG and PEG analysis
######################################

#data_type="bulk"
#designs_leucegene="30_inv16 30_t15_17"
#path_design_leucegene="leucegene3"
#data_project="leucegene3"
#for src_data in STAR_RSEM
#do
#  for design in ${designs_leucegene}
#  do
#    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design_leucegene} ${data_project} ${src_data} ${type_exec} ${num_proc} ${data_type}
#  done
#done

#designs_leucegene="30_t15_17 30_inv16"
#p_subs="0 0.01 0.02 0.04 0.06 0.08 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9"
#nums="1 2 3 4 5 6 7 8 9 10"
#path_design_leucegene="leucegene3_ss"
#data_project="leucegene3"
#for src_data in STAR_RSEM #kallisto
#do
#  for design in ${designs_leucegene}
#  do
#    for p in ${p_subs}
#    do
#      for num in ${nums}
#      do
#        design_ss="${design}/${p}/${num}"
#        bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design_ss} ${path_design_leucegene} ${data_project} ${src_data} ${type_exec} ${num_proc} ${data_type}
#      done
#    done
#  done
#done

#designs_leucegene="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"
#path_design_leucegene="leucegene3_random"
#data_project="leucegene3"
#for src_data in STAR_RSEM #kallisto
#do
#  for design in ${designs_leucegene}
#  do
#    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design_leucegene} ${data_project} ${src_data} ${type_exec} ${num_proc} ${data_type}
#  done
#done


data_type="sc"
path_design="10X"
data_project="10X"
designs_10X="1046_2 1079_1 216_11 246_10 307_9 311_8 316_7 361_6 387_5 413_4 565_3"
#for src_data in cellranger
#do
#  for design in ${designs_10X}
#  do
#    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design} ${data_project} ${src_data} ${type_exec} ${num_proc} ${data_type}
#  done
#done



designs_leucegene="1079_1 216_11"
p_subs="0 0.01 0.02 0.04 0.06 0.08 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9"
nums="1 2 3 4 5 6 7 8 9 10"
path_design_leucegene="10X_ss"
data_project="10X"
#for src_data in cellranger
#do
#  for design in ${designs_leucegene}
#  do
#    for p in ${p_subs}
#    do
#      for num in ${nums}
#      do
#        design_ss="${design}/${p}/${num}"
#        bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design_ss} ${path_design_leucegene} ${data_project} ${src_data} ${type_exec} ${num_proc} ${data_type}
#      done
#    done
#  done
#done


designs_10X="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"
path_design="10X_random"
data_project="10X"
#for src_data in cellranger
#do
#  for design in ${designs_10X}
#  do
#    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design} ${data_project} ${src_data} ${type_exec} ${num_proc} ${data_type}
#  done
#done



data_type="sc"
path_design="STAG2"
data_project="STAG2"
designs_STAG2="ok"
for src_data in SEQC
do
  for design in ${designs_STAG2}
  do
    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design} ${data_project} ${src_data} ${type_exec} ${num_proc} ${data_type}
  done
done


p_subs="0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9"
nums="1 2 3 4 5 6 7 8 9 10"
designs_STAG2="ok"
path_design_leucegene="STAG2_ss"
data_project="STAG2"
for src_data in SEQC
do
  for design in ${designs_STAG2}
  do
    for p in ${p_subs}
    do
      for num in ${nums}
      do
        design_ss="${design}/${p}/${num}"
        bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design_ss} ${path_design_leucegene} ${data_project} ${src_data} ${type_exec} ${num_proc} ${data_type}
      done
    done
  done
done


designs_STAG2="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"
path_design="STAG2_random"
data_project="STAG2"
for src_data in SEQC
do
  for design in ${designs_STAG2}
  do
    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design} ${data_project} ${src_data} ${type_exec} ${num_proc} ${data_type}
  done
done
