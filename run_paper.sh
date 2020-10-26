

type_exec="torque" # "slurm" "bash" "torque"
num_proc="4"
working_dir="/u/eaudemard/project/epcy_paper/"
#working_dir="./"

######################################
# RUN DEG and PEG analysis
######################################

data_type="bulk"
designs_leucegene="30_inv16 30_t15_17"
path_design_leucegene="leucegene3"
data_project="leucegene3"
walltime_epcy="24:00:00"
mem_epcy="50Gb"
for src_data in STAR_RSEM
do
  for design in ${designs_leucegene}
  do
    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design_leucegene} ${data_project} ${src_data} ${type_exec} ${num_proc} ${data_type} ${walltime_epcy} ${mem_epcy}
  done
done

designs_leucegene="30_t15_17 30_inv16"
p_subs="0 0.01 0.02 0.04 0.06 0.08 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9"
nums="1 2 3 4 5 6 7 8 9 10"
path_design_leucegene="leucegene3_ss"
data_project="leucegene3"
#for src_data in STAR_RSEM #kallisto
#do
#  for design in ${designs_leucegene}
#  do
#    for p in ${p_subs}
#    do
#      for num in ${nums}
#      do
#        design_ss="${design}/${p}/${num}"
#        bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design_ss} ${path_design_leucegene} ${data_project} ${src_data} ${type_exec} ${num_proc} ${data_type} ${walltime_epcy} ${mem_epcy}
#      done
#    done
#  done
#done

designs_leucegene="30_t15_17 30_inv16"
p_subs="0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9"
nums="1 2 3 4 5 6 7 8 9 10"
path_design_leucegene="leucegene3_rep_ss"
data_project="leucegene3_rep"
#for src_data in STAR_RSEM
#do
#  for design in ${designs_leucegene}
#  do
#    for p in ${p_subs}
#    do
#      for num in ${nums}
#      do
#        design_ss="${design}/${p}/${num}"
#        bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design_ss} ${path_design_leucegene} ${data_project} ${src_data} ${type_exec} ${num_proc} ${data_type} ${walltime_epcy} ${mem_epcy}
#      done
#    done
#  done
#done


designs_leucegene="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"
path_design_leucegene="leucegene3_random"
data_project="leucegene3"
#for src_data in STAR_RSEM #kallisto
#do
#  for design in ${designs_leucegene}
#  do
#    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design_leucegene} ${data_project} ${src_data} ${type_exec} ${num_proc} ${data_type} ${walltime_epcy} ${mem_epcy}
#  done
#done


num_proc="40"
data_type="sc"
path_design="10X_FACS"
data_project="10X_FACS"
walltime_epcy="200:00:00"
mem_epcy="200Gb"
designs_10X="10085_b_cells 10224_memory_t 10479_naive_t 11953_naive_cytotoxic 8385_cd56_nk 10209_cytotoxic_t 10263_regulatory_t	11213_cd4_t 2612_cd14 9232_cd34"
for src_data in cellranger
do
for design in ${designs_10X}
  do
    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design} ${data_project} ${src_data} ${type_exec} ${num_proc} ${data_type} ${walltime_epcy} ${mem_epcy}
  done
done

num_proc="20"
data_type="sc"
path_design="10X_FACS_reduce"
data_project="10X_FACS_reduce"
walltime_epcy="100:00:00"
mem_epcy="120Gb"
designs_10X="1016_cytotoxic_t 1058_b_cells 1077_naive_t 1093_memory_t 1108_regulatory_t 1217_cd4_t 1338_naive_cytotoxic 287_cd14 858_cd56_nk 948_cd34"
for src_data in cellranger
do
for design in ${designs_10X}
  do
    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design} ${data_project} ${src_data} ${type_exec} ${num_proc} ${data_type} ${walltime_epcy} ${mem_epcy}
  done
done

designs_random="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"
path_design="10X_FACS_reduce_random"
data_project="10X_FACS_reduce"
#for src_data in cellranger
#do
#  for design in ${designs_random}
#  do
#    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design} ${data_project} ${src_data} ${type_exec} ${num_proc} ${data_type} ${walltime_epcy} ${mem_epcy}
#  done
#done

num_proc="20"
data_type="sc"
path_design="STAG2_granulo"
data_project="STAG2_granulo"
walltime_epcy="100:00:00"
mem_epcy="120Gb"
designs_STAG2="2662_ko"
for src_data in SEQC
do
  for design in ${designs_STAG2}
  do
    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design} ${data_project} ${src_data} ${type_exec} ${num_proc} ${data_type} ${walltime_epcy} ${mem_epcy}
  done
done

designs_random="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"
path_design="STAG2_granulo_random"
data_project="STAG2_granulo"
#for src_data in SEQC
#do
#  for design in ${designs_random}
#  do
#    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design} ${data_project} ${src_data} ${type_exec} ${num_proc} ${data_type} ${walltime_epcy} ${mem_epcy}
#  done
#done

data_type="sc"
path_design="STAG2_ko_all"
data_project="STAG2_ko_all"
walltime_epcy="100:00:00"
mem_epcy="120Gb"
designs_STAG2="10628_ko"
#for src_data in SEQC
#do
#  for design in ${designs_STAG2}
#  do
#    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design} ${data_project} ${src_data} ${type_exec} ${num_proc} ${data_type} ${walltime_epcy} ${mem_epcy}
#  done
#done

designs_random="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"
path_design="STAG2_ko_all_random"
data_project="STAG2_ko_all"
#for src_data in SEQC
#do
#  for design in ${designs_random}
#  do
#    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design} ${data_project} ${src_data} ${type_exec} ${num_proc} ${data_type} ${walltime_epcy} ${mem_epcy}
#  done
#done
