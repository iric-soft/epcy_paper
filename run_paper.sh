

type_exec="torque" # "slurm" "bash" "torque"
num_proc="4"

num_fold_leucegene="0" # "loo" "10" "0"() ...
path_design_leucegene="../../data/design/leucegene"
designs_leucegene="28_inv16_vs_28 28_inv16 33_MLL 132_FLT3-ITD 139_NPM1_mut 5_NUP98NSD1 9_EVI1 13_Mono5 13_Tri8 18_t8_21 30_t15_17 62_CK 62_Inter 27_inv16 29_inv16"

######################################
# RUN DEG and PEG analysis
######################################

working_dir="/u/eaudemard/project/epcy_paper/"
path_design_leucegene="leucegene"
data_project="leucegene"
for src_data in STAR #kallisto
do
  for design in ${designs_leucegene}
  do
    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design_leucegene} ${data_project} ${src_data} ${type_exec} ${num_fold_leucegene} ${num_proc}
  done
done


designs_leucegene="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"
working_dir="/u/eaudemard/project/epcy_paper/"
path_design_leucegene="leucegene_random"
data_project="leucegene"
for src_data in STAR #kallisto
do
  for design in ${designs_leucegene}
  do
    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design_leucegene} ${data_project} ${src_data} ${type_exec} ${num_fold_leucegene} ${num_proc}
  done
done
