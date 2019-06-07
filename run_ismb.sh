
num_fold="loo" # "loo" "10" ...
type_exec="torque" # "slurm" "bash" "torque"
subgroups="28_inv16_vs_28" #"28_inv16 28_inv16_vs_28 33_MLL 132_FLT3-ITD 139_NPM1_mut"

######################################
# Generate cross-validation folder
######################################

cd src/pyres

for subgroup in ${subgroups}
do
  if [ ${num_fold} == "loo" ]
  then # for leave one out cross-validation
    python3 -m pyres gen_cv -p ../../data/design/${subgroup} --loo
  else # For x-fold cross-validation
    python3 -m pyres gen_cv -p ../../data/design/${subgroup} -f ${num_fold}
  fi
done

cd ../..

######################################
# RUN DEG and PEG analysis
######################################

working_dir="/u/eaudemard/project/epcy_paper/"
data_project="leucegene"
for src_data in STAR #STAR_RSEM kallisto
do
  for subgroup in ${subgroups}
  do
    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${subgroup} ${data_project} ${src_data} ${type_exec} ${num_fold}
  done
done
