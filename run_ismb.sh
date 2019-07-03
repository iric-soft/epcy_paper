

type_exec="torque" # "slurm" "bash" "torque"
num_proc="4"

num_fold_leucegene="loo" # "loo" "10" "0"() ...
path_design_leucegene="../../data/design/leucegene"
designs_leucegene="28_inv16_vs_28 28_inv16 33_MLL 132_FLT3-ITD 139_NPM1_mut"

num_fold_brca="200" # "loo" "10" "0"() ...
path_design_brca="../../data/design/TCGA_BRCA"
designs_brca="104_triple_neg"

######################################
# Generate cross-validation folder
######################################
cd src/pyres
create_cv()
{
  if [ ${num_fold} == "loo" ]
  then # for leave one out cross-validation
    python3 -m pyres gen_cv -p ${path_design}/${design} --loo
  else # For x-fold cross-validation
    python3 -m pyres gen_cv -p ${path_design}/${design} -f ${num_fold}
  fi
}

create_all_cv()
{
  for design in ${designs}
  do
    if [ -d ${path_design}/${design}/cv/fold${num_fold} ]
    then
      printf "This cross-validation directory already exist: ${path_design}/${design}/cv/fold${num_fold}\n"
      read -p "Do you want delete it ? [Y/N] " -n 1 -r
      if [[ $REPLY =~ ^[Yy]$ ]]
      then
        printf "\n-> Delete and create a new one ...\n"
        rm -fr ${path_design}/${design}/cv/fold${num_fold}
        create_cv
      else
        printf "\n-> Do nothing and move to next subgroup ...\n"
      fi
    else
      printf "\n-> Create cross-validation directory for: ${design} ...\n"
      create_cv
    fi
  done
}
################### For leucegene #############################################
#designs="28_inv16"
num_fold=${num_fold_leucegene}
path_design=${path_design_leucegene}
designs=${designs_leucegene}
create_all_cv

################### For TCGA BRCA #############################################
num_fold=${num_fold_brca}
path_design=${path_design_brca}
designs=${designs_brca}
create_all_cv

cd ../..

######################################
# RUN DEG and PEG analysis
######################################

working_dir="/u/eaudemard/project/epcy_paper/"
data_project="leucegene"
for src_data in STAR #kallisto
do
  for design in ${designs_leucegene}
  do
    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${data_project} ${src_data} ${type_exec} ${num_fold_leucegene} ${num_proc}
  done
done

data_project="TCGA_BRCA"
for src_data in htseq
do
  for design in ${designs_brca}
  do
    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${data_project} ${src_data} ${type_exec} ${num_fold_brca} ${num_proc}
  done
done
