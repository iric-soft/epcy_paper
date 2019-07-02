
num_fold="loo" # "loo" "10" "0"() ...
type_exec="torque" # "slurm" "bash" "torque"
#designs="28_inv16"
designs="28_inv16_vs_28 28_inv16 33_MLL 132_FLT3-ITD 139_NPM1_mut"
num_proc="4"

######################################
# Generate cross-validation folder
######################################

cd src/pyres

create_cv()
{
  if [ ${num_fold} == "loo" ]
  then # for leave one out cross-validation
    python3 -m pyres gen_cv -p ../../data/design/${subgroup} --loo
  else # For x-fold cross-validation
    python3 -m pyres gen_cv -p ../../data/design/${subgroup} -f ${num_fold}
  fi
}

for design in ${designs}
do
  if [ -d ../../data/design/${design}/cv/fold${num_fold} ]
  then
    printf "This cross-validation directory already exist: ../../data/design/${design}/cv/fold${num_fold}\n"
    read -p "Do you want delete it ? [Y/N] " -n 1 -r
    if [[ $REPLY =~ ^[Yy]$ ]]
    then
      printf "\n-> Delete and create a new one ...\n"
      rm -fr ../../data/design/${design}/cv/fold${num_fold}
      create_cv
    else
      printf "\n-> Do nothing and move to next subgroup ...\n"
    fi
  else
    printf "\n-> Create cross-validation directory for: ${design} ...\n"
    create_cv
  fi
done

cd ../..

######################################
# RUN DEG and PEG analysis
######################################

working_dir="/u/eaudemard/project/epcy_paper/"
data_project="leucegene"
for src_data in STAR #kallisto
do
  for design in ${designs}
  do
    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${data_project} ${src_data} ${type_exec} ${num_fold} ${num_proc}
  done
done
