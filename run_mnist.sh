

type_exec="torque" # "slurm" "bash" "torque"
num_proc="4"

num_fold_mnist="0"
path_design_mnist="mnist_fashion_train"
designs_mnist="6000_Ankle_boot 6000_Bag 6000_Coat 6000_Dress 6000_Pullover 6000_Sandal 6000_Shirt 6000_Sneaker 6000_T-shirt 6000_Trouser"
working_dir="/u/eaudemard/project/epcy_paper/"
data_project="mnist_fashion"
for src_data in train
do
  for design in ${designs_mnist}
  do
    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design_mnist} ${data_project} ${src_data} ${type_exec} ${num_fold_mnist} ${num_proc}
  done
done
