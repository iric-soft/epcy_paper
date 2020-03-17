

type_exec="torque" # "slurm" "bash" "torque"
num_proc="4"

num_fold_leucegene="0" # "loo" "10" "0"() ...
path_design_leucegene="leucegene"
designs_leucegene="28_inv16_vs_28 28_inv16 33_MLL 132_FLT3-ITD 139_NPM1_mut 62_CK 9_EVI1 62_Inter 13_Mono5 126_Normal 30_t15_17 18_t8_21 13_Tri8 5_NUP98NSD1"
designs_leucegene_cluster="cluster0 cluster1 cluster2 cluster3 cluster4 cluster5 cluster6 cluster7 cluster8 cluster9 cluster10 cluster11 cluster12"

num_fold_brca="0" # "loo" "10" "0"() ...
path_design_brca="../../data/design/TCGA_BRCA"
designs_brca="104_triple_neg"

######################################
# RUN DEG and PEG analysis
######################################

working_dir="/u/eaudemard/project/epcy_paper/"
data_project="leucegene"
#for src_data in STAR #kallisto
#do
#  #for design in ${designs_leucegene}
#  #do
#  #  bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${data_project} ${src_data} ${type_exec} ${num_fold_leucegene} ${num_proc}
#  #done
#  for design in ${designs_leucegene_cluster}
#  do
#    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design_leucegene} ${data_project} ${src_data} ${type_exec} ${num_fold_leucegene} ${num_proc}
#  done
#done

path_design_leucegene="leucegene_v2"
designs_leucegene="31_inv16 48_MLL 92_CK 13_EVI1 87_Inter 23_Mono5 275_NK 30_t15_17 18_t8_21 18_Tri8 6_NUP98NSD1"
working_dir="/u/eaudemard/project/epcy_paper/"
data_project="leucegene_v2"
for src_data in STAR #kallisto
do
  for design in ${designs_leucegene}
  do
    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design_leucegene} ${data_project} ${src_data} ${type_exec} ${num_fold_leucegene} ${num_proc}
  done
done

#path_design_leucegene="leucegene_jf"
#designs_leucegene="31_inv16 53_MLL 90_CK 13_EVI1 80_Inter 22_Mono5 267_NK 30_t15_17 18_t8_21 18_Tri8 19_NUP98"
#designs_leucegene_cluster="cluster0 cluster1 cluster2 cluster3 cluster4 cluster5 cluster6 cluster7 cluster8 cluster9 cluster10 cluster11 cluster12 cluster13 cluster14"
#working_dir="/u/eaudemard/project/epcy_paper/"
#data_project="leucegene_v2"
#for src_data in STAR #kallisto
#do
#  for design in ${designs_leucegene}
#  do
#    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design_leucegene} ${data_project} ${src_data} ${type_exec} ${num_fold_leucegene} ${num_proc}
#  done
#  for design in ${designs_leucegene_cluster}
#  do
#    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design_leucegene} ${data_project} ${src_data} ${type_exec} ${num_fold_leucegene} ${num_proc}
#  done
#done


path_design_leucegene="leucegene_comp_tt/comp"
designs_leucegene="inv16 MLL CK EVI1 Inter Mono5 NK t15_17 t8_21 Tri8 NUP98NSD1"
working_dir="/u/eaudemard/project/epcy_paper/"
data_project="leucegene_v2"
#for src_data in STAR #kallisto
#do
#  for design in ${designs_leucegene}
#  do
#    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design_leucegene} ${data_project} ${src_data} ${type_exec} ${num_fold_leucegene} ${num_proc}
#  done
#done


path_design_leucegene="leucegene_comp_tt/test_train"
designs_leucegene="inv16 MLL CK EVI1 Inter Mono5 NK t15_17 t8_21 Tri8 NUP98NSD1"
working_dir="/u/eaudemard/project/epcy_paper/"
data_project="leucegene_v2"
#for src_data in STAR #kallisto
#do
#  for design in ${designs_leucegene}
#  do
#    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design_leucegene} ${data_project} ${src_data} ${type_exec} ${num_fold_leucegene} ${num_proc}
#  done
#done



path_design_leucegene_tmp="leucegene_comp_tt_v2"
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
do
  for dataset in comp test_train
  do
    path_design_leucegene="${path_design_leucegene_tmp}/${i}/${dataset}"
    echo $path_design_leucegene
    designs_leucegene="inv16 MLL CK EVI1 Inter Mono5 NK t15_17 t8_21 Tri8 NUP98NSD1"
    working_dir="/u/eaudemard/project/epcy_paper/"
    data_project="leucegene_v2"
    for src_data in STAR #kallisto
    do
      for design in ${designs_leucegene}
      do
        bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design_leucegene} ${data_project} ${src_data} ${type_exec} ${num_fold_leucegene} ${num_proc}
      done
    done
  done
done




path_design_leucegene="TARGET_AML"
designs_leucegene="24_inv16 21_MLL 28_Other 27_Normal 17_t8_21 8_Unknown"
working_dir="/u/eaudemard/project/epcy_paper/"
data_project="TARGET_AML"
for src_data in htseq
do
  for design in ${designs_leucegene}
  do
    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design_leucegene} ${data_project} ${src_data} ${type_exec} ${num_fold_leucegene} ${num_proc}
  done
done

#num_fold_brca="0" # "loo" "10" "0"() ...
#path_design_brca="TCGA_BRCA"
#designs_brca="104_triple_neg"
#data_project="TCGA_BRCA"
#for src_data in htseq
#do
#  for design in ${designs_brca}
#  do
#    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design_brca}  ${data_project} ${src_data} ${type_exec} ${num_fold_brca} ${num_proc}
#  done
#done
