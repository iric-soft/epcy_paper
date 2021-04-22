
########################################################
# RUN DEG and PEG analysis on leucegene3 and 10X_FACS
########################################################

type_exec="torque" # "slurm" "bash"
num_proc="4"
working_dir="/u/eaudemard/project/epcy_paper/"
#working_dir="./"


#########################################################
# Run on small leucegene
#
data_type="bulk"
designs_leucegene="inv16_vs_t15_17"
path_design_leucegene="small_leucegene"
data_project="small_leucegene"
#EPCY
walltime_epcy="24:00:00"
mem_epcy="8Gb"
vmem_epcy="24Gb"
#deseq2
walltime_deseq="2:00:00"
mem_deseq="24Gb"
vmem_deseq="48Gb"
#limma voom and edger
walltime_limma="4:00:00"
mem_limma="6Gb"
vmem_limma="24Gb"
for src_data in STAR_RSEM
do
  for design in ${designs_leucegene}
  do
    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design_leucegene} ${data_project} ${src_data} ${type_exec} ${num_proc} ${data_type} ${walltime_epcy} ${mem_epcy} ${vmem_epcy} ${walltime_deseq} ${mem_deseq} ${vmem_deseq} ${walltime_limma} ${mem_limma} ${vmem_limma}
  done
done

designs_leucegene="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"
path_design_leucegene="small_leucegene_random"
data_project="small_leucegene"
for src_data in STAR_RSEM
do
  for design in ${designs_leucegene}
  do
    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design_leucegene} ${data_project} ${src_data} ${type_exec} ${num_proc} ${data_type} ${walltime_epcy} ${mem_epcy} ${vmem_epcy} ${walltime_deseq} ${mem_deseq} ${vmem_deseq} ${walltime_limma} ${mem_limma} ${vmem_limma}
  done
done


#EPCY
walltime_epcy="24:00:00"
mem_epcy="50Gb"
vmem_epcy="64Gb"
#DESeq2
walltime_deseq="2:00:00"
mem_deseq="80Gb"
vmem_deseq="90Gb"
#limma voom and edger
walltime_limma="4:00:00"
mem_limma="6Gb"
vmem_limma="24Gb"
designs_leucegene="inv16_vs_t15_17"
p_subs="0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9"
nums="1 2 3 4 5 6 7 8 9 10"
path_design_leucegene="small_leucegene_rep_ss"
data_project="small_leucegene_rep"
for src_data in STAR_RSEM
do
  for design in ${designs_leucegene}
  do
    for p in ${p_subs}
    do
      for num in ${nums}
      do
        design_ss="${design}/${p}/${num}"
        bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design_ss} ${path_design_leucegene} ${data_project} ${src_data} ${type_exec} ${num_proc} ${data_type} ${walltime_epcy} ${mem_epcy} ${vmem_epcy} ${walltime_deseq} ${mem_deseq} ${vmem_deseq} ${walltime_limma} ${mem_limma} ${vmem_limma}
      done
    done
  done
done
