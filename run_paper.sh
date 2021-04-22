
########################################################
# RUN DEG and PEG analysis on leucegene3 and 10X_FACS
########################################################

type_exec="torque" # "slurm" "bash"
num_proc="4"
working_dir="/u/eaudemard/project/epcy_paper/"
#working_dir="./"


#########################################################
# Run on leucegene3
#
data_type="bulk"
designs_leucegene="30_inv16 30_t15_17"
path_design_leucegene="leucegene3"
data_project="leucegene3"
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
path_design_leucegene="leucegene3_random"
data_project="leucegene3"
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
designs_leucegene="30_t15_17"
p_subs="0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9"
nums="1 2 3 4 5 6 7 8 9 10"
path_design_leucegene="leucegene3_rep_ss"
data_project="leucegene3_rep"
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




###################################
# Run on 10X FACS reduced version
#
num_proc="20"
data_type="sc"
path_design="10X_FACS_reduce"
data_project="10X_FACS_reduce"
#EPCY
walltime_epcy="4:00:00"
mem_epcy="48Gb"
vmem_epcy="124Gb"
#MAST
mem_mast="20Gb"
vmem_mast="24Gb"
walltime_mast="1:00:00"
#limma trend
mem_limma="12Gb"
vmem_limma="16Gb"
walltime_limma="0:30:00"

designs_10X="1016_cytotoxic_t 1058_b_cells 1077_naive_t 1093_memory_t 1108_regulatory_t 1217_cd4_t 1338_naive_cytotoxic 287_cd14 858_cd56_nk 948_cd34"
for src_data in cellranger
do
for design in ${designs_10X}
  do
    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design} ${data_project} ${src_data} ${type_exec} ${num_proc} ${data_type} ${walltime_epcy} ${mem_epcy} ${vmem_epcy} ${walltime_mast} ${mem_mast} ${vmem_mast} ${walltime_limma} ${mem_limma} ${vmem_limma}
  done
done


designs_random="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"
path_design="10X_FACS_reduce_random"
data_project="10X_FACS_reduce"
for src_data in cellranger
do
  for design in ${designs_random}
  do
    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design} ${data_project} ${src_data} ${type_exec} ${num_proc} ${data_type} ${walltime_epcy} ${mem_epcy} ${vmem_epcy} ${walltime_mast} ${mem_mast} ${vmem_mast} ${walltime_limma} ${mem_limma} ${vmem_limma}
  done
done


################################
# Run on 10X FACS
# this dataset have arround 95k cells (samples) and EPCY need a lot of
# ressources to analysed them: 40 cpu, 124Go mem, 200Go vmem and
# a minimum of 5 hours
#
#num_proc="40"
#data_type="sc"
#path_design="10X_FACS"
#data_project="10X_FACS"
##EPCY
#walltime_epcy="200:00:00"
#mem_epcy="200Gb"
#vmem_epcy="248Gb"
##MAST
#mem_mast="20Gb"
#vmem_mast="48Gb"
#walltime_mast="4:00:00"
##limma trend
#mem_limma="10Gb"
#vmem_limma="20Gb"
#walltime_limma="1:00:00"
#designs_10X="10085_b_cells 10224_memory_t 10479_naive_t 11953_naive_cytotoxic 8385_cd56_nk 10209_cytotoxic_t 10263_regulatory_t	11213_cd4_t 2612_cd14 9232_cd34"
#for src_data in cellranger
#do
#for design in ${designs_10X}
#  do
#    bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design} ${data_project} ${src_data} ${type_exec} ${num_proc} ${data_type} ${walltime_epcy} ${mem_epcy}  ${vmem_epcy} ${walltime_mast} ${mem_mast} ${vmem_mast} ${walltime_limma} ${mem_limma} ${vmem_limma}
#  done
#done
