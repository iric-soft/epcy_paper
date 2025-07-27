
########################################################
# RUN DEG and PEG analysis on leucegene3 and 10X_FACS
########################################################

type_exec="slurm" # "slurm" "bash" "torque"
working_dir="[path_to_epcy_paper_folder]" # e.g. /home/user/epcy_paper
#working_dir="./"

#########################################################
#########################################################
# Run analyses made on leucegene3
#########################################################
#########################################################

num_proc="4"
data_type="bulk"
designs_leucegene="30_inv16 30_t15_17"
path_design_leucegene="leucegene3"
data_project="leucegene3"
#EPCY
walltime_epcy="24:00:00"
mem_epcy="8Gb"
vmem_epcy="24Gb"
#DESeq2
walltime_deseq="4:00:00"
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


###############################################################
###############################################################
# Run single-cell analyses using 10X FACS dataset (~95k cells)
# reduced on average over 20 replicates to 3k, 5k,8k and 10k 
# cells (using ./src/script/other/gen_10X_reduce.r)
###############################################################
###############################################################

num_proc="20"
data_type="sc"
path_design="10X_FACS_reduce"
data_project="10X_FACS_reduce"
#EPCY
walltime_epcy="12:00:00"
mem_epcy="48Gb"
vmem_epcy="124Gb"
#MAST
mem_mast="20Gb"
vmem_mast="24Gb"
walltime_mast="4:00:00"
#limma trend
mem_limma="12Gb"
vmem_limma="16Gb"
walltime_limma="0:30:00"

num_samples="3000 5000 8000 10000"
for num_sample in ${num_samples}
do
  for num_project in {1..20}
  do
    path_design="10X_FACS_reduce_${num_sample}_${num_project}"
    data_project="10X_FACS_reduce_${num_sample}_${num_project}"

    for fullpath in ./data/design/${path_design}/*
    do
      designs_10X=$(basename "$fullpath")
      for src_data in cellranger
      do
        for design in ${designs_10X}
        do
          bash ./src/script/run/DEG_analysis.sh ${working_dir} ${design} ${path_design} ${data_project} ${src_data} ${type_exec} ${num_proc} ${data_type} ${walltime_epcy} ${mem_epcy} ${vmem_epcy} ${walltime_mast} ${mem_mast} ${vmem_mast} ${walltime_limma} ${mem_limma} ${vmem_limma}
        done
      done
    done
  done
done
