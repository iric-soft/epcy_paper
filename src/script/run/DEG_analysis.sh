#DEG_analysis.sh [working_dir] [subgroup] [dir_design] [data_project] [src_data] [type_exec] [num_fold] [num_proc]
#[subgroup]: 28_inv16 33_MLL ...
#[dir_design]: leucegene leucegene_v2 TCGA_BRCA ...
#[data_project]: leucegene ...
#[src_data]: STAR ...
#[type_exec]: bash torque ...
#[num_fold]: 10 loo ...

working_dir=$1
subgroup=$2
dir_design=$3
data_project=$4
src_data=$5
type_exec=$6
num_proc=$7
data_type=$8
path_data="${working_dir}/data"
path_jobout="${working_dir}/tmp"

###############################################################################
# wall time and memory usage for num_proc="4"
###############################################################################
# for bulk
walltime_epcy="24:00:00"
mem_epcy="8Gb"

#for sc
#walltime_epcy="200:00:00"
#mem_epcy="200Gb"

walltime_edger="3:00:00"
mem_edger="3Gb"

walltime_limma="12:00:00"
mem_limma="120Gb"

walltime_mast="12:00:00"
mem_mast="212Gb"

walltime_LDE="12:00:00"
mem_LDE="80Gb" # 2 * num_proc


echo "##################### $subgroup #########################"

path_cmd="${working_dir}/src/script/cmd"
path_matrix="${path_data}/${data_project}/${src_data}/"
path_anno="${path_data}/other/Homo_sapiens.GRCh38.84.gff3"

send2torque()
{
	eval cmd="$1"
	eval ppn="$2"
	eval mem="$3"
	eval walltime="$4"
	eval job_name="$5"
	eval path_jobout="$6"

	mkdir -p ${path_jobout}
	echo "${cmd}" | qsub -V -l nodes=1:ppn=${ppn},mem=${mem},walltime=${walltime} -j oe -N ${job_name} -o ${path_jobout}
	#FORDEBUG
	#echo "${cmd} | qsub -V -l nodes=1:ppn=${ppn},mem=${mem},walltime=${walltime} -j oe -N ${job_name} -o ${path_jobout}"
}

exec_cmd()
{
	eval type_exec="$1"
	eval cmd="$2"

	if [ $type_exec == "torque" ]
	then
		eval n_proc="$3"
		eval mem="$4"
		eval walltime="$5"
		eval job_name="$6"
		eval path_jobout="$7"
	  send2torque "\${cmd}" ${n_proc} ${mem} ${walltime} ${job_name} ${path_jobout}
	fi

	if [ $type_exec == "bash" ]
	then
	  ${cmd}
		#FORDEBUG
		#echo ${cmd}
	fi

}

epcy()
{
  eval type_run="$1"
  eval input_design="$2"
	eval path_output="$3"
	eval path_jobout="$4"

	if [ $type_run == "quant" ] && [ -f ${path_matrix}/quant.xls ]
	then
		if [ ! -f ${path_output}/quant/predictive_capability.xls ]
		then
			job_name="epcy_${subgroup}"
			path_output_epcy="${path_output}/quant"
			cmd=$(bash ${path_cmd}/epcy.sh ${input_design} ${path_matrix} ${path_output_epcy} ${num_proc})
			exec_cmd ${type_exec} "\${cmd}" ${num_proc} ${mem_epcy} ${walltime_epcy} ${job_name} ${path_jobout}
		else
			echo "epcy ${subgroup} done!"
		fi
	fi

	if [ $type_run == "count" ] && [ -f ${path_matrix}/readcounts.xls ]
	then
		if [ ! -f ${path_output}/readcounts/predictive_capability.xls ]
		then
			job_name="epcy_count_${subgroup}"
			path_output_epcy="${path_output}/readcounts"
			cmd=$(bash ${path_cmd}/epcy_count.sh ${input_design} ${path_matrix} ${path_output_epcy} ${num_proc})
			exec_cmd ${type_exec} "\${cmd}" ${num_proc} ${mem_epcy} ${walltime_epcy} ${job_name} ${path_jobout}
		else
			echo "epcy count ${subgroup} done!"
		fi
	fi

	if [ $type_run == "count_sc" ] && [ -f ${path_matrix}/readcounts.xls ]
	then
		if [ ! -f ${path_output}/readcounts/predictive_capability.xls ]
		then
			job_name="epcy_sc_count_${subgroup}"
			path_output_epcy="${path_output}/readcounts"
			cmd=$(bash ${path_cmd}/epcy_sc_count.sh ${input_design} ${path_matrix} ${path_output_epcy} ${num_proc})
			exec_cmd ${type_exec} "\${cmd}" ${num_proc} ${mem_epcy} ${walltime_epcy} ${job_name} ${path_jobout}
		else
			echo "epcy sc count ${subgroup} done!"
		fi
	fi

	if [ $type_run == "bagging" ] && [ -f ${path_matrix}/readcounts.xls ]
	then
		if [ ! -f ${path_output}/readcounts_bagging/predictive_capability.xls ]
		then
			job_name="epcy_bagging_${subgroup}"
			path_output_epcy="${path_output}/readcounts_bagging"
			cmd=$(bash ${path_cmd}/epcy_bagging.sh ${input_design} ${path_matrix} ${path_output_epcy} ${num_proc})
			exec_cmd ${type_exec} "\${cmd}" ${num_proc} ${mem_epcy} ${walltime_epcy} ${job_name} ${path_jobout}
		else
			echo "epcy bagging ${subgroup} done!"
		fi
	fi
}

deseq()
{
  eval input_design="$1"
	eval path_output="$2"
	eval path_jobout="$3"
	if [ -f ${path_matrix}/readcounts.xls ]
	then
		if [ ! -f ${path_output}/readcounts/deseq2_genes.xls ]
		then
		  job_name="DESEQ2_${subgroup}"
			path_output_deseq2="${path_output}/readcounts"
		  path_exec="${working_dir}/src/script/exec/deseq2.r"
		  cmd=$(bash ${path_cmd}/deseq2.sh ${path_exec} ${input_design} ${path_matrix} ${path_output_deseq2} ${num_proc})
			exec_cmd ${type_exec} "\${cmd}" ${num_proc} ${mem_LDE} ${walltime_LDE} ${job_name} ${path_jobout}
		else
			echo "deseq2 ${subgroup} done!"
		fi
	fi
}

edger()
{
  eval input_design="$1"
	eval path_output="$2"
	eval path_jobout="$3"

	if [ -f ${path_matrix}/readcounts.xls ]
	then
		if [ ! -f ${path_output}/readcounts/edger_genes.xls ]
		then
		  job_name="EDGER_${subgroup}"
			path_output_edger="${path_output}/readcounts"
		  path_exec="${working_dir}/src/script/exec/edger.r"
		  cmd=$(bash ${path_cmd}/edger.sh ${path_exec} ${input_design} ${path_matrix} ${path_output_edger})
			exec_cmd ${type_exec} "\${cmd}" 1 ${mem_edger} ${walltime_edger} ${job_name} ${path_jobout}
		else
			echo "edger ${subgroup} done!"
		fi
	fi
}

limma_trend()
{
  eval input_design="$1"
	eval path_output="$2"
	eval path_jobout="$3"

	if [ -f ${path_matrix}/readcounts.xls ]
	then
		if [ ! -f ${path_output}/readcounts/limma_trend_genes.xls ]
		then
		  job_name="trend_${subgroup}"
			path_output_limma="${path_output}/readcounts"
		  path_exec="${working_dir}/src/script/exec/limma_trend.r"
		  cmd=$(bash ${path_cmd}/limma.sh ${path_exec} ${input_design} ${path_matrix} ${path_output_limma})
			exec_cmd ${type_exec} "\${cmd}" 1 ${mem_limma} ${walltime_limma} ${job_name} ${path_jobout}
		else
			echo "limma trend ${subgroup} done!"
		fi
	fi
}

limma()
{
  eval input_design="$1"
	eval path_output="$2"
	eval path_jobout="$3"

	if [ -f ${path_matrix}/readcounts.xls ]
	then
		if [ ! -f ${path_output}/readcounts/limma_voom_genes.xls ]
		then
		  job_name="voom_${subgroup}"
			path_output_limma="${path_output}/readcounts"
		  path_exec="${working_dir}/src/script/exec/limma.r"
		  cmd=$(bash ${path_cmd}/limma.sh ${path_exec} ${input_design} ${path_matrix} ${path_output_limma})
			exec_cmd ${type_exec} "\${cmd}" 1 ${mem_limma} ${walltime_limma} ${job_name} ${path_jobout}
		else
			echo "limma voom ${subgroup} done!"
		fi
	fi
}

mast()
{
  eval input_design="$1"
	eval path_output="$2"
	eval path_jobout="$3"

	if [ -f ${path_matrix}/readcounts.xls ]
	then
		if [ ! -f ${path_output}/readcounts/mast_genes.xls ]
		then
		  job_name="mast_${subgroup}"
			path_output_mast="${path_output}/readcounts"
		  path_exec="${working_dir}/src/script/exec/mast.r"
		  cmd=$(bash ${path_cmd}/mast.sh ${path_exec} ${input_design} ${path_matrix} ${path_output_mast})
			exec_cmd ${type_exec} "\${cmd}" 1 ${mem_mast} ${walltime_mast} ${job_name} ${path_jobout}
		else
			echo "mast ${subgroup} done!"
		fi
	fi
}


LDE()
{
  eval input_design="$1"
	eval path_output="$2"
	eval path_jobout="$3"

	if [ -f ${path_matrix}/readcounts.xls ]
	then
		if [ ! -f ${path_output}/readcounts/deseq2_genes.xls ] && [ ! -f ${path_output}/readcounts/limma_voom_genes.xls ] && [ ! -f ${path_output}/readcounts/edger_genes.xls ]
		then
		  job_name="LDE_${subgroup}"
			path_output_lde="${path_output}/readcounts"
		  path_exec="${working_dir}/src/script/exec/LDE.r"
		  cmd=$(bash ${path_cmd}/LDE.sh ${path_exec} ${input_design} ${path_matrix} ${path_output_lde} ${num_proc})
			exec_cmd ${type_exec} "\${cmd}" ${num_proc} ${mem_LDE} ${walltime_LDE} ${job_name} ${path_jobout}
		else
			if [ ! -f ${path_output}/readcounts/limma_voom_genes.xls ]
			then
				limma ${input_design} ${path_output} ${path_jobout}
			else
				echo "limma ${subgroup} done!"
			fi

			if [ ! -f ${path_output}/readcounts/deseq2_genes.xls ]
			then
				deseq ${input_design} ${path_output} ${path_jobout}
			else
				echo "deseq2 ${subgroup} done!"
			fi


			if [ ! -f ${path_output}/readcounts/edger_genes.xls ]
			then
				edger ${input_design} ${path_output} ${path_jobout}
			else
				echo "edger ${subgroup} done!"
			fi
		fi
	fi
}

path_design="${path_data}/design/${dir_design}/${subgroup}"
path_output="${path_design}/${src_data}"
path_jobout_subgroup="${path_jobout}/${data_project}/${src_data}/${subgroup}/"

if [ $data_type == "bulk" ]
then
	#epcy "quant" ${path_design} ${path_output} ${path_jobout_subgroup}
	epcy "count" ${path_design} ${path_output} ${path_jobout_subgroup}
	epcy "bagging" ${path_design} ${path_output} ${path_jobout_subgroup}
	LDE ${path_design} ${path_output} ${path_jobout_subgroup}

fi

if [ $data_type == "sc" ]
then
  epcy "count_sc" ${path_design} ${path_output} ${path_jobout_subgroup}
	limma_trend ${path_design} ${path_output} ${path_jobout_subgroup}
	mast ${path_design} ${path_output} ${path_jobout_subgroup}
fi
