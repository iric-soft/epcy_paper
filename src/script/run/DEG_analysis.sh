#DEG_analysis.sh [working_dir] [subgroup] [data_project] [src_data] [type_exec] [num_fold]
#[subgroup]: 28_inv16 33_MLL ...
#[data_project]: leucegene ...
#[src_data]: STAR ...
#[type_exec]: bash torque ...
#[num_fold]: 10 loo ...

working_dir=$1
subgroup=$2
data_project=$3
src_data=$4
type_exec=$5
num_fold=$6
path_data="${working_dir}/data"
path_jobout="${working_dir}/tmp"

num_proc="6"

walltime_epcy="12:00:00"
mem_epcy="8Gb"

walltime_edger="12:00:00"
mem_edger="2Gb"

walltime_limma="12:00:00"
mem_limma="2Gb"

walltime_deseq2="12:00:00"
mem_deseq2="12Gb" # 2 * num_proc

path_cmd="${working_dir}/src/script/cmd"
path_design="${path_data}/design/${project_name}/${subgroup}"
path_output="${path_design}/${src_data}"
path_data="${path_data}/${data_project}/${src_data}/${project_name}"
path_jobout_subgroup="${path_jobout}/${data_project}/${src_data}/${subgroup}/"
mkdir -p ${path_jobout_subgroup}

echo $path_design

send2torque()
{
	eval cmd="$1"
	eval num_proc="$2"
	eval mem="$3"
	eval walltime="$4"
	eval job_name="$5"
	eval path_jobout="$6"

	mkdir -p ${path_jobout}
	echo "${cmd}" | qsub -V -l nodes=1:ppn=${num_proc},mem=${mem},walltime=${walltime} -j oe -N ${job_name} -o ${path_jobout}
	#FORDEBUG
	#echo "${cmd} | qsub -V -l nodes=1:ppn=${num_proc},mem=${mem},walltime=${walltime} -j oe -N ${job_name} -o ${path_jobout}"
}

exec_cmd()
{
	eval type_exec="$1"
	eval cmd="$2"

	if [ $type_exec == "torque" ]
	then
		eval num_proc="$3"
		eval mem="$4"
		eval walltime="$5"
		eval job_name="$6"
		eval path_jobout="$7"
	  send2torque "\${cmd}" ${num_proc} ${mem} ${walltime} ${job_name} ${path_jobout}
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
	eval path_output="$2"
	eval subgroup="$3"
	eval path_design="$4"
	eval path_cmd="$5"
	eval path_data="$6"
	eval num_proc="$7"
	eval type_exec="$8"
	eval mem_epcy="$9"
	eval walltime_epcy="${10}"
	eval path_jobout="${11}"

	if [ $type_run == "tpm" ]
	then
		if [ ! -f ${path_output}/tpm/prediction_capability.xls ]
		then
			job_name="epcy_${subgroup}"
			path_output_epcy="${path_output}/tpm"
			cmd=$(bash ${path_cmd}/epcy.sh ${path_design} ${path_data} ${path_output_epcy} ${num_proc})
			exec_cmd ${type_exec} "\${cmd}" ${num_proc} ${mem_epcy} ${walltime_epcy} ${job_name} ${path_jobout}
		else
			echo "epcy ${subgroup} done!"
		fi
	fi

	if [ $type_run == "count" ]
	then
		if [ ! -f ${path_output}/readcounts/prediction_capability.xls ]
		then
			job_name="epcy_count_${subgroup}"
			path_output_epcy="${path_output}/readcounts"
			cmd=$(bash ${path_cmd}/epcy_count.sh ${path_design} ${path_data} ${path_output_epcy} ${num_proc})
			exec_cmd ${type_exec} "\${cmd}" ${num_proc} ${mem_epcy} ${walltime_epcy} ${job_name} ${path_jobout}
		else
			echo "epcy count ${subgroup} done!"
		fi
	fi
}

deseq()
{
	eval path_output="$1"
	eval subgroup="$2"
	eval working_dir="$3"
	eval path_cmd="$4"
	eval path_design="$5"
	eval path_data="$6"
	eval num_proc="$7"
	eval type_exec="$8"
	eval mem_deseq2="$9"
	eval walltime_deseq2="${10}"
	eval path_jobout="${11}"

	if [ ! -f ${path_output}/readcounts/deseq2_genes.xls ]
	then
	  job_name="DESEQ2_${subgroup}"
		path_output_deseq2="${path_output}/readcounts"
	  path_exec="${working_dir}/src/script/exec/deseq2.r"
	  cmd=$(bash ${path_cmd}/deseq2.sh ${path_exec} ${path_design} ${path_data} ${path_output_deseq2} ${num_proc})
		exec_cmd ${type_exec} "\${cmd}" ${num_proc} ${mem_deseq2} ${walltime_deseq2} ${job_name} ${path_jobout}
	else
		echo "deseq2 ${subgroup} done!"
	fi
}

edger()
{
	eval path_output="$1"
	eval subgroup="$2"
	eval working_dir="$3"
	eval path_cmd="$4"
	eval path_design="$5"
	eval path_data="$6"
	eval type_exec="$7"
	eval mem_edger="$8"
	eval walltime_edger="$9"
	eval path_jobout="${10}"

	if [ ! -f ${path_output}/readcounts/edger_genes.xls ]
	then
	  job_name="EDGER_${subgroup}"
		path_output_edger="${path_output}/readcounts"
	  path_exec="${working_dir}/src/script/exec/edger.r"
	  cmd=$(bash ${path_cmd}/edger.sh ${path_exec} ${path_design} ${path_data} ${path_output_edger})
		exec_cmd ${type_exec} "\${cmd}" 1 ${mem_edger} ${walltime_edger} ${job_name} ${path_jobout}
	else
		echo "edger ${subgroup} done!"
	fi
}

limma()
{
	eval path_output="$1"
	eval subgroup="$2"
	eval working_dir="$3"
	eval path_cmd="$4"
	eval path_design="$5"
	eval path_data="$6"
	eval type_exec="$7"
	eval mem_limma="$8"
	eval walltime_limma="$9"
	eval path_jobout="${10}"

	if [ ! -f ${path_output}/readcounts/limma_voom_genes.xls ]
	then
	  job_name="limma_${subgroup}"
		path_output_limma="${path_output}/readcounts"
	  path_exec="${working_dir}/src/script/exec/limma.r"
	  cmd=$(bash ${path_cmd}/limma.sh ${path_exec} ${path_design} ${path_data} ${path_output_limma})
		exec_cmd ${type_exec} "\${cmd}" 1 ${mem_limma} ${walltime_limma} ${job_name} ${path_jobout}
	else
		echo "limma ${subgroup} done!"
	fi
}

epcy "tpm" ${path_output} ${subgroup} ${path_design} ${path_cmd} ${path_data} ${num_proc} ${type_exec} ${mem_epcy} ${walltime_epcy} ${path_jobout_subgroup}
epcy "count" ${path_output} ${subgroup} ${path_design} ${path_cmd} ${path_data} ${num_proc} ${type_exec} ${mem_epcy} ${walltime_epcy} ${path_jobout_subgroup}
deseq ${path_output} ${subgroup} ${working_dir} ${path_cmd} ${path_design} ${path_data} ${num_proc} ${type_exec} ${mem_deseq2} ${walltime_deseq2} ${path_jobout_subgroup}
edger ${path_output} ${subgroup} ${working_dir} ${path_cmd} ${path_design} ${path_data} ${type_exec} ${mem_edger} ${walltime_edger} ${path_jobout_subgroup}
limma ${path_output} ${subgroup} ${working_dir} ${path_cmd} ${path_design} ${path_data} ${type_exec} ${mem_limma} ${walltime_limma} ${path_jobout_subgroup}

if [ ! $num_fold == "0" ]
then
	cv_dir="${path_design}/cv/fold${num_fold}"
	for num_cv_dir in `ls -C1 ${cv_dir}/`
	do
		path_design_cv="${cv_dir}/${num_cv_dir}/train"
		path_output_cv="${path_design_cv}/${src_data}"
		path_jobout_subgroup_cv="${path_jobout_subgroup}/cv/${num_fold}/${num_cv_dir}"

		epcy "tpm" ${path_output_cv} ${subgroup} ${path_design_cv} ${path_cmd} ${path_data} ${num_proc} ${type_exec} ${mem_epcy} ${walltime_epcy} ${path_jobout_subgroup_cv}
		epcy "count" ${path_output_cv} ${subgroup} ${path_design_cv} ${path_cmd} ${path_data} ${num_proc} ${type_exec} ${mem_epcy} ${walltime_epcy} ${path_jobout_subgroup_cv}
		deseq ${path_output_cv} ${subgroup} ${working_dir} ${path_cmd} ${path_design_cv} ${path_data} ${num_proc} ${type_exec} ${mem_deseq2} ${walltime_deseq2} ${path_jobout_subgroup_cv}
		edger ${path_output_cv} ${subgroup} ${working_dir} ${path_cmd} ${path_design_cv} ${path_data} ${type_exec} ${mem_edger} ${walltime_edger} ${path_jobout_subgroup_cv}
		limma ${path_output_cv} ${subgroup} ${working_dir} ${path_cmd} ${path_design_cv} ${path_data} ${type_exec} ${mem_limma} ${walltime_limma} ${path_jobout_subgroup_cv}
	done
fi
