=====================================================================
EPCY Paper: Framework for Evaluating EPCY
=====================================================================

Introduction
------------

This framework was developed to compare `EPCY <https://github.com/iric-soft/epcy>`_ output with other major DEG tools using Leucegene3 and 10X FACS datasets.  
It consists of several scripts to run all analyses and generate the figures and tables used in the paper.  
**Note:** This is not a general-purpose DEG analysis framework, but a tool to reproduce the results presented in the paper.

Requirements
------------

* Python 3.11.5
  - epcy 0.2.6.4
* R >= 4.2 
  - with dependencies (see `lib_install.r <https://github.com/iric-soft/epcy_paper/blob/master/src/script/other/lib_install.r>`_)

Installation
------------

.. code:: shell

  cd [to_epcy_paper_folder]
  # Set up Python environment and install dependencies
  python3 -m venv venv
  source venv/bin/activate
  pip3 install --upgrade pip setuptools
  pip3 install wheel
  pip3 install epcy==0.2.6.4
  pip3 install -r Requirement.txt

  # Set up R library path and install R dependencies
  mkdir -p R_lib
  echo "R_LIBS=$PWD/R_lib" > .Renviron
  Rscript --vanilla ./src/script/other/lib_install.r

  # Unzip readcount matrices
  cd data/leucegene3/STAR_RSEM
  zcat readcounts.xls.gz | gzip > readcounts.xls
  cd ../../10X_FACS/cellranger
  zcat readcounts.xls.gz | gzip > readcounts.xls

Full Analyses
-------------

* To reproduce all analyses and results from the paper, it is recommended to use a cluster with a scheduler.
* The total data size is about **18 GB**, distributed as:
  - 404 analyses on Leucegene3 (101 analyses for each method)
  - 3200 analyses on 10X FACS (800 analyses for each method)
* Resource requirements for running all analyses are specified in `run_paper.sh`.  
  (These are not optimized; you may be able to use fewer resources.)
* `run_paper.sh` supports both Torque and Slurm schedulers.  
  Set the *type_exec* parameter in `run_paper.sh` to specify your scheduler.

Pipeline Overview
-----------------

The pipeline consists of three main steps:

1. **Prepare Input Files**
   
   .. code:: shell

      cd [your_epcy_paper_folder]
      bash run_init_paper.sh

   This step creates all design and readcount matrix files based on the full Leucegene3 and 10X datasets.

2. **Run All DEG and PG Analyses**
   
   .. code:: shell

      bash run_paper.sh

   - This script has a 'smart' restart mechanism: it checks if output files already exist and only reruns missing or incomplete jobs.
   - If some files are corrupted, remove them and rerun the script.

3. **Generate Figures and Tables**
   
   .. code:: shell

      bash run_res_paper.sh

   - On a MacBook Pro M3 with 8GB RAM, this step takes about **1 hour 30 minutes** and uses up to **2.81 GB** of RAM.

Alternative Approach
--------------------

* Steps 1 and 3 can be performed on a laptop or single computer.
* To skip the resource-intensive step 2, you can download the results of all analyses from Zenodo (link provided in the paper), to run step 3.
* To evaluate reproducibility of step 2, select a subset of analyses to rerun and compare your results with the downloaded ones.
