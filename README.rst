
=====================================================================
EPCY paper: framework used to evaluate EPCY
=====================================================================

-------------
Introduction:
-------------

This framework was developed to compare `EPCY <https://github.com/iric-soft/epcy>`_ ouput with other main DEG tools.

--------------
Small analysis:
--------------
This analysis is limited to analyse the smallest design (28_inv16_vs_28), to be run on an single workstation.
Full analysis follow bellow.

Requirements:
-------------

* python3 and dependency

  - pip3 install -r Requirement.txt
* R and dependency (see `lib_install.r <https://github.com/iric-soft/epcy_paper/blob/master/src/script/other/lib_install.r>`_)
* 4 multicore and 2Go
* Take 8 hours on 2GHz intel Core i5 (MacBook Pro)

Usage:
------
.. code:: shell

  cd [your_epcy_paper_folder]
  # Step 1: Create some files used as input in next script
  bash run_init.sh
  # Step 2: Create all cross-validated datasets and run DEG and EPCY analysis
  bash run_ismb_small.sh
  #Note that run_ismb_small.sh have a smart restart mechanism.
  #Rerun it, until all jobs finish their execution correctly
  # Step 3: When all job are completed run
  # Analyse DEG and EPCY outputs and create graphes
  bash run_res_ismb_small.sh

--------------
Full analysis:
--------------

Requirements:
-------------

* python3 and dependency

  - pip3 install -r Requirements.txt
* R and dependency (see `lib_install.r <https://github.com/iric-soft/epcy_paper/blob/master/src/script/other/lib_install.r>`_)
* To run full analysis it's recomanded to use a cluster with a scheduler.
* 8 multicore and

  - 80Go to run Deseq2 in parallel on 400 samples and more.
  - 4Go is enough for Limma, EdgeR and EPCY.

scheduler:
----------
* We have implemented:

  - torque

Create an issue if you need slurm, or add it in DEG_analysis.sh (see exec_cmd function).

Usage:
------
.. code:: shell

  cd [your_epcy_paper_folder]
  # Step 1: Create some files used as input in next script
  bash run_init.sh
  # Step 2: Create all cross-validated datasets and run DEG and EPCY analysis
  bash run_ismb.sh
  #Note that run_ismb.sh have a smart restart mechanism.
  #Rerun it, until all jobs finish their execution correctly
  # Step 3: When all job are completed run
  # Analyse DEG and EPCY outputs and create graphes
  bash run_res_ismb.sh
