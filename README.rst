
=====================================================================
EPCY paper: framework which regroup all analyses done to publish EPCY
=====================================================================

-------------
Introduction:
-------------

This framework was developed to analyse `_EPCY <https://github.com/iric-soft/epcy>` ouput with other main DEG tools.

-------------
Requirements:
-------------

* python3
* R and dependency (see `_lib_install.r <https://github.com/iric-soft/epcy_paper/blob/master/src/script/other/lib_install.r>)`
* A cluster with a scheduler is highly recommended to run full analysis. Using a workstation/laptop run the small version.

  - scheduler implemented: torque, slurm(soon)

-------------
Usage:
-------------

.. code:: shell

  $ cd [your_epcy_paper folder_folder]
  # will create some files used used as input in next script
  $ bash run_init.sh
  # Create all cross-validated datasets and run DEG and EPCY analysis
  $ bash run_ismb.sh
  # Analyse DEG and EPCY outputs and create graphe
  $ bash run_res_ismb.sh
