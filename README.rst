
=====================================================================
EPCY paper: framework used to evaluate EPCY 
=====================================================================

-------------
Introduction:
-------------

This framework was developed to analyse `EPCY <https://github.com/iric-soft/epcy>`_ ouput with other main DEG tools.

-------------
Requirements:
-------------

* python3 and dependency (coming soon)
* R and dependency (see `lib_install.r <https://github.com/iric-soft/epcy_paper/blob/master/src/script/other/lib_install.r>`_)
* To run full analysis it's recomanded to use a cluster with a scheduler. 

  - scheduler implemented: torque, slurm(soon)
  
* Using a workstation/laptop check the small version. (coming soon)

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
