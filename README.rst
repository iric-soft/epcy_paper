
=====================================================================
EPCY paper: framework used to evaluate EPCY
=====================================================================

-------------
Introduction:
-------------

This framework was developed to analyse `EPCY <https://github.com/iric-soft/epcy>`_ ouput with other main DEG tools.

--------------
Small analysis:
--------------

Requirements:
-------------

* python3 and dependency

  - pip install statsmodels
* R and dependency (see `lib_install.r <https://github.com/iric-soft/epcy_paper/blob/master/src/script/other/lib_install.r>`_)
* 4 multicore and 8Go

Usage:
------
(coming soon)

--------------
Full analysis:
--------------

Requirements:
-------------

* python3 and dependency

  - pip install statsmodels
* R and dependency (see `lib_install.r <https://github.com/iric-soft/epcy_paper/blob/master/src/script/other/lib_install.r>`_)
* To run full analysis it's recomanded to use a cluster with a scheduler. Without check the small version.
* 8 multicore and

  - 80Go to run Deseq2 in parallel on 400 samples and more.
  - 4Go is enough for Limma, EdgeR and EPCY.

scheduler:
----------
* We have implemented:

  - torque
  - slurm(coming soon)

Usage:
------
.. code:: shell

  cd [your_epcy_paper folder_folder]
  # will create some files used used as input in next script
  bash run_init.sh
  # Create all cross-validated datasets and run DEG and EPCY analysis
  bash run_ismb.sh
  #run_ismb.sh have a smart restart mechanism.
  #Juste rerun it, until all jobs finish their execution correctly
  #############
  # After all job are completed run
  bash run_res_ismb.sh
  # To analyse DEG and EPCY outputs and create graphes
