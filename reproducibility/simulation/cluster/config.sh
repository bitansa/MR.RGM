#!/bin/bash
###############################################################################
## config.sh  --  settings of the SLURM pipeline (the only file to edit for another cluster)
##
## Layout on the cluster: the package folder MR.RGM/ (with reproducibility/) is copied to $CLUSTER_ROOT/MR.RGM;
## the R environment and the R library live next to it; every job runs in reproducibility/simulation/.
###############################################################################
CLUSTER_ROOT="$HOME/MRRGM"                       # top folder of the study on the cluster
PKG_DIR="$CLUSTER_ROOT/MR.RGM"                   # package source (copied from the local machine)
SIM_DIR="$PKG_DIR/reproducibility/simulation"    # working folder of every job
ENV_PREFIX="$CLUSTER_ROOT/conda-env"             # conda environment: R, JAGS, compilers, system libraries
R_LIB="$CLUSTER_ROOT/Rlib"                       # R packages installed from CRAN / GitHub / r-universe
CONDA_MODULE="Anaconda3/2024.02-1"               # module that provides conda

## Partitions and resources
SETUP_PARTITION="short"                          # setup, smoke test and final check (2-hour limit)
RUN_PARTITION="long"                             # simulation and runtime tasks (7-day limit)
SIM_EXCLUDE_NODES=""                             # simulation tasks: every node (CPU and GPU nodes, 256 cores and 1 TB each)
RT_EXCLUDE_NODES="g01,g02,g03"                   # runtime tasks: the 9 CPU nodes only (identical AMD EPYC 7763 nodes)
SIM_CPUS=100                                     # cores per simulation task: the 100 replicates of a setting run at once
SIM_MEM_PER_CPU="4G"                             # up to about 2.5 GB per replicate at n = 10^6
SIM_TIME="2-00:00:00"
SIM_MAX_PARALLEL=1000                            # no limit: as many tasks at once as the free cores allow (2 per node)
RT_MEM="8G"                                      # runtime tasks: one core each
RT_TIME="3-00:00:00"
RT_MAX_PARALLEL=1000                             # every runtime task at once (one core each)
