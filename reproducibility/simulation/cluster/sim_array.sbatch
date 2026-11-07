#!/bin/bash
#SBATCH --job-name=mrrgm_sim
#SBATCH --output=cluster/logs/sim_%A_%a.out
###############################################################################
## sim_array.sbatch  --  array job, one task per line of cluster/tasks_sim.csv (submitted by submit_all.sh, which sets
## the partition, cores, memory, time and the array range)
###############################################################################
cd "$SLURM_SUBMIT_DIR"
source cluster/env.sh
Rscript cluster/run_task.R cluster/tasks_sim.csv "$SLURM_ARRAY_TASK_ID"
