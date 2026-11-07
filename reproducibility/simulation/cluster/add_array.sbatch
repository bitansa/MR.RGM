#!/bin/bash
#SBATCH --job-name=mrrgm_add
#SBATCH --output=cluster/logs/add_%A_%a.out
###############################################################################
## add_array.sbatch  --  array job, one task per line of cluster/tasks_add.csv (submitted by submit_add.sh, which sets
## the partition, cores, memory, time and the array range)
###############################################################################
cd "$SLURM_SUBMIT_DIR"
source cluster/env.sh
Rscript cluster/run_add_task.R cluster/tasks_add.csv "$SLURM_ARRAY_TASK_ID"
