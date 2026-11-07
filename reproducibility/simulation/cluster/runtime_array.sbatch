#!/bin/bash
#SBATCH --job-name=mrrgm_runtime
#SBATCH --cpus-per-task=1
#SBATCH --output=cluster/logs/runtime_%A_%a.out
###############################################################################
## runtime_array.sbatch  --  array job, one task per line of cluster/tasks_runtime.csv (one core each; submitted by
## submit_all.sh after the simulation tasks have finished, so that the runs do not share nodes with simulation tasks)
###############################################################################
cd "$SLURM_SUBMIT_DIR"
source cluster/env.sh
Rscript cluster/run_runtime_task.R cluster/tasks_runtime.csv "$SLURM_ARRAY_TASK_ID"
