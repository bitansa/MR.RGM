#!/bin/bash
#SBATCH --job-name=mrrgm_check
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=cluster/logs/check_%j.out
cd "$SLURM_SUBMIT_DIR"
source cluster/env.sh
Rscript cluster/check_results.R
Rscript cluster/collect_runtime.R
