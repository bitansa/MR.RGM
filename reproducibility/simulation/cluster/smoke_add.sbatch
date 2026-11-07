#!/bin/bash
#SBATCH --job-name=mrrgm_smoke_add
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=cluster/logs/smoke_add_%j.out
###############################################################################
## smoke_add.sbatch  --  quick check of the added methods (cluster/smoke_add.R); submitted by submit_add.sh
###############################################################################
cd "$SLURM_SUBMIT_DIR"
source cluster/env.sh
Rscript cluster/smoke_add.R
