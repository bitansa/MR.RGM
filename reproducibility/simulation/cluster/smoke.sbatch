#!/bin/bash
#SBATCH --job-name=mrrgm_smoke
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --output=cluster/logs/smoke_%j.out
cd "$SLURM_SUBMIT_DIR"
source cluster/env.sh
Rscript cluster/smoke_test.R
