#!/bin/bash
#SBATCH --job-name=mrrgm_design
#SBATCH --cpus-per-task=100
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --output=design_%j.out
###############################################################################
## check_design_matched.sbatch  --  the design-matched simulation on a SLURM cluster (run from package_check/,
## with the environment of ../simulation/cluster/):
##   sbatch -p <partition> check_design_matched.sbatch   (the #SBATCH resources above are site-specific)
## Reads ../real_data/OneK1K/output/OneK1K_candidates.rds; writes check_design_matched_OneK1K.csv / .rds and prints
## the largest supported gene set at the end of design_<job>.out.
###############################################################################
cd "$SLURM_SUBMIT_DIR"
source ../simulation/cluster/env.sh
Rscript -e 'N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")); source("check_design_matched.R")'
