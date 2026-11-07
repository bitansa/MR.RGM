#!/bin/bash
#SBATCH --job-name=mrrgm_setup
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=cluster/logs/setup_%j.out
###############################################################################
## setup.sbatch  --  builds the software environment once (run from reproducibility/simulation/):
##   sbatch -p short cluster/setup.sbatch
## 1. conda environment (conda-forge): R 4.4, JAGS + rjags, compilers and the system libraries that CRAN packages need
## 2. the remaining R packages (install_packages.R): CRAN, GitHub (MR2G, grivet), r-universe (MR2, MrDAG)
## 3. MR.RGM from the copied source; 4. checks (verify_setup.R); 5. task lists (make_tasks.R)
## Rerunning it skips what is already installed.
###############################################################################
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
source cluster/config.sh
module purge >/dev/null 2>&1
module load "$CONDA_MODULE"
source "$(conda info --base)/etc/profile.d/conda.sh"
export CONDA_PKGS_DIRS="$CLUSTER_ROOT/conda-pkgs"
if [ ! -d "$ENV_PREFIX" ]; then
  echo "[$(date '+%F %T')] creating the conda environment $ENV_PREFIX"
  conda create -y -p "$ENV_PREFIX" --override-channels -c conda-forge \
    r-base=4.4 jags r-rjags r-igraph r-rcpp r-rcpparmadillo r-proc r-ggplot2 r-scales r-microbenchmark \
    r-remotes r-future r-testthat r-glmnet r-quantreg r-robustbase r-plotly r-lme4 r-car r-gmm \
    compilers make cmake pkg-config gmp mpfr nlopt glpk gsl libxml2 zlib libcurl openssl
fi
source cluster/env.sh
mkdir -p "$R_LIB"
echo "[$(date '+%F %T')] R packages"
Rscript cluster/install_packages.R
echo "[$(date '+%F %T')] MR.RGM from $PKG_DIR"
R CMD INSTALL --preclean --library="$R_LIB" "$PKG_DIR"
echo "[$(date '+%F %T')] checks"
Rscript cluster/verify_setup.R
echo "[$(date '+%F %T')] task lists"
Rscript cluster/make_tasks.R
echo "[$(date '+%F %T')] setup finished"
