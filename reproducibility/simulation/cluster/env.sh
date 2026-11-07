#!/bin/bash
###############################################################################
## env.sh  --  sourced by every job: R environment, R library, one BLAS / OpenMP thread per process
## (each replicate is one forked R process on one core; multithreaded BLAS would oversubscribe the cores)
###############################################################################
## config.sh next to this file (found from the path of this file, so that jobs started in any folder can source it)
CLUSTER_CFG="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"
[ -f "$CLUSTER_CFG" ] || { echo "config.sh not found next to env.sh" >&2; exit 1; }
source "$CLUSTER_CFG"
module purge >/dev/null 2>&1
module load "$CONDA_MODULE"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_PREFIX"
export R_LIBS_USER="$R_LIB"
export R_LIBS_SITE=""
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 FLEXIBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
