#!/bin/bash
###############################################################################
## submit_add.sh  --  adds the methods that are new in sim_functions.R to the existing results (run from
## reproducibility/simulation/ after the package code has been copied to the cluster and setup.sbatch has installed
## any new package):
##   bash cluster/submit_add.sh
## To fit a method again for every replicate (same data, same seeds) after its implementation has changed:
##   REPLACE_METHODS=MR.RGM+ bash cluster/submit_add.sh
## (the stored results of every other method are kept; see cluster/add_methods.R)
## 1. cluster/tasks_add.csv: the result files that lack a method or hold a method to replace (make_add_tasks.R); cluster/tasks_runtime.csv:
##    regenerated (make_tasks.R), so that the runtime benchmark includes the new methods
## 2. smoke test (smoke_add.R) -> add array (one result file per task, its replicates in parallel) -> runtime array
##    (blocks already present are skipped at once) -> final check (cluster/check_results.R, collect_runtime.R)
###############################################################################
set -e
cd "$(dirname "$0")/.."
source cluster/config.sh
source cluster/env.sh
mkdir -p cluster/logs
export REPLACE_METHODS="${REPLACE_METHODS:-}"     # passed on to every job (sbatch exports the environment)
Rscript cluster/make_add_tasks.R
Rscript cluster/make_tasks.R
[ -f cluster/tasks_add.csv ] || { echo "nothing to add"; exit 0; }
NA_=$(( $(wc -l < cluster/tasks_add.csv) - 1 ))
NR=$(( $(wc -l < cluster/tasks_runtime.csv) - 1 ))
J0=$(sbatch --parsable -p "$SETUP_PARTITION" cluster/smoke_add.sbatch)
echo "smoke test:        job $J0 (one replicate per design; the add array starts only if it passes)"
J1=$(sbatch --parsable --dependency=afterok:"$J0" -p "$RUN_PARTITION" ${SIM_EXCLUDE_NODES:+--exclude="$SIM_EXCLUDE_NODES"} --cpus-per-task="$SIM_CPUS" \
       --mem-per-cpu="$SIM_MEM_PER_CPU" --time="$SIM_TIME" --array=1-"$NA_"%"$SIM_MAX_PARALLEL" cluster/add_array.sbatch)
echo "add array:         job $J1 ($NA_ tasks, $SIM_CPUS cores each)"
J2=$(sbatch --parsable --dependency=afterany:"$J1" -p "$RUN_PARTITION" ${RT_EXCLUDE_NODES:+--exclude="$RT_EXCLUDE_NODES"} --mem="$RT_MEM" \
       --time="$RT_TIME" --array=1-"$NR"%"$RT_MAX_PARALLEL" cluster/runtime_array.sbatch)
echo "runtime array:     job $J2 ($NR tasks; existing blocks are skipped)"
J3=$(sbatch --parsable --dependency=afterany:"$J1":"$J2" -p "$SETUP_PARTITION" cluster/check.sbatch)
echo "final check:       job $J3  (result: cluster/logs/check_$J3.out)"
echo "$J1 $J2 $J3" > cluster/logs/last_submission.txt
