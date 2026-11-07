#!/bin/bash
###############################################################################
## submit_all.sh  --  submits the whole study as dependent SLURM jobs (run from reproducibility/simulation/):
##   bash cluster/submit_all.sh             smoke test -> simulation array -> runtime array -> final check
##   bash cluster/submit_all.sh --no-smoke  without the smoke test (e.g. after it has passed once)
## The simulation array starts only if the smoke test passed; the runtime array starts when the simulation array has
## finished (timings on lightly loaded nodes); the check runs at the end in any case. Resubmitting is safe: finished
## settings and runtime blocks are skipped.
###############################################################################
set -e
cd "$(dirname "$0")/.."
source cluster/config.sh
[ -f cluster/tasks_sim.csv ] && [ -f cluster/tasks_runtime.csv ] || { echo "run cluster/setup.sbatch first"; exit 1; }
mkdir -p cluster/logs
NS=$(( $(wc -l < cluster/tasks_sim.csv) - 1 ))
NR=$(( $(wc -l < cluster/tasks_runtime.csv) - 1 ))
DEP=""
if [ "$1" != "--no-smoke" ]; then
  J0=$(sbatch --parsable -p "$SETUP_PARTITION" cluster/smoke.sbatch)
  echo "smoke test:        job $J0"
  DEP="--dependency=afterok:$J0"
fi
J1=$(sbatch --parsable $DEP -p "$RUN_PARTITION" ${SIM_EXCLUDE_NODES:+--exclude="$SIM_EXCLUDE_NODES"} --cpus-per-task="$SIM_CPUS" \
       --mem-per-cpu="$SIM_MEM_PER_CPU" --time="$SIM_TIME" --array=1-"$NS"%"$SIM_MAX_PARALLEL" cluster/sim_array.sbatch)
echo "simulation array:  job $J1 ($NS tasks, at most $SIM_MAX_PARALLEL at a time, $SIM_CPUS cores each)"
J2=$(sbatch --parsable --dependency=afterany:"$J1" -p "$RUN_PARTITION" ${RT_EXCLUDE_NODES:+--exclude="$RT_EXCLUDE_NODES"} --mem="$RT_MEM" \
       --time="$RT_TIME" --array=1-"$NR"%"$RT_MAX_PARALLEL" cluster/runtime_array.sbatch)
echo "runtime array:     job $J2 ($NR tasks, one core each)"
J3=$(sbatch --parsable --dependency=afterany:"$J1":"$J2" -p "$SETUP_PARTITION" cluster/check.sbatch)
echo "final check:       job $J3  (result: cluster/logs/check_$J3.out)"
echo "$J1 $J2 $J3" > cluster/logs/last_submission.txt
