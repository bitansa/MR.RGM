#!/bin/bash
###############################################################################
## status.sh  --  progress of the cluster run (run from reproducibility/simulation/):  bash cluster/status.sh
###############################################################################
cd "$(dirname "$0")/.."
echo "== jobs =="; squeue -u "$USER" -o "%.12i %.14j %.9T %.10M %.6D %.5C %R" | head -30
NS=$(( $(wc -l < cluster/tasks_sim.csv) - 1 )); NR=$(( $(wc -l < cluster/tasks_runtime.csv) - 1 ))
DS=$(ls results_{main,roles,noise,confounding,sensitivity,largen,longchain}/sim_*.rds results_strength/strength_*.rds 2>/dev/null | wc -l)
DR=$(ls "Runtime Analysis/cluster_parts/"rt_*.rds 2>/dev/null | wc -l)
echo "== progress =="; echo "simulation settings done: $DS / $NS"; echo "runtime blocks done:      $DR / $NR"
echo "== last errors =="; tail -n 3 results_*/errors.log 2>/dev/null | tail -n 12
