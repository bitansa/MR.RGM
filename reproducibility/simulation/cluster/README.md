# Running the simulation study on a SLURM cluster

The same scripts as on one machine (`run_simulations.R`, `run_instrument_strength.R`, `sim_functions.R`), run one
setting per task: every setting of every experiment and of the instrument-strength study is one task that runs its
100 replicates at once on 100 cores; the runtime benchmark is split into single-core tasks. The result files are
exactly those of a run on one machine, in the same folders, so the figures and tables are made locally afterwards
with `make_figures_tables.R`, `run_instrument_strength.R` (`STRENGTH_PLOT_ONLY <- TRUE`) and
`Runtime Analysis/RunTimeAnalysis.R` (`RUNTIME_PLOT_ONLY <- TRUE`).

| File | Purpose |
|---|---|
| `config.sh` | Paths, partitions, cores, memory, time limits (the only file to edit for another cluster) |
| `env.sh` | Sourced by every job: conda R environment, R library, one BLAS thread per process |
| `setup.sbatch` | Builds the environment once: conda (R 4.4, JAGS, compilers, libraries), `install_packages.R`, MR.RGM from source, `verify_setup.R`, `make_tasks.R` |
| `make_tasks.R` | `tasks_sim.csv` (one task per setting, 118) and `tasks_runtime.csv` (one task per p, method and block of runs, 720) |
| `smoke_test.R`, `smoke.sbatch` | Every experiment with short chains and 2 replicates; fails if a step fails or a method never succeeds |
| `run_task.R`, `sim_array.sbatch` | One simulation task (one setting, all methods, 100 replicates in parallel) |
| `run_runtime_task.R`, `runtime_array.sbatch` | One runtime task: an untimed warm-up fit, then 2 (MrDAG, GrIVET) or 10 timed runs of one method at one p on one core |
| `check_results.R`, `collect_runtime.R`, `check.sbatch` | Completeness per task and method (`check_results.csv`); runtime CSVs of `RunTimeAnalysis.R` |
| `submit_all.sh` | Smoke test → simulation array → runtime array → check, as dependent jobs |
| `status.sh` | Progress (jobs, settings and runtime blocks done, last errors) |
| `add_methods.R`, `make_add_tasks.R`, `run_add_task.R`, `add_array.sbatch`, `smoke_add.R`, `smoke_add.sbatch`, `submit_add.sh` | Adds methods to the existing result files (the second pass that fits GrIVET, or any method added to `METHODS`; only the missing methods and replicates, same data and seeds), then completes the runtime benchmark and runs the check |

## Workflow

Local machine (copy the package and its `reproducibility/` folder; build files and outputs are not copied):

```bash
rsync -av --exclude '.git' --exclude '.Rproj.user' --exclude '*.o' --exclude '*.so' --exclude '.DS_Store' \
  --exclude '.Rhistory' --exclude '.RData' --exclude 'results_*' --exclude 'figures*' --exclude 'tables*' \
  "path/to/MR.RGM/" user@cluster:MRRGM/MR.RGM/
```

Cluster (login node; everything else runs as jobs):

```bash
cd $CLUSTER_ROOT/MR.RGM/reproducibility/simulation   # CLUSTER_ROOT, partitions and resources: cluster/config.sh
mkdir -p cluster/logs
sbatch -p <SETUP_PARTITION> cluster/setup.sbatch   # about 30-60 min; check: tail cluster/logs/setup_*.out ("setup finished")
bash cluster/submit_all.sh                    # smoke test, then every method except GrIVET
REP_TIMEOUT_MIN=2880 REP_RETRIES=0 bash cluster/submit_add.sh   # second pass: GrIVET, 48 hours per fit
bash cluster/status.sh                        # progress, at any time
cat cluster/logs/check_*.out                  # at the end: "EVERYTHING COMPLETE"
```

If the smoke test fails, the other jobs stay pending ("DependencyNeverSatisfied"): `scancel -u $USER`, read
`cluster/logs/smoke_*.out`, fix, and submit again. A task that fails can be resubmitted alone
(`sbatch -p <RUN_PARTITION> --cpus-per-task=100 --mem-per-cpu=4G --time=2-00:00:00 --array=<task> cluster/sim_array.sbatch`);
finished settings are skipped, so `bash cluster/submit_all.sh --no-smoke` also completes whatever is missing.

**GrIVET (second pass).** The simulation tasks skip GrIVET (`SKIP_METHODS` in `run_task.R`), since a replicate that
exceeds its time limit is stopped with every method in it, and GrIVET does not finish within 48 hours at n = 10^6.
`submit_add.sh` then adds GrIVET to every result file, with `REP_TIMEOUT_MIN=2880` (48 hours per fit) and the 48-hour
job limit (`SIM_TIME` in `config.sh`); a fit without a result at the limit is recorded as not finished. The same
command adds any other method added to `METHODS` later (after installing its package with `setup.sbatch`; installed
packages are skipped). Only the missing methods and replicates are run, on the same data and seeds; the original
result files are kept in `_Archive before added methods/`.

Local machine (results back into the simulation folder; move earlier results of a local run to an archive first):

```bash
cd "path/to/simulation"
rsync -av user@cluster:MRRGM/MR.RGM/reproducibility/simulation/results_{main,roles,strength,noise,confounding,sensitivity,largen} ./
rsync -avs "user@cluster:MRRGM/MR.RGM/reproducibility/simulation/Runtime Analysis/runtime_*" "Runtime Analysis/"
rsync -av user@cluster:MRRGM/MR.RGM/reproducibility/simulation/cluster/check_results.csv cluster/
```

## Notes

* Seeds do not depend on where or with which other methods a replicate runs (`sim_functions.R`), so a setting gives
  the same result as in a run on one machine with the same software; numerical libraries (BLAS) can differ in the
  last digits between machines.
* One BLAS thread per process (`env.sh`): each replicate is one forked R process on one core.
* Resources (`config.sh`): simulation tasks on every node with no limit on the number running at once (100 cores and
  400 GB per task, two per node); the runtime tasks start when the simulation tasks have finished, each on one core of
  the nine identical CPU nodes (AMD EPYC 7763); 100 runs per method and p, median reported, as in `RunTimeAnalysis.R`.
