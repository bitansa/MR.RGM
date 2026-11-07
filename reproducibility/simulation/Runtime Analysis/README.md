# Runtime Analysis — `RunTimeAnalysis.R`

`RunTimeAnalysis.R` benchmarks the wall-clock runtime of the methods as the number of responses (p) increases at a
fixed sample size (n = 30,000): scale-free network with 3 instruments per response, p ∈ {2, 5, 10, 20}, the same data
generator, standardization and method settings as the simulation study (`../sim_functions.R`).

Methods: MR.RGM, MR.RGM+, MR.RGM_NoConf (50,000 iterations, 10,000 burn-in, thinning 10), MR2G (package defaults;
skipped if the package is not installed), MVMR, MendelianRandomization (10,000 bootstrap draws, all ordered pairs),
OneSampleMR, mrbayes (3 chains, 10,000 burn-in + 40,000 iterations, all ordered pairs), MrDAG (50,000 iterations,
10,000 burn-in, thinning 10; one outcome at a time with the other responses as exposures, as in the simulation designs
without roles), GrIVET (the tuning grid of the simulation study, unstandardized data). MR2 needs an exposure/outcome split with at least two responses and is timed in the `roles`
experiment instead (column `seconds` of its result files).
A method that cannot be fitted to a data set is reported with `median_sec = NA` and left out of the figure at that p.

On one machine, each method is run `RUNTIME_TIMES` (100, the number of simulation replicates) times per p with
`microbenchmark()` on the same data set; the median is reported. The timings of the paper come from the cluster
instead (`../cluster/README.md`): `cluster/run_runtime_task.R` times the runs with `proc.time()` in blocks (ten runs per
task, two for MrDAG and GrIVET), after an untimed warm-up fit, on one core of the CPU nodes; a block in which a run
stops with an error is dropped, and `cluster/collect_runtime.R` collects the blocks. `RUNTIME_PLOT_ONLY <- TRUE`
redraws the figure from that output. The MR.RGM family works from the sufficient statistics, so its cost does not
grow with n.

Usage, from the `simulation` folder: `source("Runtime Analysis/RunTimeAnalysis.R")`.
Outputs: `Runtime Analysis/runtime_vs_p_n30000.csv` (medians), `runtime_vs_p_n30000_all_runs.csv` (every run),
`runtime_versions.txt` (package versions), `figures/runtime_vs_p_n30000.pdf` (log scale). Redraw only:
`RUNTIME_PLOT_ONLY <- TRUE` before sourcing.
