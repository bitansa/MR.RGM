###############################################################################
## make_tasks.R  --  task lists of the SLURM arrays (run from reproducibility/simulation/; called by setup.sbatch)
##
## cluster/tasks_sim.csv      one task per setting of every experiment (run_simulations.R) and of the instrument-strength
##                            study (run_instrument_strength.R); each task runs the setting's 100 replicates in parallel
## cluster/tasks_runtime.csv  runtime benchmark: one task per (p, method, block of runs); the 100 runs of a method at a
##                            given p are split into blocks so that no task runs for days
## The grids are those of the scripts (experiment_grid(), STRENGTH_GRID, STRENGTH and RUNTIME_EXPRS in sim_functions.R).
###############################################################################
source("sim_functions.R")
rows <- list()
for (e in c("main", "roles", "noise", "confounding", "sensitivity", "largen", "longchain")) {
  g <- experiment_grid(e)
  for (sc in g$scen) for (p in g$p) for (n in g$n) for (st in g$settings)
    rows[[length(rows) + 1]] <- data.frame(script = "sim", experiment = e, scenario = sc, p = p, n = n,
                                           setting = st$index + 1, label = st$label, strength = "")
}
for (sc in STRENGTH_GRID$scen) for (p in STRENGTH_GRID$p) for (lv in names(STRENGTH)) for (n in STRENGTH_GRID$n)
  rows[[length(rows) + 1]] <- data.frame(script = "strength", experiment = "strength", scenario = sc, p = p, n = n,
                                         setting = NA, label = lv, strength = lv)
ts <- do.call(rbind, rows)
## largest tasks first (they start first and set the total time): p, then the number of instruments, then n
ts <- ts[order(-ts$p, -(ts$strength == "Weak10"), -ts$n), ]
ts$task <- seq_len(nrow(ts))
write.csv(ts[, c("task", setdiff(names(ts), "task"))], "cluster/tasks_sim.csv", row.names = FALSE)

RUNTIME_P    <- c(2, 5, 10, 20)
RUNTIME_RUNS <- 100
blocks <- function(m) if (m %in% c("MrDAG", "GrIVET")) 50 else 10      # the slowest methods: 2 runs per task
tr <- do.call(rbind, lapply(RUNTIME_P, function(p) do.call(rbind, lapply(names(RUNTIME_EXPRS), function(m) {
  k <- blocks(m); size <- RUNTIME_RUNS / k
  data.frame(p = p, method = m, block = seq_len(k), first_run = (seq_len(k) - 1) * size + 1, n_runs = size)
}))))
tr <- tr[order(-tr$p, -(tr$method %in% c("MrDAG", "GrIVET"))), ]
tr$task <- seq_len(nrow(tr))
write.csv(tr[, c("task", setdiff(names(tr), "task"))], "cluster/tasks_runtime.csv", row.names = FALSE)
message(sprintf("%d simulation tasks, %d runtime tasks (MaxArraySize must be at least %d)", nrow(ts), nrow(tr), max(nrow(ts), nrow(tr))))
