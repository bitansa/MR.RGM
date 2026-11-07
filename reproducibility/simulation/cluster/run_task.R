###############################################################################
## run_task.R  --  one simulation task: one setting of one experiment, all its replicates in parallel
##   Rscript cluster/run_task.R cluster/tasks_sim.csv <task>      (from reproducibility/simulation/)
## It runs the unmodified run_simulations.R / run_instrument_strength.R restricted to this setting, so the result
## files are exactly those of a run on one machine. A setting whose result file exists is skipped (rerun = resume).
###############################################################################
args <- commandArgs(trailingOnly = TRUE)
tk <- read.csv(args[1], stringsAsFactors = FALSE)
tk <- tk[tk$task == as.integer(args[2]), ]; stopifnot(nrow(tk) == 1)
tk$p <- as.numeric(tk$p); tk$n <- as.numeric(tk$n)        # doubles, as in a run on one machine
N_CORES           <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
MAX_CORES_LARGE_N <- N_CORES                  # memory is not limiting on the cluster nodes (1 TB)
REP_TIMEOUT_MIN   <- 600                      # only stops a hung replicate
SKIP_METHODS      <- "GrIVET"                 # fitted in a second pass with a 48-hour limit per fit (cluster/README.md)
message(sprintf("[%s] task %d on %s (%d cores): %s %s p = %d n = %g %s", format(Sys.time(), "%F %T"), tk$task,
                Sys.info()[["nodename"]], N_CORES, tk$experiment, tk$scenario, tk$p, tk$n, tk$label))
if (tk$script == "sim") {
  EXPERIMENT <- tk$experiment; SCENARIOS <- tk$scenario; P_VALUES <- tk$p; N_VALUES <- tk$n; SETTINGS <- tk$setting
  SIM_CHECK  <- FALSE
  source("run_simulations.R")
} else {
  ST_SCENARIOS <- tk$scenario; ST_P_VALUES <- tk$p; ST_N_VALUES <- tk$n; ST_LEVELS <- tk$strength
  STRENGTH_RUN_ONLY <- TRUE
  source("run_instrument_strength.R")
}
message(sprintf("[%s] task %d finished", format(Sys.time(), "%F %T"), tk$task))
