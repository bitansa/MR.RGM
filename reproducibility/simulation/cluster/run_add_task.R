###############################################################################
## run_add_task.R  --  one task of the add-methods array: completes one result file (cluster/add_methods.R)
##   Rscript cluster/run_add_task.R cluster/tasks_add.csv <task>     (from reproducibility/simulation/)
###############################################################################
args <- commandArgs(trailingOnly = TRUE)
tk <- read.csv(args[1], stringsAsFactors = FALSE)
tk <- tk[tk$task == as.integer(args[2]), ]; stopifnot(nrow(tk) == 1)
N_CORES         <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
## time limit per replicate in minutes (four times as long at n >= 10^6); only stops a hung fit. Default 600; a task
## whose fits need longer can be resubmitted with a larger limit:  REP_TIMEOUT_MIN=1800 sbatch ... cluster/add_array.sbatch
REP_TIMEOUT_MIN <- as.numeric(Sys.getenv("REP_TIMEOUT_MIN", "600"))
## further attempts of a replicate that fails or exceeds the limit (default 1; 0 when the limit is already long)
REP_RETRIES     <- as.integer(Sys.getenv("REP_RETRIES", "1"))
## methods whose failed fits are fitted again (comma-separated, e.g. REFIT_FAILED=GrIVET; default none)
REFIT_FAILED    <- setdiff(trimws(strsplit(Sys.getenv("REFIT_FAILED", ""), ",")[[1]]), "")
## methods that are fitted again for every replicate (comma-separated, e.g. REPLACE_METHODS=MR.RGM+; default none)
REPLACE_METHODS <- setdiff(trimws(strsplit(Sys.getenv("REPLACE_METHODS", ""), ",")[[1]]), "")
ADD_FILES       <- tk$file
message(sprintf("[%s] task %d on %s (%d cores): %s (adding or replacing %s)", format(Sys.time(), "%F %T"), tk$task,
                Sys.info()[["nodename"]], N_CORES, tk$file, tk$missing))
source("cluster/add_methods.R")
message(sprintf("[%s] task %d finished", format(Sys.time(), "%F %T"), tk$task))
