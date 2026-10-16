###############################################################################
## run_simulations.R
##
## Runs the full simulation study (all scenarios, network sizes, sample sizes,
## 20 replicates) and saves one result file per setting in results/.
##
## Usage (from this folder):
##   setwd(".../reproducibility/simulation")
##   SMOKE_TEST <- TRUE;  source("run_simulations.R")   # quick check (~minutes), writes results_smoke/
##   SMOKE_TEST <- FALSE; source("run_simulations.R")   # full study
##
## Optional before sourcing:
##   N_CORES   <- 8                                   # default: all cores - 1
##   SCENARIOS <- c("scalefree")                       # run a subset
##   P_VALUES  <- c(10); N_VALUES <- c(500, 1000)       # run a subset
##
## A setting whose result file already exists is skipped, so an interrupted run
## can simply be restarted. Delete a file to recompute that setting.
###############################################################################

source("sim_functions.R")
library(parallel)

if (!exists("SMOKE_TEST")) SMOKE_TEST <- FALSE
if (!exists("N_CORES"))    N_CORES    <- max(1, parallel::detectCores() - 1)

set <- SIM_SETTINGS
if (SMOKE_TEST) {
  ## Tiny run: checks that every method and metric works end-to-end.
  set$n_rep    <- 2
  set$mcmc     <- list(nIter = 2000, nBurnin = 500, Thin = 2)
  set$mr_boot  <- 200
  set$mrbayes  <- list(n.chains = 2, n.burn = 200, n.iter = 1000)
  run_scen <- set$scenarios; run_p <- 5; run_n <- 500
  out_dir  <- "results_smoke"
} else {
  run_scen <- if (exists("SCENARIOS")) SCENARIOS else set$scenarios
  run_p    <- if (exists("P_VALUES"))  P_VALUES  else set$p_values
  run_n    <- if (exists("N_VALUES"))  N_VALUES  else set$n_values
  out_dir  <- "results"
}
dir.create(out_dir, showWarnings = FALSE)

result_file <- function(scenario, p, n) file.path(out_dir, sprintf("sim_%s_p%d_n%d.rds", scenario, p, n))

for (scenario in run_scen) for (p in run_p) for (n in run_n) {
  f <- result_file(scenario, p, n)
  if (file.exists(f)) { message("skip (exists): ", f); next }
  message(sprintf("[%s] %s  p = %d  n = %d  (%d replicates on %d cores)",
                  format(Sys.time(), "%H:%M:%S"), scenario, p, n, set$n_rep, N_CORES))
  t0  <- Sys.time()
  res <- mclapply(seq_len(set$n_rep), run_replicate, scenario = scenario, p = p, n = n, set = set,
                  mc.cores = N_CORES, mc.preschedule = FALSE)
  bad <- vapply(res, function(r) inherits(r, "try-error") || is.null(r), logical(1))
  if (any(bad)) warning(sprintf("%d replicate(s) failed completely in %s p=%d n=%d", sum(bad), scenario, p, n))
  df  <- do.call(rbind, res[!bad])
  attr(df, "settings")     <- set
  attr(df, "session_info") <- utils::sessionInfo()
  saveRDS(df, f)
  message(sprintf("   done in %.1f min -> %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), f))
}

## Quick completeness check: replicates per method and metric family
all_res <- do.call(rbind, lapply(list.files(out_dir, pattern = "^sim_.*\\.rds$", full.names = TRUE), readRDS))
chk <- aggregate(rep ~ scenario + p + n + method + family,
                 data = unique(all_res[!is.na(all_res$value), c("scenario", "p", "n", "method", "family", "rep")]),
                 FUN = length)
names(chk)[names(chk) == "rep"] <- "replicates"
incomplete <- chk[chk$replicates < set$n_rep, ]
if (nrow(incomplete)) { message("Settings with fewer than ", set$n_rep, " valid replicates:"); print(incomplete) } else
  message("All settings complete: ", set$n_rep, " replicates for every method and metric family.")
fails <- all_res[all_res$family == "diagnostic" & all_res$value > 0, ]
if (nrow(fails)) { message("Pairwise fits that errored (kept as estimate 0, p = 1):"); print(fails[, c("scenario","p","n","rep","method","value")]) }
