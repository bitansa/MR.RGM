###############################################################################
## run_simulations.R  --  runs one experiment of the simulation study and saves one result file per setting
##
## Experiments (EXPERIMENT <- "..." before sourcing; default "main"):
##   main         Cases 1-3, p in {5, 10}, n in {500, 1000, 10000, 30000}
##   noise        error variance in {1, 4, 9} (small-world, p = 10, n in {1000, 10000})
##   confounding  confounder loadings x {0, 0.5, 1, 1.5} (small-world, p = 10, n in {1000, 10000})
##   largen       n in {1e5, 1e6} (small-world, p = 10)
##   longchain    MR.RGM+ at n = 1e6 with 250,000 iterations (same data and seeds as largen)
##   roles        4 exposures -> 4 outcomes with known roles (n in {1000, 10000})
##   sensitivity  MR.RGM family with prior inclusion Beta(3,1) / Beta(1,1), nu1 in {1e-2, 1e-3, 1e-4}
##                (small-world, p = 10, n in {1000, 10000})
## Methods: every method (METHODS in sim_functions.R) in every experiment except sensitivity (MR.RGM family only) and
## longchain (MR.RGM+ only).
##
## Usage (from this folder):
##   SMOKE_TEST <- TRUE;  source("run_simulations.R")                  # short end-to-end check of every method (~minutes)
##   EXPERIMENT <- "main"; source("run_simulations.R")                 # full experiment (results_main/)
## Optional before sourcing: N_CORES (default all cores - 1); SCENARIOS, P_VALUES, N_VALUES, SETTINGS (subsets)
##
## A setting whose result file already exists is skipped, so an interrupted run can be restarted; delete a file to recompute.
## Every result file stores the settings and the session info (package versions) as attributes. Next to each
## sim_*.rds (metrics) a fits_*.rds holds the truth and every method's estimates of every replicate, so metrics can be
## recomputed or added later without rerunning. At n = 1e6 at most MAX_CORES_LARGE_N (default 4) workers.
###############################################################################

source("sim_functions.R")
library(parallel)

if (!exists("SMOKE_TEST")) SMOKE_TEST <- FALSE
if (!exists("EXPERIMENT")) EXPERIMENT <- "main"
if (!exists("N_CORES"))    N_CORES    <- max(1, parallel::detectCores() - 1)
if (!exists("REP_TIMEOUT_MIN")) REP_TIMEOUT_MIN <- 240  # a replicate running longer is killed and restarted once (twice the limit at n = 1e6)
if (!exists("SIM_CHECK"))  SIM_CHECK  <- TRUE           # FALSE: skip the completeness check at the end (one setting per cluster task)
stopifnot(EXPERIMENT %in% setdiff(names(METHODS), "strength"))   # strength: run_instrument_strength.R

set <- SIM_SETTINGS
methods <- METHODS[[EXPERIMENT]]
if (exists("SKIP_METHODS")) methods <- setdiff(methods, SKIP_METHODS)   # methods fitted in a separate pass (cluster/add_methods.R)

## Experiment-specific grid (experiment_grid() in sim_functions.R): scenarios, p, n and the list of settings
grid <- experiment_grid(EXPERIMENT, set)

if (SMOKE_TEST) {
  ## Tiny run: checks that every method and metric works end-to-end (results_smoke/)
  set$n_rep   <- 2
  set$mcmc    <- list(nIter = 2000, nBurnin = 500, Thin = 2)
  set$mr_boot <- 200
  set$mrbayes <- list(n.chains = 3, n.burn = 200, n.iter = 800)
  set$mr2g    <- modifyList(set$mr2g, list(num_pert = 5, num_pert_M = 5))
  set$mr2     <- modifyList(set$mr2, list(niter = 600, burnin = 200, thin = 2))
  set$mrdag   <- modifyList(set$mrdag, list(niter = 600, burnin = 200, thin = 2))
  grid$n <- grid$n[1]; grid$p <- grid$p[1]; grid$settings <- grid$settings[unique(c(1, min(2, length(grid$settings)), length(grid$settings)))]   # first, second, last
  if (EXPERIMENT == "main") grid$n <- 500
  out_dir <- paste0("results_smoke_", EXPERIMENT)
} else {
  if (exists("SCENARIOS")) grid$scen <- SCENARIOS
  if (exists("P_VALUES"))  grid$p    <- P_VALUES
  if (exists("N_VALUES"))  grid$n    <- N_VALUES
  if (exists("SETTINGS"))  grid$settings <- grid$settings[SETTINGS]
  out_dir <- paste0("results_", EXPERIMENT)
}
dir.create(out_dir, showWarnings = FALSE)
options(sim_error_log = file.path(out_dir, "errors.log"))            # error messages of failed fits (see safely())

result_file <- function(scenario, p, n, label) file.path(out_dir, sprintf("sim_%s_%s_p%d_n%g_%s.rds", EXPERIMENT, scenario, p, n, label))

for (scenario in grid$scen) for (p in grid$p) for (n in grid$n) for (si in seq_along(grid$settings)) {
  ## Memory: at n = 1e6 every worker holds X, Y and their copies (about 1-1.5 GB), so at most MAX_CORES_LARGE_N workers
  n_cores <- if (n >= 1e6) min(N_CORES, if (exists("MAX_CORES_LARGE_N")) MAX_CORES_LARGE_N else 4) else N_CORES
  st <- grid$settings[[si]]                                # st$index: position in the full grid (column `setting`)
  f  <- result_file(scenario, p, n, st$label)
  if (check_existing_result(f)) { message("skip (exists): ", f); next }
  set_i <- modifyList(set, st[setdiff(names(st), c("label", "index"))])
  methods_i <- methods_for_setting(EXPERIMENT, st$label, methods)
  message(sprintf("[%s] %s: %s  p = %d  n = %g  setting = %s  (%d replicates on %d cores; methods: %s)",
                  format(Sys.time(), "%H:%M:%S"), EXPERIMENT, scenario, p, n, st$label, set_i$n_rep, n_cores, paste(methods_i, collapse = ", ")))
  t0  <- Sys.time()
  res <- run_parallel(seq_len(set_i$n_rep), run_replicate, scenario = scenario, p = p, n = n, set = set_i,
                      experiment = EXPERIMENT, setting = st$index, setting_label = st$label, methods = methods_i,
                      cores = n_cores, timeout_min = if (n >= 1e6) 2 * REP_TIMEOUT_MIN else REP_TIMEOUT_MIN,
                      label = sprintf("%s %s p=%d n=%g %s", EXPERIMENT, scenario, p, n, st$label))
  bad <- vapply(res, function(r) inherits(r, "try-error") || is.null(r), logical(1))
  if (any(bad)) warning(sprintf("%d replicate(s) failed completely in %s p=%d n=%g %s", sum(bad), scenario, p, n, st$label))
  if (all(bad)) next                                       # nothing to save; rerun this setting
  fits <- lapply(res[!bad], attr, "fits")                  # truth and estimates of every replicate and method
  df  <- do.call(rbind, lapply(res[!bad], function(r) { attr(r, "fits") <- NULL; r }))
  attr(df, "settings")     <- set_i
  attr(df, "methods")      <- methods_i
  attr(df, "sampler")      <- SAMPLER_TAG
  attr(df, "MR.RGM_build") <- unlist(utils::packageDescription("MR.RGM")[c("Version", "Packaged", "RemoteSha", "Built")])
  attr(df, "session_info") <- utils::sessionInfo()
  attr(df, "package_versions") <- package_versions()
  save_rds_atomic(list(settings = set_i, methods = methods_i, fits = fits), sub("/sim_", "/fits_", f, fixed = TRUE))   # estimates
  save_rds_atomic(df, f)                                     # metrics (written last: its presence marks the setting as done)
  message(sprintf("   done in %.1f min -> %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), f))
}

## Completeness check: valid replicates per method and metric family
if (SIM_CHECK) {
all_res <- do.call(rbind, lapply(list.files(out_dir, pattern = "^sim_.*\\.rds$", full.names = TRUE), readRDS))
if (is.null(all_res)) stop("No result file in ", out_dir, ": every setting failed or was skipped.")
chk <- aggregate(rep ~ scenario + p + n + setting_label + method + family,
                 data = unique(all_res[!is.na(all_res$value), c("scenario", "p", "n", "setting_label", "method", "family", "rep")]),
                 FUN = length)
names(chk)[names(chk) == "rep"] <- "replicates"
incomplete <- chk[chk$replicates < set$n_rep, ]
if (nrow(incomplete)) { message("Settings with fewer than ", set$n_rep, " valid replicates:"); print(incomplete) } else
  message("All settings complete: ", set$n_rep, " replicates for every method and metric family.")
fails <- all_res[all_res$family == "diagnostic" & all_res$metric == "failed_fits" & all_res$value > 0, ]
if (nrow(fails)) { message("Fits that errored or failed for some pairs (counts):"); print(fails[, c("scenario", "p", "n", "setting_label", "rep", "method", "value")]) }
if (file.exists(getOption("sim_error_log"))) { message("Error messages (", getOption("sim_error_log"), "):"); cat(tail(readLines(getOption("sim_error_log")), 20), sep = "\n") }
}
