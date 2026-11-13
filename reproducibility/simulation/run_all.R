###############################################################################
## run_all.R  --  the complete simulation study in one run, meant to be started from a Terminal (not RStudio):
##
##   cd path/to/MR.RGM/reproducibility/simulation
##   caffeinate -i Rscript run_all.R > run_all.log 2>&1 &
##
## Why a Terminal: the replicates run in forked R processes, and R's documentation discourages forking inside a GUI
## such as RStudio (forked processes may not return inside a GUI); run_parallel() in sim_functions.R additionally
## terminates and restarts a replicate that exceeds REP_TIMEOUT_MIN (default 240 minutes).
##
## Order: main, roles, instrument strength, noise, confounding, sensitivity, largen, longchain, figures and tables, runtime
## benchmark (last, so its timings are not disturbed by other work). Every step is restartable: settings whose result
## file exists are skipped (the runtime benchmark keeps finished p), so after an interruption the same command continues
## where it stopped. Progress is
## written to run_all.log (the "done in ... min" lines) and errors to results_<experiment>/errors.log.
###############################################################################

t_all <- Sys.time()
stamp <- function(...) message(sprintf("\n========== [%s] %s ==========", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
if (!exists("N_CORES")) N_CORES <- max(1, parallel::detectCores() - 1)
## GrIVET is fitted in a second pass (cluster/add_methods.R; README.md): a replicate that exceeds REP_TIMEOUT_MIN is
## stopped with every method in it, and GrIVET does not finish within 48 hours at n = 10^6
if (!exists("SKIP_METHODS")) SKIP_METHODS <- "GrIVET"

## Each step runs inside tryCatch: an error is reported (run_all.log) and the next step still runs; a rerun of the
## same command repeats only what is missing.
step <- function(label, file, experiment = NULL) {
  stamp(label)
  if (!is.null(experiment)) EXPERIMENT <<- experiment
  ## each script runs in its own environment (child of the global one), so no variable of one step leaks into the next
  ok <- tryCatch({ source(file, local = new.env(parent = globalenv())); TRUE },
                 error = function(e) { message("STEP FAILED (", label, "): ", conditionMessage(e)); FALSE })
  invisible(ok)
}
status <- c()
for (e in c("main", "roles")) status[e] <- step(paste("experiment", e), "run_simulations.R", e)
status["strength"] <- step("instrument-strength study", "run_instrument_strength.R")
for (e in c("noise", "confounding", "sensitivity", "largen", "longchain")) status[e] <- step(paste("experiment", e), "run_simulations.R", e)
status["figures"] <- step("figures and tables", "make_figures_tables.R")
status["runtime"] <- step("runtime benchmark", file.path("Runtime Analysis", "RunTimeAnalysis.R"))
stamp(sprintf("all done in %.1f hours; steps that failed: %s", as.numeric(difftime(Sys.time(), t_all, units = "hours")),
              if (all(status)) "none" else paste(names(status)[!status], collapse = ", ")))
