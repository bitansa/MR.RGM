###############################################################################
## smoke_test.R  --  short end-to-end run of every experiment, the strength study and the runtime benchmark with short
## chains (2 replicates; results_smoke_*/, results_strength_smoke/, figures_smoke/, tables_smoke/), then a check that every method
## produced successful fits. Exit status 1 if a step failed or a method never succeeded (the arrays then do not start).
##   Rscript cluster/smoke_test.R          (from reproducibility/simulation/; called by smoke.sbatch)
###############################################################################
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
## a rerun starts from scratch (the scripts skip settings whose result file exists)
unlink(c(Sys.glob("results_smoke_*"), "results_strength_smoke", "figures_smoke", "tables_smoke"), recursive = TRUE)
unlink(file.path("Runtime Analysis", "runtime_partial_smoke.rds"))
ok <- c()
run <- function(label, expr) {
  message(sprintf("\n===== [%s] smoke: %s =====", format(Sys.time(), "%T"), label))
  ok[label] <<- tryCatch({ expr; TRUE }, error = function(e) { message("FAILED: ", conditionMessage(e)); FALSE })
}
for (e in c("main", "roles", "noise", "confounding", "sensitivity", "largen"))
  run(e, local({ SMOKE_TEST <- TRUE; EXPERIMENT <- e; source("run_simulations.R", local = TRUE) }))
run("strength", local({ STRENGTH_SMOKE <- TRUE; source("run_instrument_strength.R", local = TRUE) }))
run("runtime",  local({ RUNTIME_SMOKE <- TRUE; source(file.path("Runtime Analysis", "RunTimeAnalysis.R"), local = TRUE) }))

## every method has at least one successful fit (a fit that stopped with an error has seconds = NA, failed_fits = 1)
source("sim_functions.R")
row_method <- function(m) ifelse(m %in% c("MR-IVW", "MR-SimpleMedian", "MR-WeightedMedian"), "MendelianRandomization", m)
files <- c(list.files(".", "^results_smoke_", include.dirs = TRUE), "results_strength_smoke")
d <- do.call(rbind, lapply(files, function(dir) {
  fs <- list.files(dir, "^(sim|strength)_.*\\.rds$", full.names = TRUE); if (!length(fs)) return(NULL)
  do.call(rbind, lapply(fs, function(f) { x <- as.data.frame(readRDS(f))
    x <- x[x$family == "diagnostic" & x$metric == "failed_fits", c("method", "value", "seconds")]
    if (nrow(x)) cbind(experiment = sub("^results_(smoke_)?", "", sub("_smoke$", "", dir)), x) else NULL }))
}))
d$method <- row_method(d$method)
d$fits <- 1; d$successful <- as.integer(!is.na(d$seconds)); d$partly_failed <- as.integer(d$successful == 1 & d$value > 0)
d$success <- d$successful == 1
tab <- aggregate(cbind(fits, successful, partly_failed) ~ experiment + method, data = d, FUN = sum)   # partly: some pairs/outcomes failed
print(tab[order(tab$experiment, match(tab$method, ALL_METHODS)), ], row.names = FALSE)
never <- setdiff(ALL_METHODS, unique(d$method[d$success]))
if (length(never)) message("Methods without a successful fit: ", paste(never, collapse = ", "))
for (lf in file.path(files, "errors.log")) if (file.exists(lf)) { message("--- ", lf); cat(tail(readLines(lf), 15), sep = "\n") }
if (!all(ok) || length(never)) { message("SMOKE TEST FAILED: ", paste(names(ok)[!ok], collapse = ", ")); quit(status = 1) }
message("SMOKE TEST PASSED")
