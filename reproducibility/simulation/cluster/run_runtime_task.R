###############################################################################
## run_runtime_task.R  --  one runtime task: n_runs timed fits of one method at one p, one after the other on one core
##   Rscript cluster/run_runtime_task.R cluster/tasks_runtime.csv <task>
## Same data set (runtime_data(): scale-free network, n = 30,000), same call (RUNTIME_EXPRS) and settings as
## Runtime Analysis/RunTimeAnalysis.R; wall-clock seconds of each run. A method that cannot be fitted to the data set
## is recorded as failed (median NA). Before the timed runs, one untimed
## warm-up fit (as the trial fit of RunTimeAnalysis.R): it also detects a method that cannot be fitted.
## Output: Runtime Analysis/cluster_parts/rt_p<p>_m<method number in RUNTIME_EXPRS>_<block>.rds
###############################################################################
args <- commandArgs(trailingOnly = TRUE)
tk <- read.csv(args[1], stringsAsFactors = FALSE)
tk <- tk[tk$task == as.integer(args[2]), ]; stopifnot(nrow(tk) == 1)
source("sim_functions.R")
out_dir <- file.path("Runtime Analysis", "cluster_parts"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
f <- file.path(out_dir, sprintf("rt_p%02d_m%02d_%02d.rds", tk$p, match(tk$method, names(RUNTIME_EXPRS)), tk$block))   # m = method number
if (file.exists(f)) { message("exists: ", f); quit(save = "no") }
set <- SIM_SETTINGS
rd  <- runtime_data(tk$p, 30000, set); X <- rd$X; Y <- rd$Y; D <- rd$D; Dfull <- rd$Dfull; ms <- rd$ms
X_raw <- rd$X_raw; Y_raw <- rd$Y_raw; sY <- rd$sY                   # GrIVET (unstandardized data)
expr <- RUNTIME_EXPRS[[tk$method]]
runs <- tk$first_run + seq_len(tk$n_runs) - 1
sec  <- rep(NA_real_, length(runs))
failed <- tk$method %in% names(METHOD_PACKAGE) && !has_pkg(METHOD_PACKAGE[[tk$method]])   # package missing: not timed
message(sprintf("[%s] %s, p = %d: runs %d-%d on %s", format(Sys.time(), "%F %T"), tk$method, tk$p, min(runs), max(runs), Sys.info()[["nodename"]]))
if (!failed && is.null(eval(expr))) failed <- TRUE    # warm-up fit (untimed); NULL: runtime_safe() caught an error
for (k in seq_along(runs)) {
  if (failed) break
  set.seed(runs[k])
  t0 <- proc.time()[["elapsed"]]
  v  <- eval(expr)
  el <- proc.time()[["elapsed"]] - t0
  if (is.null(v)) { failed <- TRUE; break }          # runtime_safe() returned NULL: the method cannot be fitted here
  sec[k] <- el; rm(v)
}
out <- data.frame(Method = tk$method, p = tk$p, run = runs, sec = if (failed) NA_real_ else sec, failed = failed)
cpu <- tryCatch(sub(".*: ", "", grep("model name", readLines("/proc/cpuinfo"), value = TRUE)[1]), error = function(e) NA)
attr(out, "node") <- Sys.info()[["nodename"]]; attr(out, "cpu") <- cpu; attr(out, "sampler") <- SAMPLER_TAG
attr(out, "package_versions") <- package_versions()
save_rds_atomic(out, f)
message(sprintf("[%s] done: median %.1f s", format(Sys.time(), "%F %T"), if (failed) NA else median(sec)))
