###############################################################################
## check_results.R  --  completeness of the cluster run (run from reproducibility/simulation/; called by check.sbatch)
##   - every simulation task has its result file; per method: replicates with results and failed fits
##   - every runtime block is present
## Writes cluster/check_results.csv and prints the incomplete entries.
###############################################################################
source("sim_functions.R")
row_method <- function(m) ifelse(m %in% c("MR-IVW", "MR-SimpleMedian", "MR-WeightedMedian"), "MendelianRandomization", m)
ts <- read.csv("cluster/tasks_sim.csv", stringsAsFactors = FALSE)
ts$file <- ifelse(ts$script == "sim",
  file.path(paste0("results_", ts$experiment), sprintf("sim_%s_%s_p%d_n%g_%s.rds", ts$experiment, ts$scenario, ts$p, ts$n, ts$label)),
  file.path("results_strength", sprintf("strength_%s_%s_p%d_n%d.rds", ts$scenario, ts$strength, ts$p, ts$n)))
missing <- ts[!file.exists(ts$file), ]
cat(sprintf("Simulation tasks: %d, result files present: %d\n", nrow(ts), nrow(ts) - nrow(missing)))
if (nrow(missing)) { cat("Missing (rerun these task numbers):", paste(missing$task, collapse = ","), "\n") }
if (!any(file.exists(ts$file))) stop("no simulation result file found")
out <- do.call(rbind, lapply(which(file.exists(ts$file)), function(i) {
  x <- as.data.frame(readRDS(ts$file[i])); x <- x[x$method != "data", ]
  if (!identical(attr(readRDS(ts$file[i]), "sampler"), SAMPLER_TAG)) warning("no current sampler tag: ", ts$file[i])
  x$method <- row_method(x$method)
  fl <- x[x$family == "diagnostic" & x$metric == "failed_fits", ]
  target <- if (ts$script[i] == "sim") methods_for_setting(ts$experiment[i], ts$label[i]) else METHODS$strength
  do.call(rbind, lapply(target, function(m) data.frame(task = ts$task[i], experiment = ts$experiment[i], scenario = ts$scenario[i],
    p = ts$p[i], n = ts$n[i], setting = ts$label[i], method = m, replicates = length(unique(x$rep[x$method == m])),
    failed = length(unique(fl$rep[fl$method == m & is.na(fl$seconds) & fl$value >= 1])))))
}))
write.csv(out, "cluster/check_results.csv", row.names = FALSE)
short <- out[out$replicates < SIM_SETTINGS$n_rep, ]
cat(sprintf("Method x setting entries: %d; with fewer than %d replicates: %d\n", nrow(out), SIM_SETTINGS$n_rep, nrow(short)))
if (nrow(short)) print(short, row.names = FALSE)
fails <- aggregate(failed ~ experiment + method, data = out, FUN = sum)
cat("\nFailed fits (stopped with an error) per experiment and method:\n"); print(fails[fails$failed > 0, ], row.names = FALSE)

tr <- read.csv("cluster/tasks_runtime.csv", stringsAsFactors = FALSE)
tr$file <- file.path("Runtime Analysis", "cluster_parts", sprintf("rt_p%02d_m%02d_%02d.rds", tr$p, match(tr$method, names(RUNTIME_EXPRS)), tr$block))
rm_ <- tr[!file.exists(tr$file), ]
cat(sprintf("\nRuntime tasks: %d, present: %d\n", nrow(tr), nrow(tr) - nrow(rm_)))
if (nrow(rm_)) cat("Missing runtime tasks:", paste(rm_$task, collapse = ","), "\n")
if (!nrow(missing) && !nrow(short) && !nrow(rm_)) cat("\nEVERYTHING COMPLETE\n")
