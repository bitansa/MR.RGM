###############################################################################
## smoke_add.R  --  quick check of the methods that are added to existing results (called by smoke_add.sbatch before
## the add array): one replicate of the exposure -> outcome design (n = 1,000) and of Case 2 (p = 5, n = 500), with
## the methods of the current method list that are missing from the result files (cluster/tasks_add.csv). Fails if a
## method stops with an error in the exposure -> outcome design or returns no metric.
###############################################################################
source("sim_functions.R")
ta  <- read.csv("cluster/tasks_add.csv", stringsAsFactors = FALSE)
new <- unique(unlist(strsplit(ta$missing, " ")))
message("Methods to add: ", paste(new, collapse = ", "))
ok <- TRUE
for (d in list(list(scenario = "roles", p = 8, n = 1000, experiment = "roles"), list(scenario = "smallworld", p = 5, n = 500, experiment = "main"))) {
  t0 <- proc.time()[["elapsed"]]
  r <- run_replicate(1, d$scenario, d$p, d$n, experiment = d$experiment, setting_label = "base",
                     methods = intersect(new, METHODS[[d$experiment]]))
  for (m in intersect(new, METHODS[[d$experiment]])) {
    x <- r[r$method == m, ]
    failed <- any(x$family == "diagnostic" & x$metric == "failed_fits" & x$value >= 1)
    auc <- x$value[x$family == "graph" & x$metric == "AUC"]
    message(sprintf("%s, %s (p = %d, n = %d): %s, network AUC %s, %.0f s", m, d$scenario, d$p, d$n,
                    if (failed) "FAILED" else "fitted", if (length(auc)) sprintf("%.3f", auc) else "--", proc.time()[["elapsed"]] - t0))
    if (!nrow(x) || (failed && d$scenario == "roles")) ok <- FALSE
  }
}
if (!ok) stop("smoke test of the added methods failed")
message("Smoke test of the added methods passed")
