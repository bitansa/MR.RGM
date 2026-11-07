###############################################################################
## make_add_tasks.R  --  task list of the SLURM array that adds methods to existing results (cluster/add_methods.R)
##   Rscript cluster/make_add_tasks.R        (from reproducibility/simulation/)
## One task per result file that lacks a method of the current method list (METHODS in sim_functions.R) for some
## replicate; largest settings first. Writes cluster/tasks_add.csv.
###############################################################################
source("sim_functions.R")
row_method <- function(m) ifelse(m %in% c("MR-IVW", "MR-SimpleMedian", "MR-WeightedMedian"), "MendelianRandomization", m)
files <- c(unlist(lapply(c("main", "roles", "noise", "confounding", "sensitivity", "largen"),
                         function(e) list.files(paste0("results_", e), "^sim_.*\\.rds$", full.names = TRUE))),
           list.files("results_strength", "^strength_.*\\.rds$", full.names = TRUE))
rows <- lapply(files, function(f) {
  df <- readRDS(f); strength <- grepl("^strength_", basename(f))
  exp <- if (strength) "strength" else sub("^results_", "", basename(dirname(f)))
  target <- if (strength) METHODS$strength else methods_for_setting(exp, df$setting_label[1])
  reps <- unique(df$rep)
  miss <- target[vapply(target, function(m) !all(reps %in% df$rep[row_method(df$method) == m]), logical(1))]
  if (!length(miss)) return(NULL)
  data.frame(file = f, experiment = exp, p = df$p[1], n = df$n[1], missing = paste(miss, collapse = " "))
})
ta <- do.call(rbind, rows)
if (is.null(ta)) { message("Every result file is complete."); quit(save = "no") }
ta <- ta[order(-ta$n, -ta$p), ]
ta$task <- seq_len(nrow(ta))
write.csv(ta[, c("task", "file", "experiment", "p", "n", "missing")], "cluster/tasks_add.csv", row.names = FALSE)
message(sprintf("%d result files to complete; methods to add: %s", nrow(ta), paste(unique(unlist(strsplit(ta$missing, " "))), collapse = ", ")))
