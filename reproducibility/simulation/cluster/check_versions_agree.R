###############################################################################
## check_versions_agree.R  --  confirms that the installed version of MR.RGM reproduces the stored results of the
## methods whose sampler is the same in versions 0.1.3 and 0.1.4 (MR.RGM and MR.RGM_NoConf)
##   Rscript cluster/check_versions_agree.R        (from reproducibility/simulation/)
## For three settings of the main study and three replicates each, the two methods are fitted again to the same data
## with the same seeds, and every metric is compared with the value stored in the result file. The script prints the
## largest absolute difference and stops if a stored value is not reproduced.
###############################################################################
source("sim_functions.R")
METH  <- c("MR.RGM", "MR.RGM_NoConf")
CASES <- list(list(scenario = "smallworld",       p = 10, n = 1000),
              list(scenario = "scalefree",        p = 5,  n = 500),
              list(scenario = "smallworld_pleio", p = 10, n = 10000))
REPS  <- 1:3
TOL   <- 1e-10
key   <- c("method", "family", "metric")
worst <- 0; worst_diag <- 0
for (cs in CASES) {
  f   <- file.path("results_main", sprintf("sim_main_%s_p%d_n%g_base.rds", cs$scenario, cs$p, cs$n))
  old <- readRDS(f)
  message(sprintf("%s: stored results fitted with MR.RGM %s; fitted again with MR.RGM %s", f,
                  attr(old, "package_versions")[["MR.RGM"]], as.character(packageVersion("MR.RGM"))))
  old <- as.data.frame(old)
  for (r in REPS) {
    new <- as.data.frame(run_replicate(r, cs$scenario, cs$p, cs$n, experiment = "main", setting_label = "base", methods = METH))
    m <- merge(old[old$rep == r & old$method %in% METH, c(key, "value")], new[new$method %in% METH, c(key, "value")],
               by = key, suffixes = c(".stored", ".now"))
    d <- abs(m$value.stored - m$value.now)
    d[is.na(m$value.stored) & is.na(m$value.now)] <- 0          # undefined in both
    d[is.na(d)] <- Inf                                          # defined in one only
    diag <- m$family == "diagnostic"
    worst <- max(worst, d[!diag]); worst_diag <- max(worst_diag, d[diag])
    message(sprintf("   replicate %d: %d metrics compared; largest difference %.3g (diagnostics %.3g)",
                    r, nrow(m), max(d[!diag]), if (any(diag)) max(d[diag]) else 0))
    if (any(d[!diag] > TOL)) print(m[!diag & d > TOL, ], row.names = FALSE)
  }
}
message(sprintf("Largest absolute difference over all settings: %.3g (diagnostics: %.3g)", worst, worst_diag))
if (worst > TOL) stop("the installed version does not reproduce the stored results of MR.RGM / MR.RGM_NoConf")
message("The installed version reproduces the stored results of MR.RGM and MR.RGM_NoConf.")
