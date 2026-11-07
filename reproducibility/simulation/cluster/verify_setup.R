###############################################################################
## verify_setup.R  --  checks of the installed environment (called by setup.sbatch)
###############################################################################
pk <- c("MR.RGM", "MendelianRandomization", "OneSampleMR", "mrbayes", "rjags", "MR2G", "MR2", "MrDAG", "grivet", "igraph", "pROC",
        "microbenchmark", "future", "ggplot2")
v  <- vapply(pk, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else NA_character_, "")
print(data.frame(package = pk, version = v), row.names = FALSE)
if (anyNA(v)) stop("missing: ", paste(pk[is.na(v)], collapse = ", "))
stopifnot(packageVersion("MR.RGM") >= "0.1.3", "gibbsB" %in% names(formals(MR.RGM::RGM)))
cat(R.version.string, "\n"); print(sessionInfo()$BLAS); print(extSoftVersion()[c("BLAS")])
cat("JAGS:", as.character(rjags::jags.version()), "\n")
## the simulation code loads and finds every method
setwd(Sys.getenv("SLURM_SUBMIT_DIR", "."))
source("sim_functions.R")
cat("Methods:", paste(ALL_METHODS, collapse = ", "), "\nSetup verified.\n")
