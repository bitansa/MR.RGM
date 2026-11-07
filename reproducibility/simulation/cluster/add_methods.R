###############################################################################
## add_methods.R  --  completes existing result files with methods that were added to the study later
##
## For a saved setting (results_<experiment>/sim_*.rds or results_strength/strength_*.rds) this script runs only the
## methods of the current method list (METHODS in sim_functions.R) that the file does not yet contain, for the
## replicates that lack them, on the same data and with the same seeds, and merges them into the result file and its
## estimates file (fits_*.rds). Every method is fitted after set.seed(data seed + method offset) (run_replicate()), so
## its result does not depend on which other methods run in the same replicate: a completed file is what
## run_simulations.R / run_instrument_strength.R produce with the current method list.
##
## Safety: the original files are copied once to  _Archive before added methods/  (unchanged); files are replaced
## atomically; the settings stored with a file must equal SIM_SETTINGS for every key that the existing methods use.
## Restartable: methods and replicates already present are skipped.
##
## Usage (from reproducibility/simulation/):
##   ADD_FILES <- "results_main/sim_main_smallworld_p10_n10000_base.rds"; source("cluster/add_methods.R")
##   (cluster: one result file per task, cluster/run_add_task.R; without ADD_FILES every result file is completed)
## Optional: N_CORES, REP_TIMEOUT_MIN (per replicate; four times as long at n >= 10^6), REP_RETRIES (further attempts
## of a replicate that fails or exceeds the limit; default 1), REFIT_FAILED (methods whose failed fits are fitted again).
###############################################################################

source("sim_functions.R")
if (!exists("N_CORES"))         N_CORES         <- max(1, parallel::detectCores() - 1)
if (!exists("REP_TIMEOUT_MIN")) REP_TIMEOUT_MIN <- 600
if (!exists("REP_RETRIES"))     REP_RETRIES     <- 1
## REFIT_FAILED: methods whose failed fits (stopped with an error; failed_fits = 1 without a time) are fitted again,
## for example after a correction of the method's code (default: none)
if (!exists("REFIT_FAILED"))    REFIT_FAILED    <- character(0)
BACKUP <- "_Archive before added methods"
stamp  <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))

## method name of a result row (the three estimators of MendelianRandomization are stored as MR-IVW, MR-SimpleMedian,
## MR-WeightedMedian)
row_method <- function(m) ifelse(m %in% c("MR-IVW", "MR-SimpleMedian", "MR-WeightedMedian"), "MendelianRandomization", m)
backup <- function(f) {
  dst <- file.path(BACKUP, f)
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(dst)) stopifnot(file.copy(f, dst, copy.date = TRUE))
}
## the settings stored with a result file must equal the current SIM_SETTINGS for every key of the methods it contains
## (a key that is new, such as the settings of an added method, is taken from SIM_SETTINGS)
check_settings <- function(st, f) {
  for (k in intersect(c("mcmc", "mr_boot", "mrbayes", "mr2g", "mr2", "mrdag", "alpha", "pip_thr", "n_rep"), names(st)))
    if (!identical(st[[k]], SIM_SETTINGS[[k]])) stop("setting '", k, "' of ", f, " differs from SIM_SETTINGS")
  for (k in setdiff(names(SIM_SETTINGS), names(st))) st[[k]] <- SIM_SETTINGS[[k]]
  st
}

SIM_COLS <- c("experiment", "scenario", "p", "n", "rep", "setting", "setting_label", "method", "family", "metric", "value", "seconds")
ST_COLS  <- c("scenario", "p", "n", "strength", "rep", "method", "family", "metric", "value", "seconds")

complete_file <- function(f) {
  strength <- grepl("^strength_", basename(f))
  exp <- if (strength) "strength" else sub("^results_", "", basename(dirname(f)))
  ff  <- file.path(dirname(f), if (strength) sub("^strength_", "fits_strength_", basename(f)) else sub("^sim_", "fits_", basename(f)))
  df  <- readRDS(f)
  if (!identical(attr(df, "sampler"), SAMPLER_TAG)) stop("result file without the current sampler tag: ", f)
  st <- check_settings(attr(df, "settings"), f)
  if (!identical(names(df), if (strength) ST_COLS else SIM_COLS)) stop("unexpected columns in ", f)
  if (!file.exists(ff)) stop("no estimates file ", ff)
  fo   <- readRDS(ff)
  reps <- sort(unique(df$rep))
  if (!setequal(vapply(fo$fits, function(x) as.integer(x$rep), integer(1)), reps)) stop("replicates of ", ff, " and ", f, " differ")
  target <- if (strength) METHODS$strength else methods_for_setting(exp, df$setting_label[1])
  scen <- df$scenario[1]; p <- df$p[1]; n <- df$n[1]
  if (is.null(attr(df, "package_versions_original"))) attr(df, "package_versions_original") <- attr(df, "package_versions")
  for (m in target) {
    if (m %in% REFIT_FAILED) {                                    # drop the failed fits of m, so that they are fitted again
      fail <- unique(df$rep[row_method(df$method) == m & df$family == "diagnostic" & df$metric == "failed_fits" &
                            !is.na(df$value) & df$value >= 1 & is.na(df$seconds)])
      if (length(fail)) {
        backup(f); backup(ff)
        stamp(sprintf("%s: %s  refitting %d failed fit(s) of %s: replicates %s", exp, basename(f), length(fail), m,
                      paste(fail, collapse = ", ")))
        keep_a <- attributes(df)[setdiff(names(attributes(df)), c("names", "row.names", "class"))]
        df <- df[!(row_method(df$method) == m & df$rep %in% fail), , drop = FALSE]; rownames(df) <- NULL
        for (a in names(keep_a)) attr(df, a) <- keep_a[[a]]
      }
    }
    rr <- reps[!vapply(reps, function(r) m %in% row_method(df$method[df$rep == r]), logical(1))]
    if (!length(rr)) next
    backup(f); backup(ff)
    stamp(sprintf("%s: %s  adding %s  (%d replicates on %d cores)", exp, basename(f), m, length(rr), N_CORES))
    t0 <- Sys.time()
    tmo <- if (n >= 1e6) 4 * REP_TIMEOUT_MIN else REP_TIMEOUT_MIN
    res <- if (strength) {
      run_parallel(rr, strength_replicate, scenario = scen, p = p, n = n, strength = df$strength[1], set = st, methods = m,
                   cores = N_CORES, timeout_min = tmo, retries = REP_RETRIES, label = paste(exp, basename(f), m))
    } else {
      run_parallel(rr, run_replicate, scenario = scen, p = p, n = n, set = st, experiment = exp, setting = df$setting[1],
                   setting_label = df$setting_label[1], methods = m, cores = N_CORES, timeout_min = tmo, retries = REP_RETRIES,
                   label = paste(exp, basename(f), m))
    }
    bad <- vapply(res, function(r) inherits(r, "try-error") || is.null(r), logical(1))
    if (any(bad)) message(sprintf("   %d replicate(s) did not return (rerun this task to retry them): %s",
                                  sum(bad), paste(rr[bad], collapse = ", ")))
    message(sprintf("   done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    if (all(bad)) next
    ## estimates first (merged by replicate and method name; idempotent), then the metrics (as in run_simulations.R)
    fr <- vapply(fo$fits, function(x) as.integer(x$rep), integer(1))
    add <- list()
    for (r in res[!bad]) {
      fi <- attr(r, "fits"); attr(r, "fits") <- NULL
      i  <- which(fr == as.integer(fi$rep))
      fo$fits[[i]]$methods[names(fi$methods)] <- fi$methods
      add[[length(add) + 1]] <- as.data.frame(r[r$method != "data", , drop = FALSE])   # the "data" rows are in the file
    }
    add <- do.call(rbind, add)
    done_m <- intersect(target, unique(c(row_method(df$method), row_method(add$method))))
    fo$methods <- done_m
    save_rds_atomic(fo, ff)
    keep <- attributes(df)[setdiff(names(attributes(df)), c("names", "row.names", "class"))]
    out  <- rbind(as.data.frame(df), add[, names(df)])
    out  <- out[order(out$rep), , drop = FALSE]; rownames(out) <- NULL
    for (a in names(keep)) attr(out, a) <- keep[[a]]
    attr(out, "methods") <- done_m
    attr(out, "settings") <- st
    am <- attr(out, "added_methods")
    attr(out, "added_methods") <- list(methods = union(am$methods, m), session_info = utils::sessionInfo(),
                                       package_versions = package_versions())
    save_rds_atomic(out, f)
    df <- out
  }
  invisible(TRUE)
}

files <- if (exists("ADD_FILES")) ADD_FILES else
  c(unlist(lapply(c("main", "roles", "noise", "confounding", "sensitivity", "largen"),
                  function(e) list.files(paste0("results_", e), "^sim_.*\\.rds$", full.names = TRUE))),
    list.files("results_strength", "^strength_.*\\.rds$", full.names = TRUE))
for (f in files) {
  options(sim_error_log = file.path(dirname(f), "errors.log"))
  tryCatch(complete_file(f), error = function(e) { message("   FAILED (", f, "): ", conditionMessage(e)); if (exists("ADD_FILES")) stop(e) })
}
stamp("add_methods.R finished: ", length(files), " result file(s)")
