###############################################################################
## pilot_mrdag_scaling.R  --  pilot check of the input scale of MrDAG (not part of the reported study)
##
## The study multiplies the per-instrument marginal effects by sqrt(n) before they are passed to MrDAG, whose default
## prior assumes input on the unit scale. This script runs MrDAG on the exposure -> outcome design with and without
## this scaling, on the replicates of the study (same truth and data seeds, standardized data; sim_functions.R), and
## reports the network AUC, the number of selected edges and the effect error of both versions.
##
## Output: pilot_mrdag_scaling/pilot_results.csv (one row per replicate, version and metric) and a printed summary.
## Usage (this folder):  source("pilot_mrdag_scaling.R")
## Options before sourcing: PILOT_REPS (default 1:5), PILOT_CORES (default: all cores but one).
###############################################################################

source("sim_functions.R")
if (!has_pkg("MrDAG")) stop("The package MrDAG is not installed.")
REPS  <- if (exists("PILOT_REPS"))  PILOT_REPS  else 1:5
CORES <- if (exists("PILOT_CORES")) PILOT_CORES else max(1, parallel::detectCores() - 1)
OUT   <- "pilot_mrdag_scaling"; dir.create(OUT, showWarnings = FALSE)

## MrDAG on the exposure -> outcome split, with or without the sqrt(n) input scaling; estimates as [to, from]
fit_mrdag_scale <- function(ms, D, roles, scale, set = SIM_SETTINGS) {
  p <- ncol(ms$beta); Xi <- roles$X_idx; Yi <- roles$Y_idx
  ivs <- which(colSums(D[Xi, , drop = FALSE]) > 0)
  dat <- (if (scale) sqrt(ms$n) else 1) * ms$beta[ivs, c(Yi, Xi), drop = FALSE]
  chk <- list(Y_idx = seq_along(Yi), X_idx = length(Yi) + seq_along(Xi))
  out <- MrDAG::MrDAG(data = dat, niter = set$mrdag$niter, burnin = set$mrdag$burnin, thin = set$mrdag$thin,
                      pp = set$mrdag$pp, tempMax = set$mrdag$tempMax, MrDAGcheck = chk, fileName = NULL,
                      seed = sample.int(.Machine$integer.max, 1))
  peip <- as.matrix(MrDAG::get_edgeprob(out)$edgeProb)
  ce   <- as.matrix(MrDAG::get_causaleffects(out, progress = FALSE)$causalEffects)
  ord  <- c(Yi, Xi)
  score <- matrix(NA_real_, p, p); score[ord, ord] <- t(peip)
  est <- matrix(NA_real_, p, p); est[Yi, Xi] <- t(ce)
  list(est = est, score = score, edge = !is.na(score) & score > set$pip_thr)
}

## One replicate of the exposure -> outcome design: the data of the study, then both versions of the input
one_rep <- function(job) {
  rep <- job$rep; n <- job$n; scenario <- "roles"; set <- SIM_SETTINGS
  set.seed(truth_seed(rep)); truth0 <- gen_truth(scenario, 8, set)
  set.seed(data_seed(scenario, 8, n, rep, set)); dat0 <- gen_data(truth0, n, set)
  std <- standardize_data(dat0, truth0, set); truth <- std$truth
  ms  <- marginal_stats(std$X, std$Y); attr(ms, "X") <- std$X
  rows <- list()
  for (sc in c(TRUE, FALSE)) {
    set.seed(data_seed(scenario, 8, n, rep, set) + METHOD_SEED[["MrDAG"]])        # the seed of MrDAG in the study
    t0 <- proc.time()[["elapsed"]]
    f  <- tryCatch(fit_mrdag_scale(ms, truth$D, truth$roles, sc, set), error = function(e) e)
    m  <- if (sc) "MrDAG, sqrt(n) scaling (study)" else "MrDAG, unscaled input"
    if (inherits(f, "error")) {
      rows[[length(rows) + 1]] <- data.frame(n = n, rep = rep, method = m, metric = "error", value = NA_real_,
                                             seconds = NA_real_, message = conditionMessage(f))
    } else {
      vals <- c(graph_metrics(truth$graph, f$edge, f$score), edges = sum(f$edge, na.rm = TRUE),
                effect_metrics(truth$A, f$est)["MeanAbsDev"])
      rows[[length(rows) + 1]] <- data.frame(n = n, rep = rep, method = m, metric = names(vals), value = unname(vals),
                                             seconds = proc.time()[["elapsed"]] - t0, message = NA_character_)
    }
  }
  do.call(rbind, rows)
}

jobs <- list()
for (n in c(1000, 10000)) for (r in REPS) jobs[[length(jobs) + 1]] <- list(n = n, rep = r)
message(sprintf("Pilot: %d replicates on %d cores", length(jobs), CORES))
res <- run_parallel(jobs, one_rep, cores = CORES, timeout_min = 180, label = "pilot")
ok  <- !vapply(res, inherits, logical(1), "try-error")
if (any(!ok)) message(sum(!ok), " replicate(s) failed: ", paste(unique(unlist(res[!ok])), collapse = " | "))
out <- do.call(rbind, res[ok])
write.csv(out, file.path(OUT, "pilot_results.csv"), row.names = FALSE)

## Summary: mean over the replicates
summ <- aggregate(value ~ n + method + metric, data = out[out$metric != "error", ], FUN = function(v) round(mean(v), 3))
wide <- reshape(summ, idvar = c("n", "method"), timevar = "metric", direction = "wide")
names(wide) <- sub("^value\\.", "", names(wide))
print(wide[order(wide$n, wide$method), intersect(c("n", "method", "AUC", "TPR", "FDR", "edges", "MeanAbsDev"), names(wide))], row.names = FALSE)
if (any(out$metric == "error")) print(unique(out[out$metric == "error", c("n", "method", "message")]), row.names = FALSE)
message("Saved ", file.path(OUT, "pilot_results.csv"))
