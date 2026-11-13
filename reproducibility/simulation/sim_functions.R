###############################################################################
## sim_functions.R
##
## Shared code for the simulation study of MR.RGM:
##   - simulation settings (single source of truth, reported in the paper)
##   - generation of the true network and of the data; standardization of the data
##   - the MR.RGM family: MR.RGM (known instrument map, no selection of the instrument effects),
##     MR.RGM_NoConf (same, independent errors), MR.RGM+ (D = 1, instrument effects selected)
##   - competing methods, each run with its package's documented defaults:
##       MR2G                     causal network with cycles from GWAS-type summary statistics (Lin, Pan & Xue 2026)
##       MVMR                     multivariable IVW, one outcome at a time (MendelianRandomization::mr_mvivw)
##       MR2                      multi-response MR with residual correlation (Zuber et al. 2023)
##       MrDAG                    Bayesian DAG MR for exposures -> outcomes (Zuber et al. 2025)
##       GrIVET                   causal DAG with hidden confounding from instruments (Chen, Li, Shen and Pan 2024);
##                                fitted to the unstandardized data (its tuning is not scale free; see fit_grivet())
##                                (MR2, MrDAG: with the exposure/outcome split of the "roles" design. In the designs without
##                                roles MrDAG is applied one outcome at a time with the other traits as exposures, as MVMR
##                                (by_outcome()); MR2 cannot be (its sampler conditions each response on the others, so it
##                                needs at least two responses) and is evaluated in the roles design only)
##       MR-IVW, median, mrbayes, OneSampleMR   pairwise baselines (every ordered pair of traits)
##     Every method is run in every experiment in which it can be applied (METHODS): MR2 in the roles design only; the
##     sensitivity experiment varies the MR.RGM priors and runs the MR.RGM family only (the competitors do not use these
##     priors and see the data of the main study).
##   - performance metrics
##
## Design (one replicate):
##   1. set.seed(truth seed)  -> draw graph, A and the confounder loadings C        (gen_truth)
##   2. set.seed(data seed)   -> draw U, E and X, compute Y; then standardize       (gen_data, standardize_data)
##   3. every method is fitted to the same standardized data set, each after set.seed(data seed + method offset),
##      so that no method's random stream depends on which other methods were run
##   The truth seed depends only on the replicate, so replicate r has the same true network at every n;
##   the data seed is unique to (scenario, p, n, r) and shared by the settings of an experiment (paired comparison).
##
## Experiments (run_simulations.R):
##   main         Cases 1-3 (scale-free; small-world; small-world with pleiotropy), p in {5, 10}, n in {500, ..., 30000}
##   noise        error variance 1, 4, 9 (small-world, p = 10)
##   confounding  confounder loadings scaled by 0, 0.5, 1, 1.5 (small-world, p = 10)
##   largen       n = 1e5, 1e6 (small-world, p = 10)
##   longchain    MR.RGM+ at n = 1e6 with a chain of 250,000 iterations (the data and seeds of largen)
##   sensitivity  prior inclusion Beta(3,1) vs Beta(1,1); nu1 in {1e-2, 1e-3, 1e-4} (MR.RGM family)
##   roles        exposures -> outcomes design with known roles (MR2 and MrDAG with their native exposure/outcome split)
##   strength     instrument strength (run_instrument_strength.R; strength_replicate() below)
###############################################################################

suppressPackageStartupMessages({
  library(igraph)
  library(pROC)
  library(MendelianRandomization)
  library(OneSampleMR)
  library(mrbayes)
  library(rjags)
  library(MR.RGM)
})
## Competitor packages that are not on CRAN are loaded when first used (see has_pkg()); install with
##   remotes::install_github("ZhaotongL/MR2G")
##   install.packages(c("MR2", "MrDAG"), repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
##   (or remotes::install_github("lb664/MR2"); remotes::install_github("lb664/MrDAG"))
##   remotes::install_github("chunlinli/grivet", subdir = "grivet")       # GrIVET (package name grivet)
has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)
## R package of each method that is not part of the MR.RGM package, and the versions recorded with every result file
METHOD_PACKAGE <- c(MR2G = "MR2G", MR2 = "MR2", MrDAG = "MrDAG", GrIVET = "grivet", MendelianRandomization = "MendelianRandomization",
                    OneSampleMR = "OneSampleMR", mrbayes = "mrbayes")
STUDY_PACKAGES <- c("MR.RGM", "MendelianRandomization", "OneSampleMR", "mrbayes", "rjags", "MR2G", "MR2", "MrDAG", "grivet")
package_versions <- function(pk = STUDY_PACKAGES)
  vapply(pk, function(x) if (has_pkg(x)) as.character(packageVersion(x)) else NA_character_, character(1))
## loaded once here (not first inside a seeded fit), so that every fit sees the same state whatever ran before it
invisible(vapply(c("MR2G", "MR2", "MrDAG", "grivet", "future"), has_pkg, logical(1)))

## ---------------------------------------------------------------------------
## Displayed text (tables, captions, figure labels). The result files, data frames and lists keep the method keys
## "MR.RGM", "MR.RGM+", "MR.RGM_NoConf" and the scenario keys (e.g. "smallworld_pleio"); only the text that is shown is
## converted: the methods are displayed as IV-RGM, IV-RGM+ and IV-RGM_NoConf, (horizontal) pleiotropy as invalid
## instruments, and trait(s) as response(s). References to the R package MR.RGM ("package MR.RGM", "MR.RGM >= 0.1.4",
## "MR.RGM::", ...) are kept. Character vector in, character vector out (names kept); a second application changes nothing.
## real_data_functions.R defines the same function without renaming SNPs, which the real-data tables keep.
## ---------------------------------------------------------------------------
display_text <- function(x) {
  if (!length(x)) return(x)
  nm <- names(x)
  if (!is.character(x)) x <- as.character(x)
  keep  <- c("package \\texttt{MR.RGM}", "\\texttt{MR.RGM} package", "package MR.RGM", "MR.RGM package",
             "MR.RGM 0.", "MR.RGM (version", "MR.RGM >=", "MR.RGM::", "library(MR.RGM)",
             "quantitative trait")     # package references; the genetics term "quantitative trait loci"
  token <- paste0("\001", seq_along(keep), "\002")
  for (i in seq_along(keep)) x <- gsub(keep[i], token[i], x, fixed = TRUE)
  fixed_map <- c("MR.RGM\\_NoConf" = "IV-RGM\\_NoConf", "MR.RGM_NoConf" = "IV-RGM_NoConf", "MR.RGM+" = "IV-RGM+",
                 "MR.RGM" = "IV-RGM",
                 "Horizontal pleiotropy" = "Invalid instruments", "horizontal pleiotropy" = "invalid instruments",
                 "No pleiotropy" = "Valid instruments", "no pleiotropy" = "valid instruments",
                 "Pleiotropy" = "Invalid instruments", "pleiotropy" = "invalid instruments")
  for (i in seq_along(fixed_map)) x <- gsub(names(fixed_map)[i], fixed_map[[i]], x, fixed = TRUE)
  word_map <- c("Traits" = "Responses", "traits" = "responses", "Trait" = "Response", "trait" = "response",
                "SNPs" = "instruments", "SNP" = "instrument")      # simulated instruments (the real-data tables keep "SNP")
  for (i in seq_along(word_map)) x <- gsub(paste0("\\b", names(word_map)[i], "\\b"), word_map[[i]], x, perl = TRUE)
  for (i in seq_along(keep)) x <- gsub(token[i], keep[i], x, fixed = TRUE)
  names(x) <- nm
  x
}

## ---------------------------------------------------------------------------
## Outcome of every fit of a method that may stop (used for GrIVET). A fit ends in one of four ways: a result; "peeling"
## (GrIVET's peeling step could not order every response); "numerical" (any other error, e.g. NA/NaN/Inf in a
## numerical routine); "time" (the fit did not finish within the time limit of its cluster job). fail_causes() reads
## the error logs of an experiment and returns the cause of each logged failure, using the most recent message of a
## fit. The coercion error of R >= 4.3 is removed by grivet_compat() below; a fit still recorded as failed whose most
## recent message is that error was therefore stopped by the time limit, and its cause is "time".
## add_unfinished() adds a failed-fit row (cause "time") for every replicate in which the method has no row at all,
## i.e. the fit was stopped by the time limit before it returned; the template rows are the diagnostic rows of MR.RGM.
## The added rows are returned in attr(, "unfinished"). Messages of forked workers that wrote to the log at the same
## moment can share a line; such lines are split at every time stamp before they are parsed.
## ---------------------------------------------------------------------------
fail_causes <- function(log_files, method = "GrIVET") {
  empty <- data.frame(experiment = character(0), scenario = character(0), p = numeric(0), n = numeric(0),
                      setting_label = character(0), rep = integer(0), cause = character(0), stringsAsFactors = FALSE)
  log_files <- log_files[file.exists(log_files)]
  if (!length(log_files)) return(empty)
  x <- unlist(lapply(log_files, readLines, warn = FALSE), use.names = FALSE)
  x <- gsub("(.)(\\[[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}\\])", "\\1\n\\2", x, perl = TRUE)
  x <- unlist(strsplit(x, "\n", fixed = TRUE), use.names = FALSE)
  x <- trimws(x[nzchar(trimws(x))])
  pat <- sprintf("^\\[([^]]*)\\] (\\S+) (\\S+) p=([0-9]+) n=(\\S+) setting=(\\S+) rep=([0-9]+): %s failed: (.*)$", method)
  x <- x[grepl(pat, x, perl = TRUE)]
  if (!length(x)) return(empty)
  m <- regmatches(x, regexec(pat, x, perl = TRUE))
  d <- data.frame(time = vapply(m, `[`, "", 2), experiment = vapply(m, `[`, "", 3), scenario = vapply(m, `[`, "", 4),
                  p = as.numeric(vapply(m, `[`, "", 5)), n = as.numeric(vapply(m, `[`, "", 6)),
                  setting_label = vapply(m, `[`, "", 7), rep = as.integer(vapply(m, `[`, "", 8)),
                  msg = vapply(m, `[`, "", 9), stringsAsFactors = FALSE)
  d <- d[order(d$time), , drop = FALSE]
  key <- paste(d$experiment, d$scenario, d$p, d$n, d$setting_label, d$rep, sep = "|")
  d <- d[!duplicated(key, fromLast = TRUE), , drop = FALSE]             # most recent message of each fit
  d$cause <- ifelse(grepl("peeling", d$msg), "peeling",
             ifelse(grepl("coercion to 'logical(1)'", d$msg, fixed = TRUE), "time", "numerical"))
  rownames(d) <- NULL
  d[, names(empty)]
}
add_unfinished <- function(res, method, group_cols, template_method = "MR.RGM") {
  tmpl <- res[res$method == template_method & res$family == "diagnostic" & res$metric == "failed_fits", , drop = FALSE]
  if (!nrow(tmpl)) return(res)
  g_t <- do.call(paste, c(tmpl[group_cols], sep = "|"))
  g_m <- do.call(paste, c(res[res$method == method, group_cols, drop = FALSE], sep = "|"))
  have <- paste(g_m, res$rep[res$method == method], sep = "#")
  add <- tmpl[!(paste(g_t, tmpl$rep, sep = "#") %in% have), , drop = FALSE]
  if (!nrow(add)) return(res)
  add$method <- method; add$value <- 1; add$seconds <- NA_real_
  message(sprintf("%s: %d fit(s) without a result (stopped by the time limit) recorded as failed (cause \"time\")",
                  method, nrow(add)))
  out <- rbind(res, add); attr(out, "unfinished") <- add
  out
}

## ---------------------------------------------------------------------------
## Parallel replicates with a time limit.
## Each replicate runs in its own forked process (parallel::mcparallel), at most `cores` at a time. A replicate that
## has not returned after `timeout_min` minutes is terminated and started again (same seeds, so the same result if it
## completes); after `retries` failed attempts it is recorded as failed, so that no single replicate can block a
## setting. Results are returned in the order of X; failed replicates are "try-error" objects.
## ---------------------------------------------------------------------------
run_parallel <- function(X, FUN, ..., cores = 1, timeout_min = 60, retries = 1, label = "") {
  ## Windows (no forking) or a single core: sequential, without the watchdog
  if (.Platform$OS.type == "windows" || cores <= 1)
    return(lapply(X, function(x) tryCatch(FUN(x, ...), error = function(e) structure(conditionMessage(e), class = "try-error"))))
  res <- vector("list", length(X)); queue <- seq_along(X); tries <- integer(length(X)); running <- list()
  while (length(queue) || length(running)) {
    while (length(running) < cores && length(queue)) {
      i <- queue[1]; queue <- queue[-1]; tries[i] <- tries[i] + 1
      job <- parallel::mcparallel(FUN(X[[i]], ...), silent = FALSE)
      running[[as.character(job$pid)]] <- list(i = i, start = Sys.time(), job = job)
    }
    out <- parallel::mccollect(lapply(running, `[[`, "job"), wait = FALSE, timeout = 2)
    if (!is.null(out)) for (pid in names(out)) {
      r <- running[[pid]]; if (is.null(r)) next
      res[r$i] <- list(out[[pid]]); running[[pid]] <- NULL
    }
    for (pid in names(running)) {
      r <- running[[pid]]
      if (as.numeric(difftime(Sys.time(), r$start, units = "mins")) > timeout_min) {
        tools::pskill(as.integer(pid), tools::SIGKILL); try(parallel::mccollect(r$job, wait = FALSE), silent = TRUE)
        running[[pid]] <- NULL
        msg <- sprintf("[%s] %s replicate %d killed after %g min (attempt %d)", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                       label, X[[r$i]], timeout_min, tries[r$i])
        message(msg); lf <- getOption("sim_error_log"); if (!is.null(lf)) cat(msg, "\n", file = lf, append = TRUE, sep = "")
        if (tries[r$i] <= retries) queue <- c(queue, r$i) else res[r$i] <- list(structure(msg, class = "try-error"))
      }
    }
  }
  res
}

## These scripts require MR.RGM >= 0.1.4 (MR.RGM = known map without selection; MR.RGM+ = selection of the instrument
## effects with a point mass at zero). MR.RGM and MR.RGM_NoConf give the same results under 0.1.3 and 0.1.4.
if (packageVersion("MR.RGM") < "0.1.4")
  stop("MR.RGM ", packageVersion("MR.RGM"), " is installed; these scripts require MR.RGM >= 0.1.4: ",
       "remotes::install_github(\"bitansa/MR.RGM\")")
## Every result file records the sampler it was produced with; a result file without the current tag is not reused:
## the scripts stop and ask for it to be moved out of results_*/.
SAMPLER_TAG <- "MR.RGM 0.1.3 final sampler: Gibbs B, adaptive random-walk A, least-squares start"   # identifier of the sampler; compared with identical() against the stored results
result_is_current <- function(f) identical(attr(readRDS(f), "sampler"), SAMPLER_TAG)
## The tag identifies the sampler of MR.RGM and MR.RGM_NoConf, which is the same in versions 0.1.3 and 0.1.4 of the
## package. The sampler of MR.RGM+ changed in 0.1.4 (point mass at zero for the instrument effects), so its results
## must come from version 0.1.4 or later: either the file was produced with that version, or MR.RGM+ was refitted into
## it by cluster/add_methods.R (attribute "replaced_methods"). rgm_plus_is_current() checks this for a result object.
rgm_plus_is_current <- function(x) {
  if (!"MR.RGM+" %in% unique(x$method)) return(TRUE)
  v <- attr(x, "replaced_methods")$version[["MR.RGM+"]]
  v <- if (is.null(v)) attr(x, "package_versions")["MR.RGM"] else sub("^MR\\.RGM ", "", v)
  length(v) == 1 && !is.na(v) && package_version(as.character(v)) >= "0.1.4"
}
## Results are written to a temporary file and then renamed, so that a file is either complete or absent (several
## processes, e.g. cluster tasks, may write into the same folder)
save_rds_atomic <- function(obj, f) {
  tmp <- sprintf("%s.tmp%d", f, Sys.getpid()); saveRDS(obj, tmp)
  if (!file.rename(tmp, f)) stop("could not write ", f)
  invisible(f)
}
check_existing_result <- function(f) {
  if (!file.exists(f)) return(FALSE)
  if (!result_is_current(f)) stop("Result file without the current sampler tag: ", f,
                                  ". Move it out of results_*/ before running the study again.")
  TRUE
}


## ---------------------------------------------------------------------------
## Settings
## ---------------------------------------------------------------------------
SIM_SETTINGS <- list(
  scenarios     = c("scalefree", "smallworld", "smallworld_pleio", "roles"),
  p_values      = c(5, 10, 8),           # 8 = experiment "roles" (4 exposures + 4 outcomes)
  n_values      = c(500, 1000, 10000, 30000, 1e5, 1e6),
  n_rep         = 100,               # replicates per setting, in every experiment and in the instrument-strength study
  ivs_per_trait = 3,                 # own instruments per trait
  iv_effect     = 1,                 # effect of every instrument on the generating scale (1 in the main study; smaller in run_instrument_strength.R)
  l_star        = 5,                 # latent confounders U (each N(0,1))
  conf_scale    = 1,                 # multiplier of the confounder loadings (+-1 * conf_scale; 0 = no confounding)
  err_var       = 9,                 # error variance: E ~ N(0, err_var I_p)
  a_values      = c(-0.1, 0.1),      # non-zero entries of A on the generating scale
  sw_nei        = 2,                 # small-world: neighbourhood size
  sw_rewire     = 0.1,               # small-world: rewiring probability
  roles         = list(p_x = 4, p_y = 4, a_values = c(-0.3, 0.3)),   # experiment "roles": 4 exposures -> 4 outcomes (DAG)
  standardize   = TRUE,              # center and scale X and Y before fitting (every method); truth rescaled accordingly
  ## Every MCMC method (MR.RGM family, mrbayes, MR2, MrDAG) runs 50,000 iterations in total, of which 10,000 are burn-in;
  ## thinning 10 where the package has a thinning argument (mrbayes has none). Every other tuning parameter is the package default.
  mcmc          = list(nIter = 50000, nBurnin = 10000, Thin = 10),   # MR.RGM family
  rgm_hyper     = list(aRho = 3, bRho = 1, nu1 = 0.001),   # package defaults (experiment "sensitivity" varies them)
  mr_boot       = 10000,             # MendelianRandomization: bootstrap draws for median SEs (package default)
  mrbayes       = list(n.chains = 3, n.burn = 10000, n.iter = 40000),   # 3 chains (package default); per chain 10,000 burn-in + 40,000 = 50,000
  mr2g          = list(num_pert = 100, num_pert_M = 100, IV_screen = "ScreenAug"),   # MR2G() defaults
  mr2           = list(niter = 50000, burnin = 10000, thin = 10, EVgamma = c(2, 2)),  # MR2: chain length as above (the package has no default);
                                     # EVgamma = c(2, 2) is the package default (prior mean and variance of the number of exposures per
                                     # response); with q exposures it defines a beta-binomial prior that is proper only for q >= 4
                                     # (q = 4: Beta(1, 1) on the inclusion probability; q = 3: aPi = bPi = 0; q <= 2: mean >= q).
                                     # Hence 4 exposures in the roles design
  mrdag         = list(niter = 50000, burnin = 10000, thin = 10, pp = 0.05, tempMax = 10),   # MrDAG: chain length as above (no default); pp, tempMax = package defaults
  ## GrIVET: the package has no default tuning parameters; its truncated-L1 regressions are cross-validated over these
  ## grids (5 folds, as in its own scripts). The grid of its scripts (tau 0.2-0.3, gamma 0.1-1) is tuned to effects of
  ## order 1 and selects nothing on the data of this study (pilot runs); the grid below extends it downwards so
  ## that the effects of this study can be selected, and contains a coarse version of the original one.
  grivet        = list(tau = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.3), gamma = c(1e-4, 3e-4, 1e-3, 3e-3, 0.01, 0.03, 0.1, 0.3, 1),
                       n_fold = 5),
  alpha         = 0.05,              # competitors: edge if p-value < alpha (no multiplicity correction)
  pip_thr       = 0.5                # Bayesian methods: edge / confounding link if PIP > pip_thr
)

## Methods run in each experiment: every method in every experiment in which it can be applied.
##   MR2 needs an exposure/outcome split with at least two responses (its sampler conditions each response on the
##   others, solve(R[-k, -k])), so it cannot be applied one outcome at a time: it runs in the roles design only.
##   The sensitivity experiment changes only the hyperparameters of the MR.RGM family; the competitors do not use them,
##   and its data are those of the main study (same seeds), so their results for these data are the main-study results.
ALL_METHODS <- c("MR.RGM", "MR.RGM_NoConf", "MR.RGM+", "MR2G", "MVMR", "MR2", "MrDAG", "GrIVET",
                 "MendelianRandomization", "OneSampleMR", "mrbayes")
RGM_FAMILY  <- c("MR.RGM", "MR.RGM_NoConf", "MR.RGM+")
NO_ROLES    <- setdiff(ALL_METHODS, "MR2")                   # designs without an exposure/outcome split
METHODS <- list(
  main        = NO_ROLES,
  noise       = NO_ROLES,
  confounding = NO_ROLES,
  largen      = NO_ROLES,
  longchain   = "MR.RGM+",
  sensitivity = RGM_FAMILY,
  roles       = ALL_METHODS,
  strength    = NO_ROLES
)
## Methods of a setting: the method list of its experiment
methods_for_setting <- function(experiment, setting_label, methods = METHODS[[experiment]]) methods
## Whether a setting label belongs to the current grid of its experiment (a result file of a setting that is no longer
## in the grid is left untouched and not used)
setting_in_grid <- function(experiment, setting_label)
  setting_label %in% vapply(experiment_grid(experiment)$settings, function(s) s$label, character(1))

## Grid of each experiment of run_simulations.R: scenarios, p, n and the list of settings (each a list of SIM_SETTINGS
## fields to change, a label and its position `index` in the grid, stored as column `setting` of the results)
experiment_grid <- function(experiment, set = SIM_SETTINGS) {
  g <- switch(experiment,
    main        = list(scen = c("scalefree", "smallworld", "smallworld_pleio"), p = c(5, 10), n = c(500, 1000, 10000, 30000),
                       settings = list(list(label = "base"))),
    noise       = list(scen = "smallworld", p = 10, n = c(1000, 10000),
                       settings = lapply(c(1, 4, 9), function(v) list(label = sprintf("errvar_%g", v), err_var = v))),
    confounding = list(scen = "smallworld", p = 10, n = c(1000, 10000),
                       settings = lapply(c(0, 0.5, 1, 1.5), function(v) list(label = sprintf("conf_%g", v), conf_scale = v))),
    largen      = list(scen = "smallworld", p = 10, n = c(1e5, 1e6), settings = list(list(label = "base"))),
    ## the setting n = 1e6 of largen (same data and seeds) with a chain five times as long: at this sample size the
    ## posterior is so narrow that the chain of MR.RGM+ needs more iterations to reach the sparse instrument map
    longchain   = list(scen = "smallworld", p = 10, n = 1e6,
                       settings = list(list(label = "iter250000", mcmc = list(nIter = 250000, nBurnin = 50000, Thin = 10)))),
    sensitivity = list(scen = "smallworld", p = 10, n = c(1000, 10000),
                       settings = c(lapply(list(c(3, 1), c(1, 1)), function(ab) list(label = sprintf("Beta_%g_%g", ab[1], ab[2]), rgm_hyper = modifyList(set$rgm_hyper, list(aRho = ab[1], bRho = ab[2])))),
                                    lapply(c(1e-2, 1e-4), function(v) list(label = sprintf("nu1_%g", v), rgm_hyper = modifyList(set$rgm_hyper, list(nu1 = v)))))),
    roles       = list(scen = "roles", p = 8, n = c(1000, 10000), settings = list(list(label = "base"))),
    stop("unknown experiment ", experiment))
  ## the default hyperparameters (Beta(3,1), nu1 = 1e-3) are the "base" row of the sensitivity experiment
  if (experiment == "sensitivity") g$settings <- c(list(list(label = "base")), g$settings[-1])
  for (i in seq_along(g$settings)) g$settings[[i]]$index <- i - 1
  g
}
## Grid of the instrument-strength study (run_instrument_strength.R): the cases and sizes of the main study
STRENGTH_GRID <- list(scen = c("smallworld", "smallworld_pleio"), p = c(5, 10), n = c(500, 1000, 10000, 30000))
## Offsets added to the data seed before each method (independent random streams). Data seeds are below 2e6, so the
## offsets (multiples of 1e7) never give a seed used for data generation or by another method.
METHOD_SEED <- c(MendelianRandomization = 1e7, OneSampleMR = 2e7, mrbayes = 3e7, MVMR = 4e7, MR2G = 5e7, MR2 = 6e7, MrDAG = 7e7,
                 MR.RGM = 8e7, MR.RGM_NoConf = 9e7, "MR.RGM+" = 1e8, GrIVET = 1.1e8)

## Seeds -----------------------------------------------------------------------
truth_seed <- function(rep) as.integer(rep)
## Unique for every (scenario, p, n, replicate) as long as rep < 1000; shared by the settings of an experiment (paired design)
data_seed  <- function(scenario, p, n, rep, set = SIM_SETTINGS) {
  si <- match(scenario, set$scenarios); pi <- match(p, set$p_values); ni <- match(n, set$n_values)
  if (anyNA(c(si, pi, ni))) stop("data_seed: scenario, p or n not in SIM_SETTINGS")
  if (rep < 1 || rep >= 1000) stop("data_seed: rep must be in 1..999")
  as.integer(1e6 + 1e5 * si + 1e4 * pi + 1e3 * ni + rep)
}

## ---------------------------------------------------------------------------
## Truth: graph, A (p x p, A[to, from]), B (p x k), D (instrument map), C (p x l)
## ---------------------------------------------------------------------------
gen_truth <- function(scenario, p, set = SIM_SETTINGS) {

  if (scenario == "scalefree") {
    ## Barabasi-Albert tree; every edge made reciprocal (directed 2-cycles)
    g <- sample_pa(n = p, directed = FALSE)
    g <- as_directed(g, mode = "mutual")
    adj <- as.matrix(as_adjacency_matrix(g, sparse = FALSE)); dimnames(adj) <- NULL   # adj[from, to]
    A <- matrix(sample(set$a_values, p^2, replace = TRUE), p, p) * t(adj)             # A[to, from]
    roles <- NULL
  } else if (scenario %in% c("smallworld", "smallworld_pleio")) {
    ## Watts-Strogatz small-world graph; each edge gets a random direction
    sw <- sample_smallworld(dim = 1, size = p, nei = set$sw_nei, p = set$sw_rewire)
    el <- as_edgelist(sw, names = FALSE)
    el <- t(apply(el, 1, function(e) if (runif(1) > 0.5) e else rev(e)))
    g  <- make_graph(as.vector(t(el)), n = p, directed = TRUE)
    adj <- as.matrix(as_adjacency_matrix(g, sparse = FALSE)); dimnames(adj) <- NULL
    A <- matrix(sample(set$a_values, p^2, replace = TRUE), p, p) * t(adj)
    roles <- NULL
  } else if (scenario == "roles") {
    ## Known roles: p_x exposures (1..p_x) and p_y outcomes (p_x+1..p); edges exposure -> outcome (each outcome
    ## receives one or two exposures) and a random DAG among the outcomes (each later outcome receives at most one
    ## earlier outcome); no edge into an exposure. Effects +-0.3 (a_values of set$roles).
    px <- set$roles$p_x; py <- set$roles$p_y; p <- px + py
    A <- matrix(0, p, p)
    for (j in (px + 1):p) {
      from <- sample(px, sample(1:2, 1))
      A[j, from] <- sample(set$roles$a_values, length(from), replace = TRUE)
      if (j > px + 1 && runif(1) < 0.5) { cand <- (px + 1):(j - 1); k0 <- cand[sample.int(length(cand), 1)]   # sample.int: sample(cand, 1) would misbehave when cand has length 1
                                           A[j, k0] <- sample(set$roles$a_values, 1) }
    }
    roles <- list(X_idx = 1:px, Y_idx = (px + 1):p)
  } else stop("unknown scenario ", scenario)
  diag(A) <- 0

  ## Instruments: ivs_per_trait own instruments per trait, effect iv_effect
  k0 <- set$ivs_per_trait * p
  B  <- matrix(0, p, k0)
  for (i in 1:p) B[i, ((i - 1) * set$ivs_per_trait + 1):(i * set$ivs_per_trait)] <- 1
  D  <- B
  iv_eff <- if (is.null(set$iv_effect)) 1 else set$iv_effect

  if (scenario == "smallworld_pleio") {
    ## One shared (pleiotropic) instrument per consecutive trait pair (1,2), (3,4), ...
    ## (if p is odd, an extra one for (1, p)). The shared instrument affects both traits;
    ## the instrument map D given to the methods keeps only the first trait of each pair.
    starts <- if (p %% 2 == 0) seq(1, p, by = 2) else seq(1, p - 1, by = 2)
    pleio  <- sapply(starts, function(i) { v <- rep(0, p); v[c(i, i + 1)] <- 1; v })
    if (p %% 2 == 1) { v <- rep(0, p); v[c(1, p)] <- 1; pleio <- cbind(pleio, v) }
    pleio  <- matrix(pleio, nrow = p)
    first_only <- pleio
    for (j in seq_len(ncol(pleio))) {
      r1 <- which(pleio[, j] == 1)
      if (length(r1) > 1) first_only[r1[-1], j] <- 0
    }
    B <- cbind(B, pleio)
    D <- cbind(D, first_only)
  }
  B <- B * iv_eff                                        # D (the instrument map) stays 0/1

  ## Latent confounders: l_star factors, loadings +-1 (times conf_scale), half of them set to 0
  C <- matrix(sample(c(-1, 1), p * set$l_star, replace = TRUE), p, set$l_star)
  C[sample(which(C != 0), floor(length(which(C != 0)) / 2))] <- 0
  C <- C * set$conf_scale

  rho <- max(Mod(eigen(A)$values))
  if (rho >= 1) stop("spectral radius of A is ", round(rho, 3), " >= 1")

  list(scenario   = scenario, p = p,
       A          = A,
       graph      = (A != 0) * 1,
       B          = B,
       D          = D,
       C          = C,
       roles      = roles,
       Sigma_star = C %*% t(C) + set$err_var * diag(p))   # covariance of C U + E
}

## ---------------------------------------------------------------------------
## Data: Y = (I - A)^{-1} (B X + C U + E), rows = individuals.
## U and E are drawn before X, so designs that differ only in the instruments (run_instrument_strength.R)
## share the confounders and errors of a replicate.
## ---------------------------------------------------------------------------
gen_data <- function(truth, n, set = SIM_SETTINGS) {
  p <- truth$p; k <- ncol(truth$B); l <- ncol(truth$C)
  U <- matrix(rnorm(n * l), n, l)
  E <- matrix(rnorm(n * p, sd = sqrt(set$err_var)), n, p)
  X <- matrix(rnorm(n * k), n, k)
  M <- solve(diag(p) - truth$A)
  Y <- (X %*% t(truth$B) + U %*% t(truth$C) + E) %*% t(M)
  list(X = X, Y = Y)
}

## Center and scale X and Y (as recommended for the SSSL prior and done for the real data) and express the
## truth on the standardized scale: a*_ij = a_ij * s_j / s_i, b*_ih = b_ih * t_h / s_i, Sigma* -> S^-1 Sigma* S^-1,
## where s = sd(Y), t = sd(X). The graph and the instrument map are unchanged. The confounding truth used by conf_metrics
## (pairs whose |Sigma*_jh| exceeds the mean over pairs) is kept on the generating scale (Sigma_star_gen), because the
## rescaling changes the ranking of the entries; the fitted indicators Z are scale free.
standardize_data <- function(dat, truth, set = SIM_SETTINGS) {
  if (!isTRUE(set$standardize)) return(list(X = dat$X, Y = dat$Y, truth = truth))
  sY <- apply(dat$Y, 2, sd); sX <- apply(dat$X, 2, sd)
  Y <- scale(dat$Y); X <- scale(dat$X)
  attr(Y, "scaled:center") <- attr(Y, "scaled:scale") <- NULL
  attr(X, "scaled:center") <- attr(X, "scaled:scale") <- NULL
  t2 <- truth
  t2$A <- truth$A * outer(1 / sY, sY)          # A[to, from] * s_from / s_to
  t2$B <- truth$B * outer(1 / sY, sX)
  t2$Sigma_star <- truth$Sigma_star * outer(1 / sY, 1 / sY)
  t2$Sigma_star_gen <- truth$Sigma_star
  t2$scale <- list(sY = sY, sX = sX)
  list(X = X, Y = Y, truth = t2)
}

## ---------------------------------------------------------------------------
## Metrics
## ---------------------------------------------------------------------------
## TPR / FDR / MCC (TPR = 0 if there are no true positives, FDR = 0 if there are no discoveries)
metrics_binary_graph <- function(true_graph, estimate_graph) {
  stopifnot(length(true_graph) == length(estimate_graph))
  estimate_graph[is.na(estimate_graph)] <- 0
  TP <- sum(true_graph == 1 & estimate_graph == 1)
  FP <- sum(true_graph == 0 & estimate_graph == 1)
  FN <- sum(true_graph == 1 & estimate_graph == 0)
  TN <- sum(true_graph == 0 & estimate_graph == 0)
  TPR <- TP / max(1, sum(true_graph))
  FDR <- FP / max(1, (TP + FP))
  denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  MCC <- if (denom == 0) 0 else (TP * TN - FP * FN) / denom
  c(TPR = TPR, FDR = FDR, MCC = MCC)
}

## AUC with a fixed direction: larger score = more evidence for an edge
auc_safe <- function(truth, score) {
  ok <- !is.na(score)
  if (length(unique(truth[ok])) < 2 || sum(ok) < 2) return(NA_real_)
  tryCatch(as.numeric(pROC::auc(pROC::roc(truth[ok], score[ok], levels = c(0, 1),
                                          direction = "<", quiet = TRUE))),
           error = function(e) NA_real_)
}

off_diag <- function(M) M[row(M) != col(M)]

## Graph recovery; all matrices oriented [to, from] like A. `sel` restricts the entries (e.g. exposure -> outcome).
graph_metrics <- function(truth_graph, decision, score, sel = NULL) {
  if (is.null(sel)) sel <- row(truth_graph) != col(truth_graph)
  tg <- truth_graph[sel]
  c(AUC = auc_safe(tg, score[sel]), metrics_binary_graph(tg, decision[sel] * 1))
}

## Causal effect estimation; A_hat oriented [to, from]; NA entries (not estimated) are ignored
effect_metrics <- function(A, A_hat, sel = NULL) {
  if (is.null(sel)) sel <- row(A) != col(A)
  d <- abs(A - A_hat)[sel]
  if (all(is.na(d))) return(c(MaxAbsDev = NA_real_, MeanAbsDev = NA_real_, MeanSqDev = NA_real_))
  c(MaxAbsDev = max(d, na.rm = TRUE), MeanAbsDev = mean(d, na.rm = TRUE), MeanSqDev = mean(d^2, na.rm = TRUE))
}

## Effect error by the truth, computed from the saved estimates (fits_*.rds; used by make_figures_tables.R and
## run_instrument_strength.R): the mean absolute deviation over the true edges (A != 0) and over the non-edges (A = 0:
## the mean absolute effect estimated where there is none), and the reference "All-zero": the effect metrics of the
## estimate A_hat = 0, i.e. of a method that reports no effect at all (its MeanAbsDev is mean |A| over the
## off-diagonal entries). The MAD over all entries is dominated by the many zero entries of a sparse network, so a
## method that shrinks every effect to zero attains the All-zero value without finding any edge.
## Returns one row per (replicate, method, family, metric): families effect_edges, effect_nonedges (every method) and
## effect (method "All-zero").
effect_split_rows <- function(fits_file) {
  fo <- readRDS(fits_file)
  avg <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  do.call(rbind, lapply(fo$fits, function(fi) {
    A <- fi$truth$A; off <- row(A) != col(A); e <- off & A != 0; z <- off & A == 0
    rows <- lapply(names(fi$methods), function(m) {
      x <- fi$methods[[m]]; Ah <- if (!is.null(x$AEst)) x$AEst else x$est
      if (is.null(Ah)) return(NULL)
      data.frame(method = m, family = c("effect_edges", "effect_nonedges"), metric = "MeanAbsDev",
                 value = c(avg(abs(A - Ah)[e]), avg(abs(Ah)[z])), stringsAsFactors = FALSE) })
    z0 <- effect_metrics(A, 0 * A)
    rows[[length(rows) + 1]] <- data.frame(method = "All-zero", family = "effect", metric = names(z0), value = unname(z0),
                                           stringsAsFactors = FALSE)
    out <- do.call(rbind, rows); out$rep <- fi$rep; out
  }))
}

## Confounding structure: truth = |Sigma*| (upper triangle) / max, thresholded at its mean
conf_metrics <- function(Sigma_star, Z_hat, thr = SIM_SETTINGS$pip_thr) {
  ut <- upper.tri(Sigma_star)
  s  <- abs(Sigma_star[ut])
  if (max(s) == 0) return(c(AUC = NA_real_, TPR = NA_real_, FDR = NA_real_, MCC = NA_real_))   # no confounding
  s  <- s / max(s)
  tg <- (s > mean(s)) * 1
  if (length(unique(tg)) < 2) return(c(AUC = NA_real_, TPR = NA_real_, FDR = NA_real_, MCC = NA_real_))
  z  <- Z_hat[ut]
  c(AUC = auc_safe(tg, z), metrics_binary_graph(tg, (z > thr) * 1))
}

## Instrument-trait selection (MR.RGM+): truth = B != 0, score = PhiEst (posterior inclusion probability of b_jh)
iv_metrics <- function(B, Phi_hat) c(AUC = auc_safe(as.vector((B != 0) * 1), as.vector(Phi_hat)))

## Weak-identification diagnostics of an MR.RGM fit
rgm_diagnostics <- function(fit, Syy) {
  v <- if (is.matrix(fit$SigmaEst)) diag(fit$SigmaEst) else as.numeric(fit$SigmaEst)
  ratio <- v / diag(Syy)
  GP <- fit$GammaPst; m <- dim(GP)[3]; h <- m %/% 2; off <- row(fit$GammaEst) != col(fit$GammaEst)
  halves <- suppressWarnings(cor(apply(GP[, , 1:h, drop = FALSE], c(1, 2), mean)[off],
                                 apply(GP[, , (h + 1):m, drop = FALSE], c(1, 2), mean)[off]))
  hd <- max(abs(apply(GP[, , 1:h, drop = FALSE], c(1, 2), mean)[off] - apply(GP[, , (h + 1):m, drop = FALSE], c(1, 2), mean)[off]))
  w <- if (is.null(fit$warnings)) character(0) else fit$warnings
  c(failed_fits = 0, err_ratio_median = median(ratio), err_ratio_max = max(ratio), warning_fired = as.numeric(max(ratio) > 5),
    warn_errvar = as.numeric(any(grepl("error variances exceed", w))), warn_weakIV = as.numeric(any(grepl("Weak instruments", w))),
    warn_spectral = as.numeric(any(grepl("spectral radius", w))), warn_halves = as.numeric(any(grepl("two halves", w))),
    warn_accept = as.numeric(any(grepl("acceptance rate", w))), warn_any = as.numeric(length(w) > 0),
    pip_halves = halves, halves_maxdiff = hd, spectral_radius = max(Mod(eigen(fit$AEst)$values)),
    accept_A = if (is.null(fit$AccptA)) NA_real_ else fit$AccptA,
    min_F = if (is.null(fit$diagnostics$instrumentF)) NA_real_ else suppressWarnings(min(fit$diagnostics$instrumentF, na.rm = TRUE)))
}

## Population instrument strength of the simulated traits: R^2 of the regression of Y_j on the instruments that the
## map D assigns to it (reduced form Pi = (I - A)^{-1} B, Var(Y) = M (B B' + Sigma*) M', X independent with unit
## variance). Scale invariant, so the generating-scale truth can be used; returns median, min and max over traits.
## Also returned: the implied mu^2 = n R^2 / (1 - R^2) and the F statistic for k = number of own instruments.
pop_r2 <- function(truth, n = NULL) {
  M  <- solve(diag(truth$p) - truth$A); Pi <- M %*% truth$B
  vy <- diag(M %*% (tcrossprod(truth$B) + truth$Sigma_star) %*% t(M))
  r2 <- vapply(seq_len(truth$p), function(j) { own <- truth$D[j, ] == 1; sum(Pi[j, own]^2) / vy[j] }, numeric(1))
  k  <- rowSums(truth$D)
  out <- c(R2_median = median(r2), R2_min = min(r2), R2_max = max(r2))
  if (!is.null(n)) out <- c(out, mu2_median = median(n * r2 / (1 - r2)),
                            F_median = median((r2 / k) / ((1 - r2) / (n - k - 1))))
  out
}

## ---------------------------------------------------------------------------
## Per-SNP marginal summary statistics (GWAS-type input of the summary-level competitors), computed once
## per replicate from the same sample: beta[h, j] and se[h, j] of the regression of Y_j on X_h.
## ---------------------------------------------------------------------------
marginal_stats <- function(X, Y) {
  n <- nrow(X); k <- ncol(X); p <- ncol(Y)
  Xc <- scale(X, center = TRUE, scale = FALSE); Yc <- scale(Y, center = TRUE, scale = FALSE)
  sxx <- colSums(Xc^2)
  beta <- crossprod(Xc, Yc) / sxx                                   # k x p
  resid_ss <- matrix(colSums(Yc^2), k, p, byrow = TRUE) - beta^2 * sxx
  se <- sqrt(resid_ss / (n - 2) / sxx)
  snp <- paste0("snp", seq_len(k)); trait <- paste0("Y", seq_len(p))
  dimnames(beta) <- dimnames(se) <- list(snp, trait)
  list(beta = beta, se = se, n = n, snp = snp, trait = trait)
}

## ---------------------------------------------------------------------------
## Competing methods. Pairwise: exposure i -> outcome j using the instruments of
## trait i (D[i, ]). Returned matrices are [exposure, outcome]; transpose to
## compare with A[to, from].
## ---------------------------------------------------------------------------
snp_betas <- function(y, G) {
  t(sapply(seq_len(ncol(G)), function(m) {
    co <- summary(lm(y ~ G[, m]))$coefficients
    c(beta = co[2, 1], se = co[2, 2])
  }))
}

## A pair whose fit errors is recorded as NA and counted in `failed`; NA scores are left out of the AUC (auc_safe) and
## an NA p-value is not an edge (the number of such pairs is stored as failed_fits)
fit_mr_package <- function(X, Y, D, set = SIM_SETTINGS) {
  p <- ncol(Y)
  meth <- c(SimpleMedian = "Simple median", WeightedMedian = "Weighted median", IVW = "IVW")
  est <- lapply(meth, function(m) matrix(NA_real_, p, p))
  pv  <- lapply(meth, function(m) matrix(NA_real_, p, p))
  failed <- 0
  for (i in 1:p) {
    ivs <- which(D[i, ] == 1)
    if (length(ivs) < 3) next
    G  <- X[, ivs, drop = FALSE]
    bx <- snp_betas(Y[, i], G)
    for (j in 1:p) {
      if (i == j) next
      by  <- snp_betas(Y[, j], G)
      obj <- mr_input(bx = bx[, 1], bxse = bx[, 2], by = by[, 1], byse = by[, 2],
                      snps = paste0("rs", ivs))
      res <- tryCatch(mr_allmethods(obj, method = "main", iterations = set$mr_boot),
                      error = function(e) NULL)
      if (is.null(res)) { failed <- failed + 1; next }
      v <- res@Values
      for (m in names(meth)) {
        est[[m]][i, j] <- v[v$Method == meth[[m]], "Estimate"]
        pv[[m]][i, j]  <- v[v$Method == meth[[m]], "P-value"]
      }
    }
  }
  list(est = est, pval = pv, failed = failed)
}

fit_onesamplemr <- function(X, Y, D) {
  p <- ncol(Y)
  est <- matrix(NA_real_, p, p); pv <- matrix(NA_real_, p, p); failed <- 0
  for (i in 1:p) {
    ivs <- which(D[i, ] == 1)
    if (length(ivs) < 1) next
    G <- X[, ivs, drop = FALSE]
    nm <- paste0("z", ivs); colnames(G) <- nm
    for (j in 1:p) {
      if (i == j) next
      ## OneSampleMR >= 0.1.8 fits TSRI by GMM (gmm::gmm). On centered data both intercepts are exactly 0, the numerical
      ## gradient of gmm then has a zero column and the coefficient covariance is singular ("The covariance matrix of
      ## the coefficients is singular": every SE Inf, every p-value 1). The identity-link TSRI is invariant to shifts of
      ## the outcome and the exposure (the causal slope and its SE are unchanged), so constants are added to give
      ## nonzero intercepts.
      df  <- data.frame(y = Y[, j] + 3, x = Y[, i] + 1, G)
      fit <- tryCatch(suppressWarnings(tsri(y ~ x, reformulate(nm), df, link = "identity")), error = function(e) NULL)
      ok <- FALSE
      if (!is.null(fit) && !is.null(fit$fit)) {
        co <- summary(fit$fit)$coefficients
        if ("x" %in% rownames(co) && is.finite(co["x", "Std. Error"])) {   # a singular covariance is a failed fit
          est[i, j] <- co["x", "Estimate"]
          pv[i, j]  <- co["x", "Pr(>|t|)"]
          ok <- TRUE
        }
      }
      if (!ok) failed <- failed + 1
    }
  }
  list(est = est, pval = pv, failed = failed)
}

fit_mrbayes <- function(X, Y, D, set = SIM_SETTINGS) {
  p <- ncol(Y)
  est <- matrix(NA_real_, p, p); pv <- matrix(NA_real_, p, p); failed <- 0
  has_seed <- "seed" %in% names(formals(mrbayes::mr_ivw_rjags))
  for (i in 1:p) {
    ivs <- which(D[i, ] == 1)
    if (length(ivs) < 2) next
    G  <- X[, ivs, drop = FALSE]
    bx <- snp_betas(Y[, i], G)
    for (j in 1:p) {
      if (i == j) next
      by <- snp_betas(Y[, j], G)
      dat <- mr_format(rsid = paste0("SNP", ivs), xbeta = bx[, 1], ybeta = by[, 1],
                       xse = bx[, 2], yse = by[, 2])
      args <- list(object = dat, prior = "default", n.chains = set$mrbayes$n.chains,
                   n.burn = set$mrbayes$n.burn, n.iter = set$mrbayes$n.iter)
      if (has_seed) args$seed <- sample.int(.Machine$integer.max, set$mrbayes$n.chains)   # distinct per pair and chain
      fit <- tryCatch(do.call(mr_ivw_rjags, args), error = function(e) NULL)
      if (is.null(fit)) { failed <- failed + 1; next }
      s <- do.call(rbind.data.frame, fit$samples)$Estimate
      est[i, j] <- mean(fit$CausalEffect)
      pv[i, j]  <- 2 * min(mean(s > 0), mean(s < 0))
    }
  }
  list(est = est, pval = pv, failed = failed)
}

## ---------------------------------------------------------------------------
## MVMR (MendelianRandomization::mr_mvivw, package defaults): one outcome at a time, all other traits as
## exposures, instruments = the instruments of the exposures (rows of D). Returned matrices are [to, from]
## (direct effect of each exposure on the outcome, conditional on the other exposures).
## ---------------------------------------------------------------------------
fit_mvmr <- function(ms, D) {
  p <- ncol(ms$beta)
  est <- matrix(NA_real_, p, p); pv <- matrix(NA_real_, p, p); failed <- 0
  for (j in 1:p) {
    ex  <- setdiff(1:p, j)
    ivs <- which(colSums(D[ex, , drop = FALSE]) > 0)
    if (length(ivs) <= length(ex)) { failed <- failed + 1; next }     # MVMR needs more instruments than exposures
    obj <- mr_mvinput(bx = ms$beta[ivs, ex, drop = FALSE], bxse = ms$se[ivs, ex, drop = FALSE],
                      by = ms$beta[ivs, j], byse = ms$se[ivs, j], snps = ms$snp[ivs])
    res <- tryCatch(mr_mvivw(obj), error = function(e) NULL)
    if (is.null(res)) { failed <- failed + 1; next }
    est[j, ex] <- res@Estimate
    pv[j, ex]  <- res@Pvalue
  }
  list(est = est, pval = pv, failed = failed)
}

## ---------------------------------------------------------------------------
## MR2G (Lin, Pan & Xue 2026; GitHub ZhaotongL/MR2G), run with the defaults of MR2G(): IV_screen = "ScreenAug",
## num_pert = num_pert_M = 100. Inputs built from the same sample: per-SNP marginal effects and standard errors
## (b_mat, se_mat: SNPs x traits), per-trait sample sizes, the trait correlation matrix rho_mat (full sample overlap:
## the phenotypic correlation of the traits), the candidate instruments of each trait pair (IV_list: the union of the
## instruments of the two traits according to the map D, in the order (1,2), (1,3), ..., (p-1,p)) and the LD matrix of
## all SNPs (R_list, one block). MR2G screens the candidates itself (genome-wide significance, p < 5e-8, hard-coded in
## the package), so with weak instruments it may have fewer than 3 instruments for a pair; the package warns and the
## result is kept (NA if the fit fails). Output: Direct_graph_Est / _pval, rows = outcomes, columns = exposures
## (as stated in the package README), i.e. the [to, from] orientation of A.
##
## Compatibility with igraph >= 2.0 (MR2G 0.1.0): in the ScreenAug step MR2G() builds a graph from
## `Direct_graph_pval < pthres`, whose diagonal is NaN (the direct effect of a trait on itself is 0 in every
## perturbation, so its z-statistic is 0/0); igraph >= 2.0 refuses an adjacency matrix with NA ("Cannot create a graph
## object because the adjacency matrix contains NAs"). mr2g_screenaug() below is MR2G() with IV_screen = "ScreenAug"
## line by line (same internal functions, same defaults: pthres = 0.05 / (p^2 - p), random_start = 0, check = TRUE,
## num_pert_M perturbations in both phases), except that NA entries of that adjacency matrix are set to FALSE (no
## self-loops), which is what the step intends. Nothing else is changed.
## ---------------------------------------------------------------------------
mr2g_screenaug <- function(b_mat, se_mat, n_vec, rho_mat, IV_list, R_list, trait_vec, num_pert_M = 100, seed = 0,
                           random_start = 0, check = TRUE) {
  ns <- asNamespace("MR2G")
  pthres <- 0.05 / ((length(n_vec))^2 - length(n_vec))                       # MR2G() default
  SP_res <- ns$ScreenPair(b_mat = b_mat, se_mat = se_mat, pval_mat = NULL, n_mat = NULL, n_vec = n_vec,
                          IV_list = IV_list, R_list = R_list, rho_mat = rho_mat)
  SM_res <- ns$ScreenMax(b_mat = b_mat, se_mat = se_mat, n_vec = n_vec, pre_screen = SP_res$IJ_snp_list)
  SM_res$DP_mat_list <- SP_res$DP_mat_list
  dp_list_SM <- ns$Graph_Perturb_fast(screen_res = SM_res, n_vec = n_vec, rho_mat = rho_mat, random_start = random_start,
                                      num_pert = num_pert_M, seed = seed, trait_vec = trait_vec)
  resMax <- ns$subset_Graph(dp_list = dp_list_SM, keep_trait = "all", check = check)
  B <- resMax$dp_res$Direct_graph_pval < pthres
  B[is.na(B)] <- FALSE                                                       # the only change (see above)
  gg <- igraph::graph_from_adjacency_matrix(B * 1, mode = "directed")
  SA_res <- ns$ScreenAug(b_mat = b_mat, se_mat = se_mat, n_vec = n_vec, pre_screen_pair = SP_res$IJ_snp_list,
                         pre_screen_max = SM_res$IJ_snp_list, graph = gg)
  SA_res$DP_mat_list <- SP_res$DP_mat_list
  dp_list_SA <- ns$Graph_Perturb_fast(screen_res = SA_res, n_vec = n_vec, rho_mat = rho_mat, random_start = random_start,
                                      num_pert = num_pert_M, seed = seed, trait_vec = trait_vec)
  resAug <- ns$subset_Graph(dp_list = dp_list_SA, keep_trait = "all", check = check)
  list(ScreenAug = lapply(resAug$dp_res, t), ScreenMax = lapply(resMax$dp_res, t))
}

fit_mr2g <- function(ms, Y, D, set = SIM_SETTINGS) {
  if (!has_pkg("MR2G")) stop("package MR2G not installed")
  p <- ncol(ms$beta)
  b_mat <- ms$beta; se_mat <- ms$se
  colnames(b_mat) <- colnames(se_mat) <- ms$trait
  n_vec <- rep(ms$n, p)
  rho_mat <- cor(Y); dimnames(rho_mat) <- list(ms$trait, ms$trait)
  IV_list <- list(); k <- 1
  for (i in 1:(p - 1)) for (j in (i + 1):p) { IV_list[[k]] <- ms$snp[D[i, ] == 1 | D[j, ] == 1]; k <- k + 1 }
  R <- cor(as.matrix(attr(ms, "X")))
  dimnames(R) <- list(ms$snp, ms$snp)
  R_list <- list(list(R = R, snp = ms$snp))
  if (has_pkg("future")) future::plan(future::sequential)        # no nested parallelism inside the replicate workers
  seed <- sample.int(.Machine$integer.max, 1)                     # MR2G calls set.seed(seed) itself (default 0)
  out <- if (set$mr2g$IV_screen == "ScreenAug") {
    mr2g_screenaug(b_mat, se_mat, n_vec, rho_mat, IV_list, R_list, trait_vec = ms$trait,
                   num_pert_M = set$mr2g$num_pert_M, seed = seed)$ScreenAug
  } else {
    MR2G::MR2G(b_mat = b_mat, se_mat = se_mat, n_vec = n_vec, rho_mat = rho_mat, IV_list = IV_list, R_list = R_list,
               trait_vec = ms$trait, IV_screen = set$mr2g$IV_screen, num_pert = set$mr2g$num_pert,
               num_pert_M = set$mr2g$num_pert_M, seed = seed)
  }
  est <- as.matrix(out$Direct_graph_Est); pv <- as.matrix(out$Direct_graph_pval)
  if (!is.null(rownames(est)) && all(ms$trait %in% rownames(est))) { est <- est[ms$trait, ms$trait]; pv <- pv[ms$trait, ms$trait] }
  dimnames(est) <- dimnames(pv) <- NULL
  diag(est) <- NA; diag(pv) <- NA
  list(est = est, pval = pv, failed = 0)                            # [to, from]
}

## ---------------------------------------------------------------------------
## MR2 (Zuber et al. 2023; package MR2): multi-response MR with exposures X_idx and responses Y_idx (roles design;
## designs without roles: by_outcome() below).
## Input: per-SNP marginal effects of the instruments of the exposures (rows of D of the exposures) on the exposures
## (betaHat_X) and on the responses (betaHat_Y). Chain length as for every MCMC method of the study (SIM_SETTINGS$mr2:
## 50,000 iterations, 10,000 burn-in, thinning 10; the package has no default), EVgamma = c(2, 2) (package default);
## everything else default. Output (PostProc): gammaPost = marginal
## posterior inclusion probability and thetaPost = posterior mean of the direct effect of each exposure on each
## response (responses x exposures), RPost = posterior mean residual correlation between the responses.
##
## Scale: the default priors of MR2 assume summary statistics of unit order (hyper-inverse-Wishart prior on the residual
## covariance with identity scale; the package's simulators draw effects of order 1). The marginal effects on the
## standardized data are of order 1 / sqrt(n) for their noise, so with the defaults the prior swamps the data (pilot:
## every inclusion probability ~ 0.6 and every effect ~ 0). Both matrices are therefore multiplied by sqrt(n), which
## puts the noise on unit scale (approximately z-statistics, since every SE is about 1 / sqrt(n) on standardized data)
## and leaves the direct effects theta unchanged (the model is linear in the common scale).
## ---------------------------------------------------------------------------
fit_mr2 <- function(ms, D, roles, set = SIM_SETTINGS) {
  if (!has_pkg("MR2")) stop("package MR2 not installed")
  p <- ncol(ms$beta); Xi <- roles$X_idx; Yi <- roles$Y_idx
  if (length(Xi) < 4) stop("MR2 with its default EVgamma = c(2, 2) needs at least 4 exposures (", length(Xi), " given)")
  if (length(Yi) < 2) stop("MR2 needs at least 2 responses (", length(Yi), " given)")
  ivs <- which(colSums(D[Xi, , drop = FALSE]) > 0)
  sc <- sqrt(ms$n)
  bY <- sc * ms$beta[ivs, Yi, drop = FALSE]; bX <- sc * ms$beta[ivs, Xi, drop = FALSE]
  out <- MR2::MR2(bY, bX, EVgamma = set$mr2$EVgamma, niter = set$mr2$niter, burnin = set$mr2$burnin,
                  thin = set$mr2$thin, monitor = set$mr2$niter + 1, seed = sample.int(.Machine$integer.max, 1))
  pp  <- MR2::PostProc(out, bY, bX)
  g <- as.matrix(pp$gammaPost); th <- as.matrix(pp$thetaPost)
  ## PostProc() returns gammaPost and thetaPost as (exposures x responses) (package source; confirmed by pilot_competitors.R), so transpose
  g <- t(g); th <- t(th)                                             # (responses x exposures) = [to, from]
  est <- matrix(NA_real_, p, p); score <- matrix(NA_real_, p, p)
  est[Yi, Xi] <- th; score[Yi, Xi] <- g
  Rres <- NULL                                                       # residual correlation: defined with >= 2 responses
  if (length(Yi) > 1 && !is.null(pp$RPost)) {
    R <- as.matrix(pp$RPost); R[upper.tri(R)] <- t(R)[upper.tri(R)]; diag(R) <- 1   # PostProc() fills only the lower triangle
    Rres <- matrix(NA_real_, p, p); Rres[Yi, Yi] <- R
  }
  list(est = est, score = score, resid_cor = Rres, failed = 0)        # [to, from]; score = PIP
}

## ---------------------------------------------------------------------------
## MrDAG (Zuber et al. 2025; package MrDAG): Bayesian DAG learning with exposures X_idx and outcomes Y_idx
## (roles design; designs without roles: by_outcome() below); input = per-SNP marginal effects of the exposures' instruments on all traits, outcomes first.
## Chain length as for every MCMC method (SIM_SETTINGS$mrdag: 50,000 / 10,000 / 10); pp = 0.05 and tempMax = 10
## (package defaults). Output: get_edgeprob() = posterior edge inclusion probabilities (PEIP) among the outcomes and
## among the exposures; get_causaleffects() = posterior genetically predicted effects of the exposures on the outcomes.
## Scale: as for MR2, the marginal effects are multiplied by sqrt(n) (default Wishart prior U = I / nIV assumes unit-order
## data); the causal effects (ratios of covariances) and the graph are unchanged by a common scale.
## ---------------------------------------------------------------------------
fit_mrdag <- function(ms, D, roles, set = SIM_SETTINGS) {
  if (!has_pkg("MrDAG")) stop("package MrDAG not installed")
  p <- ncol(ms$beta); Xi <- roles$X_idx; Yi <- roles$Y_idx
  ivs <- which(colSums(D[Xi, , drop = FALSE]) > 0)
  dat <- sqrt(ms$n) * ms$beta[ivs, c(Yi, Xi), drop = FALSE]
  chk <- list(Y_idx = seq_along(Yi), X_idx = length(Yi) + seq_along(Xi))
  out <- MrDAG::MrDAG(data = dat, niter = set$mrdag$niter, burnin = set$mrdag$burnin, thin = set$mrdag$thin,
                      pp = set$mrdag$pp, tempMax = set$mrdag$tempMax, MrDAGcheck = chk, fileName = NULL,
                      seed = sample.int(.Machine$integer.max, 1))
  peip <- as.matrix(MrDAG::get_edgeprob(out)$edgeProb)              # traits in the order (outcomes, exposures)
  ce   <- MrDAG::get_causaleffects(out, progress = FALSE)$causalEffects
  ce   <- as.matrix(ce)
  ord  <- c(Yi, Xi)
  ## edgeProb[i, j] is the PEIP of the edge i -> j (adjacency convention, row = from; the package sets the forbidden
  ## [outcome, exposure] block to NA), so score[to, from] uses the transpose. causalEffects is (exposures x outcomes)
  ## (package source; confirmed by pilot_competitors.R), so it is transposed too.
  score <- matrix(NA_real_, p, p)
  score[ord, ord] <- t(peip)
  est <- matrix(NA_real_, p, p)
  est[Yi, Xi] <- t(ce)
  list(est = est, score = score, failed = 0)                          # [to, from]; score = PEIP (NA between exposures and outcomes)
}

## ---------------------------------------------------------------------------
## MrDAG in the designs without known roles (every trait may affect every other, feedback loops): it needs a split of
## the traits into exposures and outcomes, with effects only from exposures to outcomes. It is applied one outcome at a
## time, exactly as MVMR (MR2 cannot be: it needs at least two responses): for each trait j, j is the single outcome and the other
## p - 1 traits are the exposures, with the instruments of the exposures (rows of D); row j of the [to, from] matrices
## (effects and inclusion probabilities of the edges into j) is taken from that fit. Every edge i -> j of the network is
## thus estimated, conditional on the other traits, by the fit in which j is the outcome. An outcome whose fit errors is
## left as NA and counted in `failed`; if every fit errors the method has failed (the first error is passed on).
## MrDAG's edges among the exposures of one fit are not used (each edge is taken from the fit in which its target is
## the outcome).
## MrDAG's effects are its causal effects of the exposures on the outcome (get_causaleffects, effects under intervention
## within its DAG, so they include paths through other exposures), as in the roles design; its edge probabilities are
## those of the direct edges. The number of failed outcomes is stored as failed_fits (as for MVMR and the pairwise methods).
## ---------------------------------------------------------------------------
by_outcome <- function(fit_fun, ms, D, set = SIM_SETTINGS) {
  p <- ncol(ms$beta)
  est <- matrix(NA_real_, p, p); score <- matrix(NA_real_, p, p); failed <- 0; first_err <- NULL
  for (j in 1:p) {
    r <- tryCatch(fit_fun(ms, D, list(X_idx = setdiff(1:p, j), Y_idx = j), set),
                  error = function(e) {
                    if (is.null(first_err)) first_err <<- conditionMessage(e)
                    lf <- getOption("sim_error_log")                 # one line per failed outcome
                    if (!is.null(lf)) cat(sprintf("[%s] by_outcome: outcome %d of %d failed: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                                                  j, p, conditionMessage(e)), file = lf, append = TRUE)
                    NULL })
    if (is.null(r)) { failed <- failed + 1; next }
    est[j, -j] <- r$est[j, -j]; score[j, -j] <- r$score[j, -j]
  }
  if (failed == p) stop(first_err)
  list(est = est, score = score, failed = failed)                     # [to, from]
}

## ---------------------------------------------------------------------------
## MR.RGM family (summary-level input of the standardized data)
##   MR.RGM        : D = own-instrument map, selectB = FALSE (no selection of the allowed instrument effects)
##   MR.RGM_NoConf : as MR.RGM with SigmaStarModel = "diagonal"
##   MR.RGM+       : D = 1, selectB = TRUE (every instrument effect selected)
## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------
## GrIVET (Chen, Li, Shen and Pan 2024, JASA; package grivet, GitHub chunlinli/grivet): causal discovery for a directed
## ACYCLIC graph with hidden confounding, from instruments, on individual-level data. Its data-driven pipeline (that of
## its own discovery simulations), with no information on the truth and no instrument map (every instrument is a
## candidate for every trait):
##   1. cv.intdag.pmle.diff.aic()  instrument -> trait effects V (truncated-L1 regressions, cross-validated, AIC refit)
##   2. topological_order(V)       peeling: repeatedly removes a trait whose instrument acts on no other remaining
##                                 trait, giving the ancestral relations; on a directed cycle there is no such trait,
##                                 and a trait that the peeling does not reach makes the fit stop (a failed fit)
##   3. cv.intdag.coe.penalized()  causal effects given the ancestral relations (U[from, to])
## Input: the UNSTANDARDIZED data. Its penalties are not scale free: on the standardized data it selects (almost) no
## instrument and no edge, even with the extended grid (pilot runs), whereas on the unstandardized data it
## recovers the exposure -> outcome design; every other method is fitted to the standardized data. The estimated effects
## are converted to the standardized scale exactly as the truth is (standardize_data(): a[to, from] * s_from / s_to),
## so the effect metrics are comparable. Edge score = |effect| (GrIVET gives no p-values or probabilities), edge = effect
## != 0 (its own selection).
## ---------------------------------------------------------------------------
## Compatibility correction for grivet under R >= 4.3. In the internal functions cv.tlpreg0.aic and cv.tlpreg0, the
## tuning parameters are chosen by
##   try(locations <- <index of the smallest cross-validation error>, { warning(...); locations <- <middle of the grid> })
## where the second argument, intended as the fallback when every cross-validation error is missing, is the argument
## `silent` of try(). Before R 4.3, try() evaluated it on an error (so the fallback ran) and used its first element as
## `silent`; R >= 4.3 stops with "'length = 2' in coercion to 'logical(1)'". The correction replaces this statement by
##   locations <- tryCatch(<index of the smallest cross-validation error>, error = function(e) { <the same fallback> })
## which gives the result of the package under R < 4.3 (the fallback block is the package's own code). The functions
## are patched in the loaded namespace only; the installed package is not modified. Idempotent.
grivet_compat <- function() {
  ns <- asNamespace("grivet")
  for (fn in c("cv.tlpreg0.aic", "cv.tlpreg0")) {
    if (!exists(fn, envir = ns, inherits = FALSE)) next
    f <- get(fn, envir = ns); b <- body(f)
    k <- which(vapply(as.list(b), function(s) is.call(s) && identical(s[[1]], as.name("try")), logical(1)))
    if (length(k) != 1) next                                       # already corrected, or a version without the construct
    tr <- b[[k]]
    if (length(tr) != 3 || !is.call(tr[[2]]) || !identical(tr[[2]][[1]], as.name("<-")) ||
        !identical(tr[[2]][[2]], as.name("locations"))) next
    b[[k]] <- bquote(locations <- tryCatch(.(tr[[2]][[3]]), error = function(e) .(tr[[3]])))
    body(f) <- b
    utils::assignInNamespace(fn, f, ns = "grivet")
  }
  invisible(TRUE)
}

fit_grivet <- function(X_raw, Y_raw, sY, set = SIM_SETTINGS) {
  if (!has_pkg("grivet")) stop("package grivet not installed")
  grivet_compat()
  g  <- set$grivet
  s1 <- grivet::cv.intdag.pmle.diff.aic(X_raw, Y_raw, g$tau, g$gamma, g$n_fold)
  s2 <- suppressWarnings(grivet::topological_order(s1$V))         # warnings: min() of an empty set when the peeling ends
  reached <- colSums(s2$in_mat) > 0
  if (!all(reached)) stop(sprintf("GrIVET peeling reached %d of %d traits", sum(reached), ncol(Y_raw)))
  s3 <- grivet::cv.intdag.coe.penalized(X_raw, Y_raw, s2$an_mat, s2$in_mat, s2$iv_mat, g$tau, g$gamma, g$n_fold)
  est <- t(s3$U) * outer(1 / sY, sY)                               # [to, from], standardized scale
  diag(est) <- NA
  list(est = est, score = abs(est), edge = est != 0, failed = 0)
}

fit_rgm <- function(X, Y, D, model = c("SSSL", "diagonal"), set = SIM_SETTINGS, selectB = NULL, hyper = set$rgm_hyper) {
  model <- match.arg(model)
  n <- nrow(Y)
  if (is.null(selectB)) selectB <- all(D == 1)
  ## the package's warnings are recorded (fit$warnings; counted in rgm_diagnostics) instead of printed by the workers
  ws <- character(0)
  fit <- withCallingHandlers(
    RGM(Syy = crossprod(Y) / n, Syx = crossprod(Y, X) / n, Sxx = crossprod(X) / n,
        D = D, n = n, selectB = selectB,
        nIter = set$mcmc$nIter, nBurnin = set$mcmc$nBurnin, Thin = set$mcmc$Thin,
        prior = "Spike and Slab", SigmaStarModel = model,
        aRho = hyper$aRho, bRho = hyper$bRho, nu1 = hyper$nu1),
    warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
  fit$warnings <- ws
  fit
}

## ---------------------------------------------------------------------------
## One replicate: the methods of `experiment` on the same data; returns a long data.frame
##   setting: integer index of the experiment-specific setting (noise level, confounding scale, hyperparameter set)
##   set:     SIM_SETTINGS with the experiment-specific fields already changed by run_simulations.R
## ---------------------------------------------------------------------------
run_replicate <- function(rep, scenario, p, n, set = SIM_SETTINGS, experiment = "main", setting = 0, setting_label = "",
                          methods = METHODS[[experiment]]) {
  ## setting / setting_label only label the output rows (index and name of the experiment-specific setting)

  set.seed(truth_seed(rep));                                truth0 <- gen_truth(scenario, p, set)
  ## The data seed does not depend on the setting: within an experiment the settings (error variance, confounder
  ## scale, hyperparameters) are compared on the same underlying draws (paired), and the "base" setting of every
  ## experiment reproduces the data of the main study for the same (scenario, p, n, replicate).
  dseed <- data_seed(scenario, p, n, rep, set);             set.seed(dseed)
  dat0  <- gen_data(truth0, n, set)
  std   <- standardize_data(dat0, truth0, set)
  X <- std$X; Y <- std$Y; truth <- std$truth; D <- truth$D
  ms <- marginal_stats(X, Y); attr(ms, "X") <- X
  Syy <- crossprod(Y) / n

  out <- list()
  fits <- list()                       # estimates of every method, saved next to the metrics (new metrics without a rerun)
  keep <- function(name, ...) fits[[name]] <<- list(...)
  add <- function(method, family, vals, secs = NA) {
    out[[length(out) + 1]] <<- data.frame(experiment = experiment, scenario = scenario, p = p, n = n, rep = rep,
                                          setting = setting, setting_label = setting_label,
                                          method = method, family = family,
                                          metric = names(vals), value = unname(vals),
                                          seconds = secs, stringsAsFactors = FALSE)
  }
  timed <- function(expr) { t0 <- proc.time()[["elapsed"]]; v <- expr; list(v = v, s = proc.time()[["elapsed"]] - t0) }
  sel_xy <- if (!is.null(truth$roles)) { s <- matrix(FALSE, p, p); s[truth$roles$Y_idx, truth$roles$X_idx] <- TRUE; s } else NULL
  sel_yy <- if (!is.null(truth$roles)) { s <- matrix(FALSE, p, p); s[truth$roles$Y_idx, truth$roles$Y_idx] <- TRUE; s & row(s) != col(s) } else NULL
  ## confounding: in the roles design the comparison is restricted to the outcomes (the block MR2 models)
  Sg <- truth0$Sigma_star                                         # confounding truth on the generating scale
  conf_block <- function(Z) if (is.null(truth$roles)) conf_metrics(Sg, Z, set$pip_thr) else
    conf_metrics(Sg[truth$roles$Y_idx, truth$roles$Y_idx], Z[truth$roles$Y_idx, truth$roles$Y_idx], set$pip_thr)

  ## competitors returning p-values, matrices [to, from] after `orient`
  add_pval <- function(name, est, pval, failed, secs) {
    keep(name, est = est, pval = pval, failed = failed, seconds = secs)
    add(name, "graph",  graph_metrics(truth$graph, pval < set$alpha, 1 - pval), secs)
    add(name, "effect", effect_metrics(truth$A, est), secs)
    if (!is.null(sel_xy)) {
      add(name, "graph_xy",  graph_metrics(truth$graph, pval < set$alpha, 1 - pval, sel_xy), secs)
      add(name, "effect_xy", effect_metrics(truth$A, est, sel_xy), secs)
      add(name, "graph_yy",  graph_metrics(truth$graph, pval < set$alpha, 1 - pval, sel_yy), secs)
    }
    add(name, "diagnostic", c(failed_fits = failed), secs)
  }
  ## Bayesian competitors returning inclusion probabilities, matrices [to, from]
  add_pip <- function(name, est, score, failed, secs, conf_score = NULL) {
    keep(name, est = est, score = score, conf_score = conf_score, failed = failed, seconds = secs)
    add(name, "graph",  graph_metrics(truth$graph, score > set$pip_thr, score), secs)
    add(name, "effect", effect_metrics(truth$A, est), secs)
    if (!is.null(sel_xy)) {
      add(name, "graph_xy",  graph_metrics(truth$graph, score > set$pip_thr, score, sel_xy), secs)
      add(name, "effect_xy", effect_metrics(truth$A, est, sel_xy), secs)
      add(name, "graph_yy",  graph_metrics(truth$graph, score > set$pip_thr, score, sel_yy), secs)
    }
    if (!is.null(conf_score)) add(name, "confounding", conf_block(abs(conf_score)), secs)
    add(name, "diagnostic", c(failed_fits = failed), secs)
  }
  ## competitors returning their own edge selection and effect estimates (GrIVET), matrices [to, from]; score = |effect|
  add_sel <- function(name, est, score, edge, failed, secs) {
    keep(name, est = est, score = score, edge = edge, failed = failed, seconds = secs)
    add(name, "graph",  graph_metrics(truth$graph, edge, score), secs)
    add(name, "effect", effect_metrics(truth$A, est), secs)
    if (!is.null(sel_xy)) {
      add(name, "graph_xy",  graph_metrics(truth$graph, edge, score, sel_xy), secs)
      add(name, "effect_xy", effect_metrics(truth$A, est, sel_xy), secs)
      add(name, "graph_yy",  graph_metrics(truth$graph, edge, score, sel_yy), secs)
    }
    add(name, "diagnostic", c(failed_fits = failed), secs)
  }
  add_rgm <- function(name, fit, secs, conf = FALSE, iv = FALSE) {
    keep(name, AEst = fit$AEst, GammaEst = fit$GammaEst, ZEst = fit$ZEst, BEst = fit$BEst, PhiEst = fit$PhiEst,
         SigmaEst = fit$SigmaEst, AccptA = fit$AccptA, AccptB = fit$AccptB, PropVarA = fit$PropVarA,
         diagnostics = fit$diagnostics, warnings = fit$warnings, gibbsB = fit$settings$gibbsB, seconds = secs,
         halves = { GP <- fit$GammaPst; h <- dim(GP)[3] %/% 2
                    list(first = apply(GP[, , 1:h, drop = FALSE], c(1, 2), mean), second = apply(GP[, , (h + 1):dim(GP)[3], drop = FALSE], c(1, 2), mean)) })
    add(name, "graph",  graph_metrics(truth$graph, fit$GammaEst > set$pip_thr, fit$GammaEst), secs)
    add(name, "effect", effect_metrics(truth$A, fit$AEst), secs)
    if (!is.null(sel_xy)) {
      add(name, "graph_xy",  graph_metrics(truth$graph, fit$GammaEst > set$pip_thr, fit$GammaEst, sel_xy), secs)
      add(name, "effect_xy", effect_metrics(truth$A, fit$AEst, sel_xy), secs)
      add(name, "graph_yy",  graph_metrics(truth$graph, fit$GammaEst > set$pip_thr, fit$GammaEst, sel_yy), secs)
    }
    if (conf) add(name, "confounding", conf_block(fit$ZEst), secs)
    if (iv)   add(name, "iv_selection", iv_metrics(truth$B, fit$PhiEst), secs)
    add(name, "diagnostic", rgm_diagnostics(fit, Syy), secs)
  }
  ## A method that errors is recorded as a failed fit and its error message is appended to the log file
  ## getOption("sim_error_log") (set by run_simulations.R), because messages of forked workers are not shown.
  ## For the competitors, a failed fit is stored as "no network found" (on_fail: every effect 0, every p-value 1 or
  ## inclusion probability 0) with seconds = NA and failed_fits = 1. make_figures_tables.R and run_instrument_strength.R
  ## exclude these fits: a method is summarised over its successful fits, and the number is stated in the captions.
  ## Typical cause: MR2G needs genome-wide significant instruments (p < 5e-8, fixed in the package) for every trait,
  ## which weak instruments or a small n may not provide.
  null_pval <- function(name) { est <- matrix(0, p, p); diag(est) <- NA; pv <- matrix(1, p, p); diag(pv) <- NA
    add_pval(name, est, pv, failed = 1, secs = NA) }
  null_pip <- function(name, conf = FALSE) {
    if (is.null(truth$roles)) {                                    # by_outcome(): every off-diagonal entry estimated
      est <- matrix(0, p, p); diag(est) <- NA
      return(add_pip(name, est, est, failed = 1, secs = NA)) }
    Xi <- truth$roles$X_idx; Yi <- truth$roles$Y_idx
    est <- matrix(NA_real_, p, p); est[Yi, Xi] <- 0
    score <- matrix(NA_real_, p, p); score[Yi, Xi] <- 0
    if (name == "MrDAG") { score[Yi, Yi] <- 0; score[Xi, Xi] <- 0; diag(score) <- NA }
    cs <- NULL; if (conf) { cs <- matrix(NA_real_, p, p); cs[Yi, Yi] <- 0; diag(cs) <- 1 }
    add_pip(name, est, score, failed = 1, secs = NA, conf_score = cs) }
  null_sel <- function(name) { est <- matrix(0, p, p); diag(est) <- NA
    add_sel(name, est, est, est != 0, failed = 1, secs = NA) }
  safely <- function(name, expr_fun, on_fail = NULL) {
    set.seed(dseed + METHOD_SEED[[name]])
    r <- tryCatch(expr_fun(), error = function(e) {
      msg <- sprintf("[%s] %s %s p=%d n=%g setting=%s rep=%d: %s failed: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                     experiment, scenario, p, n, setting_label, rep, name, conditionMessage(e))
      message(msg)
      lf <- getOption("sim_error_log"); if (!is.null(lf)) cat(msg, "\n", file = lf, append = TRUE, sep = "")
      if (is.null(on_fail)) add(name, "diagnostic", c(failed_fits = 1)) else on_fail()
      NULL })
    invisible(r)
  }
  Dfull <- matrix(1, nrow = p, ncol = ncol(X))

  for (m in methods) switch(m,
    "MendelianRandomization" = safely(m, function() {
      r <- timed(fit_mr_package(X, Y, D, set))
      for (mm in names(r$v$est)) add_pval(paste0("MR-", mm), t(r$v$est[[mm]]), t(r$v$pval[[mm]]), r$v$failed, r$s) },
      on_fail = function() for (mm in c("SimpleMedian", "WeightedMedian", "IVW")) null_pval(paste0("MR-", mm))),
    "OneSampleMR" = safely(m, function() { r <- timed(fit_onesamplemr(X, Y, D)); add_pval(m, t(r$v$est), t(r$v$pval), r$v$failed, r$s) }),
    "mrbayes"     = safely(m, function() { r <- timed(fit_mrbayes(X, Y, D, set));  add_pval(m, t(r$v$est), t(r$v$pval), r$v$failed, r$s) }),
    "MVMR"        = safely(m, function() { r <- timed(fit_mvmr(ms, D));           add_pval(m, r$v$est, r$v$pval, r$v$failed, r$s) },
                           on_fail = function() null_pval(m)),
    "MR2G"        = safely(m, function() { r <- timed(suppressWarnings(fit_mr2g(ms, Y, D, set))); add_pval(m, r$v$est, r$v$pval, r$v$failed, r$s) },
                           on_fail = function() null_pval(m)),
    ## MR2: roles design only (METHODS); MrDAG: native split in the roles design, one outcome at a time otherwise
    "MR2"         = safely(m, function() { if (is.null(truth$roles)) stop("MR2 needs the exposure/outcome split of the roles design")
      r <- timed(fit_mr2(ms, D, truth$roles, set)); add_pip(m, r$v$est, r$v$score, r$v$failed, r$s, conf_score = r$v$resid_cor) },
      on_fail = function() null_pip(m, conf = !is.null(truth$roles))),
    "MrDAG"       = safely(m, function() {
      r <- timed(if (is.null(truth$roles)) by_outcome(fit_mrdag, ms, D, set) else fit_mrdag(ms, D, truth$roles, set))
      add_pip(m, r$v$est, r$v$score, r$v$failed, r$s) },
      on_fail = function() null_pip(m)),
    "GrIVET"      = safely(m, function() {
      r <- timed(fit_grivet(dat0$X, dat0$Y, if (is.null(truth$scale)) rep(1, p) else truth$scale$sY, set))
      add_sel(m, r$v$est, r$v$score, r$v$edge, r$v$failed, r$s) },
      on_fail = function() null_sel(m)),
    "MR.RGM"        = safely(m, function() { r <- timed(fit_rgm(X, Y, D, "SSSL", set, selectB = FALSE));     add_rgm(m, r$v, r$s, conf = TRUE) }),
    "MR.RGM_NoConf" = safely(m, function() { r <- timed(fit_rgm(X, Y, D, "diagonal", set, selectB = FALSE)); add_rgm(m, r$v, r$s) }),
    "MR.RGM+"       = safely(m, function() { r <- timed(suppressWarnings(fit_rgm(X, Y, Dfull, "SSSL", set, selectB = TRUE)));
      add_rgm(m, r$v, r$s, conf = TRUE, iv = TRUE) }),
    stop("unknown method ", m))

  ## instrument strength of the simulated traits (population R^2 of each trait on its own instruments)
  add("data", "instruments", pop_r2(truth0, n))
  res <- do.call(rbind, out)
  ## Truth and estimates of this replicate (saved by run_simulations.R in fits_*.rds): every metric of the paper can be
  ## recomputed, and new metrics added, from these files without rerunning any method.
  attr(res, "fits") <- list(experiment = experiment, scenario = scenario, p = p, n = n, rep = rep, setting_label = setting_label,
                            data_seed = dseed,
                            truth = list(A = truth$A, graph = truth$graph, B = truth$B, D = truth$D, roles = truth$roles,
                                         Sigma_star = truth$Sigma_star, Sigma_star_gen = truth0$Sigma_star, scale = truth$scale),
                            methods = fits)
  res
}

## ---------------------------------------------------------------------------
## Instrument-strength study (run_instrument_strength.R): one replicate of run_replicate() with the instrument effect
## and the number of instruments per trait of the strength level. Truth and data seeds are those of the main study and
## U, E are drawn before X (gen_data), so within a replicate the graph, the confounders and the errors are identical
## across strength levels: only the instrument effects differ (paired comparison). Returns the long data.frame with the
## column `strength` (instead of experiment / setting) and the estimates in attr(, "fits").
## ---------------------------------------------------------------------------
STRENGTH <- c(Strong = 1, Moderate = 0.5, Weak = 0.3, Weak10 = sqrt(0.027))   # effect of each instrument
IVS      <- c(Strong = 3, Moderate = 3,   Weak = 3,   Weak10 = 10)            # instruments per trait
## "Weak10": the same total instrument strength as "Weak" (10 x 0.027 = 3 x 0.09) spread over 10 SNPs per trait
strength_replicate <- function(rep, scenario, p, n, strength, set = SIM_SETTINGS, methods = METHODS$strength) {
  set$iv_effect <- STRENGTH[[strength]]; set$ivs_per_trait <- IVS[[strength]]
  r <- run_replicate(rep, scenario, p, n, set = set, experiment = "strength", setting_label = strength, methods = methods)
  f <- attr(r, "fits")
  out <- data.frame(scenario = scenario, p = p, n = n, strength = strength, rep = rep, method = r$method, family = r$family,
                    metric = r$metric, value = r$value, seconds = r$seconds, stringsAsFactors = FALSE)
  attr(out, "fits") <- list(scenario = scenario, p = p, n = n, strength = strength, rep = rep, data_seed = f$data_seed,
                            truth = f$truth, methods = f$methods)
  out
}

## ---------------------------------------------------------------------------
## Runtime benchmark (Runtime Analysis/RunTimeAnalysis.R): the data set for p traits (scale-free network, n = 30,000,
## replicate-1 truth, data seed 9e6 + p) and the call of every method, evaluated by microbenchmark() in the
## environment that holds X, Y, D, Dfull, ms and set (runtime_data()).
## ---------------------------------------------------------------------------
runtime_data <- function(p, n = 30000, set = SIM_SETTINGS) {
  set.seed(truth_seed(1));         truth0 <- gen_truth("scalefree", p, set)
  set.seed(as.integer(9e6 + p));   dat0   <- gen_data(truth0, n, set)
  std <- standardize_data(dat0, truth0, set)
  ms  <- marginal_stats(std$X, std$Y); attr(ms, "X") <- std$X
  list(X = std$X, Y = std$Y, D = std$truth$D, Dfull = matrix(1, p, ncol(std$X)), ms = ms,
       X_raw = dat0$X, Y_raw = dat0$Y, sY = std$truth$scale$sY)          # GrIVET: unstandardized data (fit_grivet())
}
## a failing competitor returns NULL instead of stopping the benchmark (the benchmark design has no exposure/outcome
## split, so MR2 is not timed; its seconds per fit in the roles design are in the results of that experiment)
runtime_safe <- function(expr) tryCatch(expr, error = function(e) { message("   failed: ", conditionMessage(e)); NULL })
RUNTIME_EXPRS <- alist(
  `MR.RGM`                 = fit_rgm(X, Y, D,     "SSSL",     set, selectB = FALSE),
  `MR.RGM+`                = suppressWarnings(fit_rgm(X, Y, Dfull, "SSSL", set, selectB = TRUE)),
  `MR.RGM_NoConf`          = fit_rgm(X, Y, D,     "diagonal", set, selectB = FALSE),
  `MVMR`                   = runtime_safe(fit_mvmr(ms, D)),
  `MR2G`                   = runtime_safe(suppressWarnings(fit_mr2g(ms, Y, D, set))),
  `MendelianRandomization` = runtime_safe(fit_mr_package(X, Y, D, set)),
  `OneSampleMR`            = runtime_safe(fit_onesamplemr(X, Y, D)),
  `mrbayes`                = runtime_safe(fit_mrbayes(X, Y, D, set)),
  `MrDAG`                  = runtime_safe(by_outcome(fit_mrdag, ms, D, set)),
  `GrIVET`                 = runtime_safe(fit_grivet(X_raw, Y_raw, sY, set)))
