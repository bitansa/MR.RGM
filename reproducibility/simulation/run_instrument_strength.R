###############################################################################
## run_instrument_strength.R  --  every method under strong and weak instruments
##                                 (online appendix; complements the three main cases)
##
## Question: when does estimating the instrument map (MR.RGM+, D = 1 with selection) cost more than it gains over
## using a known map without selection (MR.RGM)? The main study uses strong instruments (generating effect 1,
## population R^2 of a trait on its own instruments about 0.2). Real cis-eQTL instruments are much weaker
## (R^2 about 0.01-0.08 per gene). Here the instrument effect is varied while everything else is as in Case 2
## (small-world, no pleiotropy) and Case 3 (small-world, pleiotropy):
##   strength:  strong (effect 1), moderate (0.5), weak (0.3), all with 3 SNPs per trait, and weak with
##              10 instruments per response (same total strength as weak; cis-eQTL-like)
##   cases:     small-world without / with horizontal pleiotropy, p = 5 and 10 (as in the main study)
##   n:         500, 1,000, 10,000, 30,000 (as in the main study);  100 replicates;  every method of the main study
##              (METHODS$strength in sim_functions.R; same MCMC settings and package defaults)
## One replicate = strength_replicate() in sim_functions.R, i.e. run_replicate() of the main study with the instrument
## effect and number of the strength level. Truth and data seeds are those of the main study and U, E are drawn before
## X, so for a given replicate the graph, the confounders and the errors are identical across strength levels: only the
## instrument effects differ (paired comparison). Data are standardized before fitting, as in the main study.
##
## Reported (same metrics as the main study, plus the real-data diagnostics):
##   causal network AUC, MAD of A, confounding AUC, instrument-selection AUC (MR.RGM+),
##   population R^2 of each trait on its own instruments, mu^2 = n R^2 / (1 - R^2) and F (medians over traits),
##   error-variance ratio = estimated error variance / second moment of the trait (median and max over traits;
##   RGM() warns when the max exceeds 5), the correlation of the causal PIPs from the two halves of the chain,
##   and the spectral radius of the estimated A.
##
## Usage (from the simulation folder):
##   source("run_instrument_strength.R")      # several days on 9 cores; restartable (finished settings skipped)
##   STRENGTH_PLOT_ONLY <- TRUE; source("run_instrument_strength.R")   # redraw from results_strength/
## Optional before sourcing: ST_SCENARIOS, ST_P_VALUES, ST_N_VALUES, ST_LEVELS (subsets of the grid);
##   STRENGTH_RUN_ONLY <- TRUE: run (and save) only, without the summaries, tables and figures (one setting per cluster task)
## Outputs: results_strength/*.rds, figures/strength_*.pdf (+ png), tables/strength_*.tex (network AUC and effect MAD
##          of every method, confounding/IV/error ratio of the MR.RGM family, and strength_diag_p=*.tex with the
##          chain-half correlation and the fraction of fits that trigger the package's error-variance warning),
##          tables/strength_summary.csv
###############################################################################

source("sim_functions.R")
suppressPackageStartupMessages({ library(parallel); library(ggplot2); library(scales) })
if (!exists("N_CORES")) N_CORES <- max(1, parallel::detectCores() - 1)

## Strength levels: STRENGTH (effect of each instrument) and IVS (instruments per trait) in sim_functions.R
ST_LABEL <- c(Strong = "Strong", Moderate = "Moderate", Weak = "Weak", Weak10 = "Weak, 10 instruments")
## "Weak, 10 instruments": the same total instrument strength as "Weak" (10 x 0.027 = 3 x 0.09) spread over 10 instruments
## per trait, as for cis-eQTLs (several weak SNPs per gene); MR.RGM+ then has many more SNP-trait effects to estimate.
ST_SCEN  <- STRENGTH_GRID$scen                                # cases, network sizes and sample sizes of the main study
ST_P     <- STRENGTH_GRID$p                                   # (STRENGTH_GRID in sim_functions.R)
ST_N     <- STRENGTH_GRID$n
ST_LEV   <- names(STRENGTH)
out_dir  <- "results_strength"
ST_METHODS <- METHODS$strength                                 # every method (sim_functions.R)
if (exists("SKIP_METHODS")) ST_METHODS <- setdiff(ST_METHODS, SKIP_METHODS)   # fitted in a separate pass (cluster/add_methods.R)
FIG_DIR  <- "figures"; TAB_DIR <- "tables"
## STRENGTH_SMOKE <- TRUE: 2 replicates, short chains, p = 5, n = 500, separate output folders (end-to-end check)
SMOKE_ST <- exists("STRENGTH_SMOKE") && isTRUE(STRENGTH_SMOKE)
if (SMOKE_ST) { ST_P <- 5; ST_N <- 500; out_dir <- "results_strength_smoke"; FIG_DIR <- "figures_smoke"; TAB_DIR <- "tables_smoke" }
dir.create(out_dir, showWarnings = FALSE); dir.create(file.path(FIG_DIR, "png"), recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, showWarnings = FALSE)
options(sim_error_log = file.path(out_dir, "errors.log"))      # error messages of failed fits

## ---------------------------------------------------------------------------
## Run (settings whose result file exists are skipped)
## ---------------------------------------------------------------------------
if (!(exists("STRENGTH_PLOT_ONLY") && isTRUE(STRENGTH_PLOT_ONLY))) {
  set <- SIM_SETTINGS
  if (SMOKE_ST) { set$n_rep <- 2; set$mcmc <- list(nIter = 2000, nBurnin = 500, Thin = 2)
                  set$mr2g <- modifyList(set$mr2g, list(num_pert = 5, num_pert_M = 5)); set$mr_boot <- 200
                  set$mrbayes <- list(n.chains = 3, n.burn = 200, n.iter = 800)
                  set$mr2 <- modifyList(set$mr2, list(niter = 600, burnin = 200, thin = 2))
                  set$mrdag <- modifyList(set$mrdag, list(niter = 600, burnin = 200, thin = 2)) }
  run_scen <- if (exists("ST_SCENARIOS")) ST_SCENARIOS else ST_SCEN; run_p <- if (exists("ST_P_VALUES")) ST_P_VALUES else ST_P
  run_n    <- if (exists("ST_N_VALUES")) ST_N_VALUES else ST_N;      run_lev <- if (exists("ST_LEVELS")) ST_LEVELS else ST_LEV
  for (sc in run_scen) for (p_ in run_p) for (st in run_lev) for (n in run_n) {
    f <- file.path(out_dir, sprintf("strength_%s_%s_p%d_n%d.rds", sc, st, p_, n))
    if (check_existing_result(f)) next
    message(sprintf("[%s] %s  %s instruments  p = %d  n = %d  (%d replicates on %d cores; methods: %s)",
                    format(Sys.time(), "%H:%M:%S"), sc, st, p_, n, set$n_rep, N_CORES, paste(ST_METHODS, collapse = ", ")))
    t0  <- Sys.time()
    res <- run_parallel(seq_len(set$n_rep), strength_replicate, scenario = sc, p = p_, n = n, strength = st,
                        set = set, methods = ST_METHODS, cores = N_CORES,
                        timeout_min = if (exists("REP_TIMEOUT_MIN")) REP_TIMEOUT_MIN else 240,
                        label = sprintf("strength %s %s p=%d n=%d", sc, st, p_, n))
    bad <- vapply(res, function(r) inherits(r, "try-error") || is.null(r), logical(1))
    if (any(bad)) warning(sprintf("%d replicate(s) failed in %s %s n=%d", sum(bad), sc, st, n))
    if (all(bad)) next
    fits <- lapply(res[!bad], attr, "fits")
    df <- do.call(rbind, lapply(res[!bad], function(r) { attr(r, "fits") <- NULL; r }))
    attr(df, "settings") <- c(set, list(strength = STRENGTH, ivs = IVS))
    attr(df, "methods")  <- ST_METHODS
    save_rds_atomic(list(settings = attr(df, "settings"), methods = ST_METHODS, fits = fits), sub("/strength_", "/fits_strength_", f, fixed = TRUE))
    attr(df, "session_info") <- utils::sessionInfo()
    attr(df, "package_versions") <- package_versions()
    attr(df, "sampler") <- SAMPLER_TAG
    attr(df, "MR.RGM_build") <- unlist(utils::packageDescription("MR.RGM")[c("Version", "Packaged", "RemoteSha", "Built")])
    save_rds_atomic(df, f)                                       # written last: its presence marks the setting as done
    message(sprintf("   done in %.1f min -> %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), f))
  }
}

## ---------------------------------------------------------------------------
## Summary, tables and figures
## ---------------------------------------------------------------------------
if (!(exists("STRENGTH_RUN_ONLY") && isTRUE(STRENGTH_RUN_ONLY))) {
files <- list.files(out_dir, pattern = "^strength_.*\\.rds$", full.names = TRUE)
if (!length(files)) stop("No results in ", normalizePath(out_dir, mustWork = FALSE), "; run the study first.")
## each result file plus the effect error on the true edges / non-edges and the All-zero reference, computed from the
## estimates saved next to it (effect_split_rows() in sim_functions.R)
res <- do.call(rbind, lapply(files, function(f) {
  d  <- as.data.frame(readRDS(f)); for (a in c("settings", "session_info", "package_versions", "methods", "sampler", "MR.RGM_build")) attr(d, a) <- NULL
  ff <- sub("/strength_", "/fits_strength_", f, fixed = TRUE)
  if (!file.exists(ff)) return(d)
  s  <- effect_split_rows(ff)
  k  <- d[rep(1, nrow(s)), c("scenario", "p", "n", "strength"), drop = FALSE]
  s  <- cbind(k, s, seconds = NA_real_); rownames(s) <- NULL
  rbind(d, s[, names(d)]) }))
res <- res[!is.na(res$value), , drop = FALSE]
## GrIVET fits whose cluster job ended before they returned are added as failed fits with cause "time", and the cause
## of every failed GrIVET fit is read from the error log (add_unfinished() and fail_causes() in sim_functions.R)
UNFINISHED <- NULL
if (any(res$method == "GrIVET")) { res <- add_unfinished(res, "GrIVET", c("scenario", "p", "n", "strength")); UNFINISHED <- attr(res, "unfinished"); rownames(res) <- NULL }
GRIVET_CAUSE <- fail_causes(file.path(out_dir, "errors.log"), "GrIVET")
fit_cause <- function(d) {                                       # d: failed-fit rows; "error": no message in the log
  gc <- setNames(GRIVET_CAUSE$cause, paste(GRIVET_CAUSE$scenario, GRIVET_CAUSE$p, GRIVET_CAUSE$n, GRIVET_CAUSE$setting_label, GRIVET_CAUSE$rep, sep = "|"))
  uf <- if (is.null(UNFINISHED)) character(0) else paste(UNFINISHED$scenario, UNFINISHED$p, UNFINISHED$n, UNFINISHED$strength, UNFINISHED$rep, sep = "|")
  k  <- paste(d$scenario, d$p, d$n, d$strength, d$rep, sep = "|")
  ifelse(d$method == "GrIVET", ifelse(k %in% uf, "time", ifelse(k %in% names(gc), gc[k], "error")), ifelse(d$method == "MR2G", "significance", "error"))
}
## A fit that stops with an error (MR2G when a trait has no genome-wide significant instrument; any other method) is
## stored as "no network found" with seconds = NA and failed_fits = 1; as in make_figures_tables.R, every table and
## figure summarises the successful fits only (entries over fewer fits are marked; the counts are in
## tables/strength_failures.tex and tables/strength_figure_notes.txt).
## With fewer than MIN_FITS successful fits a method is not summarised in that setting: the table entry reads
## "failed (k/100)" and the method is left out of that panel of the figures (listed in tables/strength_figure_notes.txt).
fit_key  <- function(d) paste(d$scenario, d$p, d$n, d$strength, d$method, d$rep, sep = "|")
fail_row <- res$family == "diagnostic" & res$metric == "failed_fits" & is.na(res$seconds) & res$value >= 1
fail_df  <- res[fail_row, , drop = FALSE]                        # one row per failed fit
res$failed_fit <- fit_key(res) %in% unique(fit_key(fail_df))
res <- res[!res$failed_fit | res$family == "diagnostic", , drop = FALSE]
if (nrow(fail_df)) message(nrow(fail_df), " failed fits excluded from the summaries")
MIN_FITS   <- 10
N_REP_ST   <- max(table(paste(res$scenario, res$p, res$n, res$strength, res$method)[res$method == "MR.RGM" & res$family == "graph" & res$metric == "AUC"]))
st_key     <- function(d) paste(d$scenario, d$p, d$n, d$strength, d$method, sep = "|")
FAIL_COUNT <- table(st_key(fail_df))
n_failed   <- function(key) if (key %in% names(FAIL_COUNT)) as.integer(FAIL_COUNT[[key]]) else 0L
fail_df$cause <- if (nrow(fail_df)) fit_cause(fail_df) else character(0)
if (any(fail_df$method == "GrIVET")) {
  cs <- fail_df$cause[fail_df$method == "GrIVET"]
  message("GrIVET failed fits by cause: ", paste(names(table(cs)), table(cs), sep = " ", collapse = ", "))
}
TIME_COUNT <- table(st_key(fail_df)[fail_df$cause == "time"])
n_time     <- function(key) if (key %in% names(TIME_COUNT)) as.integer(TIME_COUNT[[key]]) else 0L
low_keys   <- names(FAIL_COUNT)[N_REP_ST - as.integer(FAIL_COUNT) < MIN_FITS]
res <- res[!(st_key(res) %in% low_keys) | res$family == "diagnostic", , drop = FALSE]
message(sprintf("Read %d result files (%d rows).", length(files), nrow(res)))
key <- paste(res$scenario, res$p, res$strength, res$n, res$method, res$family, res$metric, sep = "|")
summ <- do.call(rbind, lapply(split(res, key), function(d)
  data.frame(d[1, c("scenario", "p", "strength", "n", "method", "family", "metric")],
             mean = mean(d$value), sd = if (nrow(d) > 1) sd(d$value) else NA_real_, median = median(d$value), reps = nrow(d))))
rownames(summ) <- NULL
write.csv(summ, file.path(TAB_DIR, "strength_summary.csv"), row.names = FALSE)

## case_lab is displayed text only (tables, facet labels, figure notes): converted with display_text() (sim_functions.R)
case_lab <- display_text(c(smallworld = "No pleiotropy", smallworld_pleio = "Pleiotropy"))
res$Case     <- factor(case_lab[res$scenario], levels = case_lab)
res$Strength <- factor(ST_LABEL[res$strength], levels = ST_LABEL[names(STRENGTH)])
res$N        <- factor(res$n, levels = ST_N, labels = paste0("n = ", comma(ST_N)))
## Order, labels and colours of the methods: identical to make_figures_tables.R
method_order <- c("MR.RGM", "MR.RGM+", "MR.RGM_NoConf", "MR2G", "MVMR", "MR2", "MrDAG", "GrIVET",
                  "MR-SimpleMedian", "MR-WeightedMedian", "MR-IVW", "OneSampleMR", "mrbayes", "All-zero")
method_col   <- c("MR.RGM" = "#0072B2", "MR.RGM+" = "#E69F00", "MR.RGM_NoConf" = "#56B4E9",
                  "MR2G" = "#D55E00", "MVMR" = "#009E73", "MR2" = "#CC79A7", "MrDAG" = "#F0E442", "GrIVET" = "#882255",
                  "MR-SimpleMedian" = "#999999", "MR-WeightedMedian" = "#666666", "MR-IVW" = "#BBBBBB",
                  "OneSampleMR" = "#000000", "mrbayes" = "#8C510A", "All-zero" = "#DDDDDD")   # All-zero: reference (estimate 0)
method_tex   <- c("MR.RGM" = "\\texttt{MR.RGM}", "MR.RGM+" = "\\texttt{MR.RGM+}", "MR.RGM_NoConf" = "\\texttt{MR.RGM\\_NoConf}",
                  "MR2G" = "\\texttt{MR2G}", "MVMR" = "\\texttt{MVMR}", "MR2" = "\\texttt{MR2}", "MrDAG" = "\\texttt{MrDAG}", "GrIVET" = "\\texttt{GrIVET}",
                  "MR-SimpleMedian" = "\\texttt{MR-SimpleMedian}", "MR-WeightedMedian" = "\\texttt{MR-WeightedMedian}",
                  "MR-IVW" = "\\texttt{MR-IVW}", "OneSampleMR" = "\\texttt{OneSampleMR}", "mrbayes" = "\\texttt{mrbayes}",
                  "All-zero" = "\\textit{All-zero}")
method_shape <- c("MR.RGM" = 16, "MR.RGM+" = 17, "MR.RGM_NoConf" = 15, "MR2G" = 18, "MVMR" = 8, "MR2" = 7, "MrDAG" = 3,
                  "GrIVET" = 4, "MR-SimpleMedian" = 1, "MR-WeightedMedian" = 2, "MR-IVW" = 0, "OneSampleMR" = 5,
                  "mrbayes" = 6, "All-zero" = 13)
method_lty   <- c("MR.RGM" = "solid", "MR.RGM+" = "solid", "MR.RGM_NoConf" = "solid", "MR2G" = "longdash", "MVMR" = "longdash",
                  "MR2" = "longdash", "MrDAG" = "dotdash", "GrIVET" = "dotdash", "MR-SimpleMedian" = "dashed",
                  "MR-WeightedMedian" = "dashed", "MR-IVW" = "dashed", "OneSampleMR" = "dotted", "mrbayes" = "dotted",
                  "All-zero" = "twodash")
st_methods   <- intersect(method_order, unique(res$method))                    # every method in the result files
theme_st <- theme_bw(base_size = 16) +
  theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 22, face = "bold"), axis.text = element_text(size = 18),
        strip.text = element_text(size = 18, face = "bold"), legend.title = element_blank(),
        legend.text = element_text(size = 18), legend.position = "top", panel.grid.minor = element_blank())
## Print-size style: off unless a figure driver has defined
## THESIS_STYLE (then: theme_bw(base_size = 9), sizes from the driver's CSV, cairo_pdf, no png). The default figures are unaffected.
STYLE  <- if (exists("THESIS_STYLE")) THESIS_STYLE else NULL
THESIS <- !is.null(STYLE)
if (THESIS) theme_st <- theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"), legend.title = element_blank(), legend.position = "top",
        legend.margin = margin(0, 0, 0, 0), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))     # four strength levels per small panel
save_fig <- function(g, file, w = 18, h = 9) {
  if (THESIS) {                                  # narrow print-size panels: vertical tick labels, so they do not run into the next panel
    if (STYLE$size(file, w, h)[1] < 3.5)
      g <- suppressMessages(g + scale_x_discrete(labels = identity)) +    # one-line labels ("Weak, 10 instruments")
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    return(STYLE$save(g, file.path(FIG_DIR, file), w, h))
  }
  ggsave(file.path(FIG_DIR, file), g, width = w, height = h, units = "in")
  ggsave(file.path(FIG_DIR, "png", sub("\\.pdf$", ".png", file)), g, width = w, height = h, units = "in", dpi = 300)
}
ST_FIG_NOTES <- character()
## caption notes of a figure: methods over fewer than N_REP_ST fits, and those left out (fewer than MIN_FITS successful fits)
st_notes <- function(p_, fam, met, file, methods) {
  for (m in methods) {
    keys <- unique(paste(sub("\\|[^|]*$", "", st_key(res[res$p == p_ & res$family == fam & res$metric == met, ])), m, sep = "|"))
    k    <- vapply(keys, n_failed, integer(1)); keys <- keys[k > 0]; k <- k[k > 0]
    if (length(keys)) {
      z <- do.call(rbind, strsplit(keys, "|", fixed = TRUE))
      txt <- sprintf("%s, n = %s: failed %d/%d%s", case_lab[z[, 1]], format(as.numeric(z[, 3]), big.mark = ",", scientific = FALSE), k, N_REP_ST,
                     ifelse(N_REP_ST - k < MIN_FITS, " (not shown)", sprintf(" (summarized over %d fits)", N_REP_ST - k)))
      txt <- paste0(ST_LABEL[z[, 4]], ", ", txt)
      ST_FIG_NOTES <<- c(ST_FIG_NOTES, paste0(file, ": ", m, " did not return a result in some replicates; ", paste(txt, collapse = "; "), "."))
    }
  }
}
## boxplots per strength level (the MR.RGM family: confounding, instrument selection, error-variance ratio)
st_box <- function(p_, fam, met, ylab, title, file, ref = NULL, log_y = FALSE, methods = c("MR.RGM", "MR.RGM+")) {
  d <- res[res$p == p_ & res$family == fam & res$metric == met & res$method %in% methods, ]
  if (!nrow(d)) return(invisible(NULL))
  methods <- intersect(methods, unique(d$method))
  d$Method <- factor(d$method, levels = methods, labels = display_text(methods))     # labels: displayed method names
  st_notes(p_, fam, met, file, methods)
  g <- ggplot(d, aes(Strength, value, fill = Method))
  if (!is.null(ref)) g <- g + geom_hline(yintercept = ref, linetype = "dashed", colour = "grey50")
  g <- g + geom_boxplot(colour = "grey15", outlier.size = 1.5, position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = setNames(method_col[methods], display_text(methods))) +
    facet_grid(Case ~ N, labeller = labeller(Case = label_wrap_gen(12))) +                 # case labels on two lines
    scale_x_discrete(labels = function(x) sub(", ", ",\n", x, fixed = TRUE)) +             # two-line "Weak,\n10 instruments"
    labs(title = display_text(title), x = "Instrument Strength", y = display_text(ylab)) + theme_st
  if (log_y) g <- g + scale_y_log10()
  if (length(methods) == 1) g <- g + guides(fill = "none")
  save_fig(g, file)
}
## every method: mean (+- SD) over the successful fits, one column of points per method (so that no two methods
## overlap), the four instrument designs side by side within a method (colour and symbol); panels: sample size
## (rows) by case (columns); the MR.RGM family is shaded
st_col   <- c(Strong = "#0072B2", Moderate = "#009E73", Weak = "#E69F00", Weak10 = "#D55E00")
st_shape <- c(Strong = 16, Moderate = 17, Weak = 15, Weak10 = 18)
st_line <- function(p_, fam, met, ylab, title, file, ref = NULL, log_y = FALSE) {
  d <- summ[summ$p == p_ & summ$family == fam & summ$metric == met & summ$method %in% st_methods & !is.na(summ$mean), ]
  if (!nrow(d)) return(invisible(NULL))
  st_notes(p_, fam, met, file, st_methods)                             # also methods left out at every n
  lv         <- intersect(st_methods, unique(d$method))
  d$Method   <- factor(d$method, levels = lv, labels = display_text(lv))          # labels: displayed method names
  d$Case     <- factor(case_lab[d$scenario], levels = case_lab)
  d$Strength <- factor(ST_LABEL[d$strength], levels = ST_LABEL[names(STRENGTH)])
  d$N        <- factor(d$n, levels = ST_N, labels = paste0("n = ", comma(ST_N)))
  d$lo <- d$mean - d$sd
  if (log_y) d$lo <- pmax(d$lo, d$mean / 10)                           # log scale: lower whisker kept positive
  n_fam <- sum(levels(d$Method) %in% display_text(c("MR.RGM", "MR.RGM+", "MR.RGM_NoConf")))
  g <- ggplot(d, aes(Method, mean, colour = Strength, shape = Strength)) + scale_x_discrete() +
    annotate("rect", xmin = 0.5, xmax = n_fam + 0.5, ymin = if (log_y) 0 else -Inf, ymax = Inf, fill = "grey92")
  if (!is.null(ref)) g <- g + geom_hline(yintercept = ref, linetype = "dashed", colour = "grey50")
  g <- g + geom_pointrange(aes(ymin = lo, ymax = mean + sd), position = position_dodge(width = 0.75),
                           size = if (THESIS) 0.18 else 0.4, linewidth = if (THESIS) 0.35 else 0.6) +
    scale_colour_manual(values = setNames(st_col, ST_LABEL[names(st_col)]), name = NULL) +
    scale_shape_manual(values = setNames(st_shape, ST_LABEL[names(st_shape)]), name = NULL) +
    facet_grid(N ~ Case) + labs(title = display_text(title), x = NULL, y = display_text(ylab)) + theme_st +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.x = element_blank())
  if (log_y) g <- g + scale_y_log10()
  if (THESIS) return(STYLE$save(g, file.path(FIG_DIR, file), 18, 12))     # method labels: no vertical relabelling
  save_fig(g, file, 18, 12)
}
for (p_ in ST_P) {
  st_line(p_, "graph",      "AUC",           "Area Under Curve (AUC)",        sprintf("Causal network recovery by instrument strength (p = %d)", p_), sprintf("strength_auc_p=%d.pdf", p_), ref = 0.5)
  st_line(p_, "effect",     "MeanAbsDev",    "Mean Absolute Deviation (MAD)", sprintf("Causal effect estimation by instrument strength (p = %d)", p_), sprintf("strength_mad_p=%d.pdf", p_), log_y = TRUE)   # log scale: MR.RGM_NoConf diverges with weak instruments
  st_box(p_, "confounding", "AUC",           "Area Under Curve (AUC)",        sprintf("Confounding recovery by instrument strength (p = %d)", p_), sprintf("strength_conf_p=%d.pdf", p_), ref = 0.5)
  st_box(p_, "iv_selection","AUC",           "Area Under Curve (AUC)",        sprintf("Instrument selection by MR.RGM+ by instrument strength (p = %d)", p_), sprintf("strength_iv_p=%d.pdf", p_), ref = 0.5, methods = "MR.RGM+")
  st_box(p_, "diagnostic",  "err_ratio_max", "Error Variance Ratio (max over traits)", sprintf("Weak-identification diagnostic (p = %d)", p_), sprintf("strength_errvar_p=%d.pdf", p_), ref = c(1, 5), log_y = TRUE)
}

## LaTeX tables: mean (SD) over replicates
fmt <- function(m, s) ifelse(is.na(m), "--", sprintf("%.3f (%s)", m, ifelse(is.na(s), "--", sprintf("%.3f", s))))
cell_val <- function(sc, p_, st, n_, meth, fam, met, what = "mean") {   # entries over fewer fits are marked $^\\dagger$
  k <- n_failed(paste(sc, p_, n_, st, meth, sep = "|"))
  if (what == "mean" && N_REP_ST - k < MIN_FITS)
    return(sprintf(if (n_time(paste(sc, p_, n_, st, meth, sep = "|")) == k) "not finished (%d/%d)" else "failed (%d/%d)", k, N_REP_ST))
  r <- summ[summ$scenario == sc & summ$p == p_ & summ$strength == st & summ$n == n_ & summ$method == meth & summ$family == fam & summ$metric == met, ]
  if (!nrow(r)) "--" else paste0(fmt(r$mean, r$sd), if (r$reps < N_REP_ST) "$^\\dagger$" else "")
}
## Every table at a readable font size; adjustbox only guards against a table wider than the text
tex_table <- function(file, caption, label, colspec, hdr, rows, size = "\\footnotesize") {
  L <- c("\\begin{table}[htbp]", "\\centering", size, "\\setlength{\\tabcolsep}{3.5pt}", sprintf("\\caption{%s}", caption),
         sprintf("\\label{%s}", label), "\\begin{adjustbox}{max width=\\textwidth, max totalheight=0.78\\textheight}", sprintf("\\begin{tabular}{%s}", colspec),
         "\\hline", paste(hdr, "\\\\"), "\\hline", "\\hline", rows, "\\end{tabular}", "\\end{adjustbox}", "\\end{table}")
  L <- gsub("(?<=[0-9]),(?=[0-9]{3}(?![0-9]))", "{,}", L, perl = TRUE)   # thousands separator that also renders correctly in math mode
  writeLines(display_text(L), file.path(TAB_DIR, file))          # displayed method names and wording (sim_functions.R)
}
## Two-row header: column labels "Quantity Method" are grouped under a common top cell, either by quantity
## (by = "quantity": top = quantity, bottom = method) or by method (by = "method": top = method, bottom = quantity);
## lead = leading columns without a group. The case (valid / invalid instruments) is a row above each block.
two_row_header <- function(lead, labels, by = c("quantity", "method")) {
  by <- match.arg(by)
  qty <- sub("\\\\$", "", sub(" (MR\\.RGM\\+?)$", "", labels))   # drop the "\" of a "\ " space before the method name
  meth <- sub("^.* (MR\\.RGM\\+?)$", "\\1", labels)
  top_lab <- if (by == "quantity") qty else meth; bot_lab <- if (by == "quantity") meth else qty
  r <- rle(top_lab)
  top <- vapply(seq_along(r$lengths), function(i) if (r$lengths[i] == 1) sprintf("\\textbf{%s}", r$values[i]) else
    sprintf("\\multicolumn{%d}{c|}{\\textbf{%s}}", r$lengths[i], r$values[i]), character(1))
  end <- length(lead) + cumsum(r$lengths); start <- end - r$lengths + 1
  cl  <- paste(sprintf("\\cline{%d-%d}", start, end)[r$lengths > 1], collapse = "")   # rule under each grouped cell
  paste(paste(c(rep("", length(lead)), top), collapse = " & "), "\\\\", cl,
        paste(c(sprintf("\\textbf{%s}", lead), sprintf("\\textbf{%s}", bot_lab)), collapse = " & "))
}
case_row <- function(sc, ncol) c(sprintf("\\multicolumn{%d}{|l|}{\\textit{%s}} \\\\", ncol, case_lab[[sc]]), "\\hline")
## rows = (instrument design, n) within each case; columns = the MR.RGM family quantities
write_table <- function(p_, cols, file, caption, label, size = "\\scriptsize") {
  hdr <- two_row_header(c("Instruments", "$n$"), names(cols))
  rows <- c()
  for (sc in ST_SCEN) { rows <- c(rows, case_row(sc, 2 + length(cols)))
    for (st in names(STRENGTH)) { for (n_ in ST_N) {
      rows <- c(rows, paste(c(if (n_ == ST_N[1]) ST_LABEL[[st]] else "", format(n_, big.mark = ","),
                              vapply(cols, function(cl) cell_val(sc, p_, st, n_, cl[1], cl[2], cl[3], if (length(cl) > 3) cl[4] else "mean"), character(1))), collapse = " & "), "\\\\") }
      rows <- c(rows, "\\hline") } }
  tex_table(file, sprintf("%s Entries are mean (SD) over %d replicates.", caption, max(summ$reps)), label,
            sprintf("|l|r|%s", paste(rep("c|", length(cols)), collapse = "")), hdr, rows, size)
}
## every method, one metric: one table per (p, case), rows = (instrument design, method), columns = n
write_method_table <- function(p_, sc, fam, met, stub, what) {
  tm  <- intersect(st_methods, unique(summ$method[summ$family == fam & summ$metric == met]))
  hdr <- paste(c("\\textbf{Instruments}", "\\textbf{Method}", sprintf("\\textbf{$n = %s$}", format(ST_N, big.mark = ","))), collapse = " & ")
  rows <- c()
  for (st in names(STRENGTH)) {
    r2 <- summ$mean[summ$scenario == sc & summ$p == p_ & summ$strength == st & summ$n == ST_N[1] & summ$metric == "R2_median"]
    for (k in seq_along(tm)) rows <- c(rows, paste(c(if (k == 1) sprintf("%s ($R^2 \\approx %.2f$)", ST_LABEL[[st]], if (length(r2)) r2 else NA) else "", method_tex[[tm[k]]],
      vapply(ST_N, function(n_) cell_val(sc, p_, st, n_, tm[k], fam, met), character(1))), collapse = " & "), "\\\\")
    rows <- c(rows, "\\hline") }
  tex_table(sprintf("strength_%s_%s_p=%d.tex", stub, scen_label[[sc]], p_),
            sprintf("Instrument-strength study, %s, small-world network, $p = %d$: %s. Entries are mean (SD) over %d replicates; entries marked $^\\dagger$ are over the successful fits only (Table~\\ref{tab:strength_failures} gives the number of failed fits); where fewer than %d fits succeeded the entry gives the number of failed fits instead.",
                    tolower(case_lab[[sc]]), p_, what, max(summ$reps), MIN_FITS),
            sprintf("tab:strength_%s_%s_p=%d", stub, scen_label[[sc]], p_), sprintf("|l|l|%s", paste(rep("c|", length(ST_N)), collapse = "")), hdr, rows)
}
scen_label <- c(smallworld = "nopleio", smallworld_pleio = "pleio")
for (p_ in ST_P) {
  for (sc in ST_SCEN) for (tb in list(list("graph", "AUC", "network_auc", "causal network recovery (AUC)"),
                                      list("effect", "MeanAbsDev", "effect_mad", "causal effect estimation (mean absolute deviation, MAD, over all off-diagonal entries; All-zero: the estimate 0)"),
                                      list("effect_edges", "MeanAbsDev", "effect_mad_edges", "causal effect estimation on the true edges (MAD)"),
                                      list("effect_nonedges", "MeanAbsDev", "effect_mad_nonedges", "mean absolute estimated effect over the true non-edges")))
    write_method_table(p_, sc, tb[[1]], tb[[2]], tb[[3]], tb[[4]])
  write_table(p_, list("Conf.\\ AUC MR.RGM" = c("MR.RGM", "confounding", "AUC"), "Conf.\\ AUC MR.RGM+" = c("MR.RGM+", "confounding", "AUC"),
                       "IV AUC MR.RGM+" = c("MR.RGM+", "iv_selection", "AUC"),
                       "Err.\\ ratio MR.RGM" = c("MR.RGM", "diagnostic", "err_ratio_max"), "Err.\\ ratio MR.RGM+" = c("MR.RGM+", "diagnostic", "err_ratio_max")),
              sprintf("strength_conf_iv_diag_p=%d.tex", p_),
              sprintf("MR.RGM and MR.RGM+ by instrument strength (small-world network, $p = %d$): confounding recovery (AUC), instrument selection by MR.RGM+ (AUC), and the largest ratio of estimated error variance to the variance of a trait (RGM() warns above 5).", p_),
              sprintf("tab:strength_conf_iv_diag_p=%d", p_))
}

## Convergence / weak-identification diagnostics table: chain-half PIP correlation (mean), fraction of fits that
## trigger the package's error-variance warning (max error-variance ratio > 5), and the spectral radius of the estimated A (mean)
diag_cell <- function(sc, p_, st, n_, meth, met, f = function(x) sprintf("%.2f", mean(x))) {
  x <- res$value[res$scenario == sc & res$p == p_ & res$strength == st & res$n == n_ & res$method == meth & res$family == "diagnostic" & res$metric == met]
  if (!length(x)) "--" else f(x)
}
for (p_ in ST_P) {
  hdr  <- two_row_header(c("Instruments", "$n$", "$R^2$", "$F$"),
                         c("Halves corr.\\ MR.RGM", "Halves corr.\\ MR.RGM+", "Err.-var.\\ warning (\\%) MR.RGM", "Err.-var.\\ warning (\\%) MR.RGM+",
                           "Spectral radius MR.RGM", "Spectral radius MR.RGM+"))
  rows <- c()
  for (sc in ST_SCEN) { rows <- c(rows, case_row(sc, 10)); for (st in names(STRENGTH)) { for (n_ in ST_N) {
    r2 <- summ$mean[summ$scenario == sc & summ$p == p_ & summ$strength == st & summ$n == n_ & summ$metric == "R2_median"]
    Fm <- summ$mean[summ$scenario == sc & summ$p == p_ & summ$strength == st & summ$n == n_ & summ$metric == "F_median"]
    pct <- function(x) sprintf("%.0f", 100 * mean(x))
    rows <- c(rows, paste(c(if (n_ == ST_N[1]) ST_LABEL[[st]] else "", format(n_, big.mark = ","),
                            if (length(r2)) sprintf("%.3f", r2) else "--", if (length(Fm)) format(round(Fm), big.mark = ",", scientific = FALSE) else "--",
                            diag_cell(sc, p_, st, n_, "MR.RGM", "pip_halves"), diag_cell(sc, p_, st, n_, "MR.RGM+", "pip_halves"),
                            diag_cell(sc, p_, st, n_, "MR.RGM", "warning_fired", pct), diag_cell(sc, p_, st, n_, "MR.RGM+", "warning_fired", pct),
                            diag_cell(sc, p_, st, n_, "MR.RGM", "spectral_radius"), diag_cell(sc, p_, st, n_, "MR.RGM+", "spectral_radius")),
                          collapse = " & "), "\\\\") }
    rows <- c(rows, "\\hline") } }
  tex_table(sprintf("strength_diag_p=%d.tex", p_),
            sprintf("Diagnostics of MR.RGM and MR.RGM+ by instrument strength (small-world network, $p = %d$): population $R^2$ of a trait on its own instruments and the implied $F$ statistic (medians over traits, mean over replicates); mean correlation of the causal posterior inclusion probabilities between the two halves of the chain; percentage of the %d replicates in which the largest error-variance ratio exceeds 5 (the package's error-variance warning); mean spectral radius of the estimated $A$.", p_, max(summ$reps)),
            sprintf("tab:strength_diag_p=%d", p_), "|l|r|c|c|c|c|c|c|c|c|", hdr, rows, "\\scriptsize")
}
## Warnings issued by RGM(): weak instruments (first-stage F < 10; computed when the map has exclusions, i.e. for
## MR.RGM), error variances, and any warning; percentage of the replicates
for (p_ in ST_P) {
  pct  <- function(x) sprintf("%.0f", 100 * mean(x))
  hdr  <- two_row_header(c("Instruments", "$n$"), c("Weak-IV (\\%) MR.RGM", "Err.-var.\\ (\\%) MR.RGM", "Any (\\%) MR.RGM",
                                                      "Err.-var.\\ (\\%) MR.RGM+", "Any (\\%) MR.RGM+"), by = "method")
  rows <- c()
  for (sc in ST_SCEN) { rows <- c(rows, case_row(sc, 7)); for (st in names(STRENGTH)) { for (n_ in ST_N) {
    rows <- c(rows, paste(c(if (n_ == ST_N[1]) ST_LABEL[[st]] else "", format(n_, big.mark = ","),
                            diag_cell(sc, p_, st, n_, "MR.RGM", "warn_weakIV", pct), diag_cell(sc, p_, st, n_, "MR.RGM", "warn_errvar", pct),
                            diag_cell(sc, p_, st, n_, "MR.RGM", "warn_any", pct),
                            diag_cell(sc, p_, st, n_, "MR.RGM+", "warn_errvar", pct), diag_cell(sc, p_, st, n_, "MR.RGM+", "warn_any", pct)),
                          collapse = " & "), "\\\\") }
    rows <- c(rows, "\\hline") } }
  tex_table(sprintf("strength_warn_p=%d.tex", p_),
            sprintf("Warnings issued by \\texttt{RGM()} by instrument strength (small-world network, $p = %d$): percentage of the %d replicates with the weak-instrument warning (first-stage $F < 10$; computed when the instrument--trait map has exclusions, that is, for \\texttt{MR.RGM}), the error-variance warning, and any warning.", p_, max(summ$reps)),
            sprintf("tab:strength_warn_p=%d", p_), "|l|r|c|c|c|c|c|", hdr, rows, "\\scriptsize")
}
## Failed fits per method and setting (every table and figure above summarises the successful fits only)
if (length(FAIL_COUNT)) {
  z  <- do.call(rbind, strsplit(names(FAIL_COUNT), "|", fixed = TRUE))       # scenario, p, n, strength, method
  fr <- data.frame(scenario = z[, 1], p = as.numeric(z[, 2]), n = as.numeric(z[, 3]), strength = z[, 4], method = z[, 5],
                   k = as.integer(FAIL_COUNT), stringsAsFactors = FALSE)
  fr <- fr[order(match(fr$method, method_order), fr$scenario, fr$p, match(fr$strength, names(STRENGTH)), fr$n), ]
  cause_lab <- c(significance = "no significant instrument", peeling = "peeling", numerical = "numerical", time = "time", error = "error")
  fr$cause <- vapply(seq_len(nrow(fr)), function(i) {             # e.g. "peeling 3, numerical 2"
    x <- fail_df$cause[st_key(fail_df) == paste(fr$scenario[i], fr$p[i], fr$n[i], fr$strength[i], fr$method[i], sep = "|")]
    tb <- table(factor(x, levels = names(cause_lab))); tb <- tb[tb > 0]
    paste(sprintf("%s %d", cause_lab[names(tb)], as.integer(tb)), collapse = ", ") }, character(1))
  body <- sprintf("%s & %s & %d & %s & %s & %d / %d & %d & %s \\\\", case_lab[fr$scenario], ST_LABEL[fr$strength], fr$p,
                  format(fr$n, big.mark = ","), method_tex[fr$method], fr$k, N_REP_ST, N_REP_ST - fr$k, fr$cause)
  tex_table("strength_failures.tex",
            "Instrument-strength study: fits without a result (number out of the replicates), the number of successful fits over which the method is summarized in every strength table and figure, and the cause (MR2G: no genome-wide significant instrument; GrIVET: peeling, its peeling step could not order every trait; numerical, it stopped with a numerical error; time, it did not finish within 48 hours).",
            "tab:strength_failures", "|l|l|r|r|l|c|c|l|",
            "\\textbf{Case} & \\textbf{Instruments} & \\textbf{$p$} & \\textbf{$n$} & \\textbf{Method} & \\textbf{Failed} & \\textbf{Successful} & \\textbf{Cause}",
            c(body, "\\hline"), "\\scriptsize")
}
writeLines(display_text(if (length(ST_FIG_NOTES)) ST_FIG_NOTES else "Every method has all replicates in every strength figure."),
           file.path(TAB_DIR, "strength_figure_notes.txt"))
message("Instrument-strength study: ", nrow(summ), " summary rows; figures/strength_*.pdf and tables/strength_*.tex written.")
}
