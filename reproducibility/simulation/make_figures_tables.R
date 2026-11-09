###############################################################################
## make_figures_tables.R
##
## Reads the result files of run_simulations.R (results_<experiment>/sim_*.rds) and writes
##   figures/  every simulation figure of the paper (file names as in the manuscript's \includegraphics)
##   tables/   every simulation table as a .tex file (labels and captions included), to be \input{}
##   tables/simulation_summary.csv          all means, SDs and replicate counts (every experiment)
##   tables/simulation_all_replicates.csv   every replicate value (long format)
##   figures/png/  300-dpi PNG copies of every figure (slides); the PDFs are vector graphics
##
## Experiments read (folders that do not exist are skipped):
##   results_main         Cases 1-3, p in {5, 10}, n in {500, ..., 30000}                      -> <case>_auc/_mad/_conf/_iv_p=*.pdf,
##                                                                                               main_summary_p=*.pdf, conf_recovery_p=*.pdf, feedback_*.tex
##   results_noise        error variance 1 / 4 / 9                                            -> noise_*.pdf, noise_*.tex
##   results_confounding  confounder loadings x 0 / 0.5 / 1 / 1.5                             -> confscale_*.pdf, confscale_*.tex
##   results_largen       n = 1e5, 1e6                                                        -> largen.pdf (with the main n), largen.tex, largen_seconds.tex
##   results_sensitivity  prior inclusion Beta(a,b), nu1, nu2                                 -> sensitivity_*.tex
##   results_roles        4 exposures -> 4 outcomes (MR2, MrDAG with their exposure/outcome split) -> roles_*.pdf, roles_*.tex
##   timing table over all experiments (mean seconds per fit)                                 -> timing.tex
##
## Every number is computed from the result files; every method present in a result file is reported.
## Usage (from this folder):  source("make_figures_tables.R")
## Optional before sourcing: EXPERIMENTS (subset of the folders above, without "results_")
###############################################################################

suppressPackageStartupMessages({ library(ggplot2); library(scales) })
source("sim_functions.R")                                       # effect_split_rows(): effect error by the truth, from fits_*.rds

## The whole script runs inside tryCatch, which reports which step failed and re-raises the error (run_all.R catches it
## and continues with the next step).
.mft_ok <- tryCatch({

if (!exists("EXPERIMENTS")) EXPERIMENTS <- c("main", "noise", "confounding", "largen", "sensitivity", "roles")
## RESULTS_PREFIX / FIG_DIR / TAB_DIR: defaults "results_", "figures", "tables"; the smoke test of this script uses
## RESULTS_PREFIX <- "results_smoke_"; FIG_DIR <- "figures_smoke"; TAB_DIR <- "tables_smoke"
if (!exists("RESULTS_PREFIX")) RESULTS_PREFIX <- "results_"
if (!exists("FIG_DIR")) FIG_DIR <- "figures"
if (!exists("TAB_DIR")) TAB_DIR <- "tables"
dir.create(FIG_DIR, showWarnings = FALSE); dir.create(TAB_DIR, showWarnings = FALSE)
dir.create(file.path(FIG_DIR, "png"), showWarnings = FALSE)

read_exp <- function(exp) {
  files <- list.files(paste0(RESULTS_PREFIX, exp), pattern = "^sim_.*\\.rds$", full.names = TRUE)
  if (!length(files)) return(NULL)
  strip <- function(r) { for (a in c("settings", "methods", "session_info", "package_versions")) attr(r, a) <- NULL; r }
  ## each result file plus the effect error on the true edges / non-edges and the All-zero reference, computed from
  ## the estimates saved next to it (fits_*.rds), in the same format
  d <- do.call(rbind, lapply(files, function(f) {
    r  <- as.data.frame(strip(readRDS(f))); ff <- sub("/sim_", "/fits_", f, fixed = TRUE)
    if (!file.exists(ff)) return(r)
    s  <- effect_split_rows(ff)
    k  <- r[rep(1, nrow(s)), c("experiment", "scenario", "p", "n", "setting", "setting_label"), drop = FALSE]
    s  <- cbind(k, s, seconds = NA_real_); rownames(s) <- NULL
    rbind(r, s[, names(r)]) }))
  list(data = d, n_rep = attr(readRDS(files[1]), "settings")$n_rep)
}
exp_data <- Filter(Negate(is.null), setNames(lapply(EXPERIMENTS, read_exp), EXPERIMENTS))
if (!length(exp_data)) stop("No result files found in results_<experiment>/; run run_simulations.R first.")
N_REP <- exp_data[[1]]$n_rep
res <- do.call(rbind, lapply(exp_data, `[[`, "data"))
rownames(res) <- NULL
message("Experiments read: ", paste(names(exp_data), collapse = ", "), " (", nrow(res), " rows)")
## GrIVET fits whose cluster job ended before they returned have no row in the result files; they are added as
## failed fits with cause "time" (add_unfinished() in sim_functions.R), in every experiment in which GrIVET was run,
## so that each GrIVET fit ends as a result, a peeling failure, a numerical failure or a fit stopped by the time limit.
g_exps <- unique(res$experiment[res$method == "GrIVET"])
if (length(g_exps)) {
  inx <- res$experiment %in% g_exps
  g_res <- add_unfinished(res[inx, , drop = FALSE], "GrIVET", c("experiment", "scenario", "p", "n", "setting_label"))
  UNFINISHED <- attr(g_res, "unfinished")
  res <- rbind(res[!inx, , drop = FALSE], g_res); rownames(res) <- NULL
}
if (!exists("UNFINISHED")) UNFINISHED <- NULL
## Cause of each failed GrIVET fit: "time" for the fits added above, otherwise the cause logged in the error logs
## (fail_causes() in sim_functions.R); a failed fit whose message is missing from the logs is reported as "error"
GRIVET_CAUSE <- fail_causes(file.path(paste0(RESULTS_PREFIX, names(exp_data)), "errors.log"), "GrIVET")
cause_key    <- function(d) paste(d$experiment, d$scenario, d$p, d$n, d$setting_label, d$rep, sep = "|")
fit_cause    <- function(d) {                                    # d: failed-fit rows
  gc <- setNames(GRIVET_CAUSE$cause, cause_key(GRIVET_CAUSE))
  uf <- if (is.null(UNFINISHED)) character(0) else cause_key(UNFINISHED)
  k  <- cause_key(d)
  ifelse(d$method == "GrIVET", ifelse(k %in% uf, "time", ifelse(k %in% names(gc), gc[k], "error")),
         ifelse(d$method == "MR2G", "significance", "error"))
}
local({ fr <- res[res$method == "GrIVET" & res$family == "diagnostic" & res$metric == "failed_fits" & !is.na(res$value) & res$value >= 1 & is.na(res$seconds), ]
        if (nrow(fr)) { cs <- fit_cause(fr); message("GrIVET failed fits by cause: ", paste(names(table(cs)), table(cs), sep = " ", collapse = ", "))
                        if (any(cs == "error")) message("  without a logged message: ", paste(cause_key(fr)[cs == "error"], collapse = "; ")) } })
## Failed fits (a method that stopped with an error, e.g. MR2G when a trait has no genome-wide significant
## instrument): the result files store them as "no network found" with seconds = NA and failed_fits = 1. Every
## table and figure summarises the SUCCESSFUL fits of a method only; replicate counts below N_REP are marked
## ($^\dagger$) and listed in the caption, and tables/figure_notes.txt gives the counts for every figure.
## (keys are built with paste(); interaction() would form the full cartesian product of the factor levels)
fit_key   <- function(d) paste(d$experiment, d$scenario, d$p, d$n, d$setting_label, d$method, d$rep, sep = "|")
fail_row  <- res$family == "diagnostic" & res$metric == "failed_fits" & is.na(res$seconds) & !is.na(res$value) & res$value >= 1
fail_keys <- unique(fit_key(res[fail_row, ]))
res$failed_fit <- fit_key(res) %in% fail_keys
write.csv(res, file.path(TAB_DIR, "simulation_all_replicates.csv"), row.names = FALSE)   # every row, failed fits flagged
res_all <- res
res <- res[!res$failed_fit | res$family == "diagnostic", ]
if (length(fail_keys)) message(length(fail_keys), " failed fits excluded from the summaries (see tables/failures.tex)")
## Number of failed fits per (experiment, scenario, p, n, setting, method). A method with fewer than MIN_FITS successful
## fits in a setting is not summarised there: its table entry reads "failed (k/N_REP)" and it is left out of the figures.
MIN_FITS    <- 10
set_key     <- function(d) paste(d$experiment, d$scenario, d$p, d$n, d$setting_label, d$method, sep = "|")
fd_         <- res_all[res_all$family == "diagnostic" & res_all$metric == "failed_fits" & res_all$failed_fit, ]
FAIL_COUNT  <- if (nrow(fd_)) table(set_key(fd_)) else table(character(0))
n_failed    <- function(key) if (key %in% names(FAIL_COUNT)) as.integer(FAIL_COUNT[[key]]) else 0L
TIME_COUNT  <- if (nrow(fd_)) table(set_key(fd_)[fit_cause(fd_) == "time"]) else table(character(0))
n_time      <- function(key) if (key %in% names(TIME_COUNT)) as.integer(TIME_COUNT[[key]]) else 0L
low_keys    <- names(FAIL_COUNT)[N_REP - as.integer(FAIL_COUNT) < MIN_FITS]
res <- res[!(set_key(res) %in% low_keys) | res$family == "diagnostic", ]

## ---------------------------------------------------------------------------
## Labels and order (identical in figures, tables and text)
## ---------------------------------------------------------------------------
method_order <- c("MR.RGM", "MR.RGM+", "MR.RGM_NoConf", "MR2G", "MVMR", "MR2", "MrDAG", "GrIVET",
                  "MR-SimpleMedian", "MR-WeightedMedian", "MR-IVW", "OneSampleMR", "mrbayes", "All-zero")
method_tex <- c("MR.RGM" = "\\texttt{MR.RGM}", "MR.RGM+" = "\\texttt{MR.RGM+}", "MR.RGM_NoConf" = "\\texttt{MR.RGM\\_NoConf}",
                "MR2G" = "\\texttt{MR2G}", "MVMR" = "\\texttt{MVMR}", "MR2" = "\\texttt{MR2}", "MrDAG" = "\\texttt{MrDAG}", "GrIVET" = "\\texttt{GrIVET}",
                "MR-SimpleMedian" = "\\texttt{MR-SimpleMedian}", "MR-WeightedMedian" = "\\texttt{MR-WeightedMedian}",
                "MR-IVW" = "\\texttt{MR-IVW}", "OneSampleMR" = "\\texttt{OneSampleMR}", "mrbayes" = "\\texttt{mrbayes}",
                "All-zero" = "\\textit{All-zero}")
## One colour per method, identical in every figure (Okabe-Ito, colour-blind safe, plus greys)
method_col <- c("MR.RGM" = "#0072B2", "MR.RGM+" = "#E69F00", "MR.RGM_NoConf" = "#56B4E9",
                "MR2G" = "#D55E00", "MVMR" = "#009E73", "MR2" = "#CC79A7", "MrDAG" = "#F0E442", "GrIVET" = "#882255",
                "MR-SimpleMedian" = "#999999", "MR-WeightedMedian" = "#666666", "MR-IVW" = "#BBBBBB",
                "OneSampleMR" = "#000000", "mrbayes" = "#8C510A", "All-zero" = "#4D4D4D")   # All-zero: reference (estimate 0), dark grey so that its dashed line is visible in print
## One plotting symbol and one line type per method, so that lines of similar colour remain distinguishable
method_shape <- c("MR.RGM" = 16, "MR.RGM+" = 17, "MR.RGM_NoConf" = 15, "MR2G" = 18, "MVMR" = 8, "MR2" = 7, "MrDAG" = 3,
                  "GrIVET" = 4, "MR-SimpleMedian" = 1, "MR-WeightedMedian" = 2, "MR-IVW" = 0, "OneSampleMR" = 5,
                  "mrbayes" = 6, "All-zero" = 13)
method_lty   <- c("MR.RGM" = "solid", "MR.RGM+" = "solid", "MR.RGM_NoConf" = "solid", "MR2G" = "longdash", "MVMR" = "longdash",
                  "MR2" = "longdash", "MrDAG" = "dotdash", "GrIVET" = "dotdash", "MR-SimpleMedian" = "dashed",
                  "MR-WeightedMedian" = "dashed", "MR-IVW" = "dashed", "OneSampleMR" = "dotted", "mrbayes" = "dotted",
                  "All-zero" = "twodash")
## The same palettes keyed by the displayed method names (display_text() in sim_functions.R), which are the factor
## labels built by meth_f(); used in every manual scale whose aesthetic maps such a factor
method_col_d   <- setNames(method_col,   display_text(names(method_col)))
method_shape_d <- setNames(method_shape, display_text(names(method_shape)))
method_lty_d   <- setNames(method_lty,   display_text(names(method_lty)))
## Colours and symbols of the sample sizes (orange: small n, blue: large n) in the figures with one column per method
n_col   <- c("#E69F00", "#D55E00", "#56B4E9", "#0072B2", "#009E73", "#000000")
n_shape <- c(16, 17, 15, 18, 8, 4)
scale_methods <- function(...) list(scale_colour_manual(values = method_col_d, ...), scale_shape_manual(values = method_shape_d, ...),
                                    scale_linetype_manual(values = method_lty_d, ...))
scen_text <- c(scalefree        = "a scale-free network with feedback loops and unmeasured confounding",
               smallworld       = "a small-world network with feedback loops and unmeasured confounding",
               smallworld_pleio = "a small-world network with feedback loops, unmeasured confounding, and horizontal pleiotropy",
               roles            = "the exposure--outcome design (4 exposures, 4 outcomes, unmeasured confounding)")
scen_label <- c(scalefree = "scalefree", smallworld = "smallworld", smallworld_pleio = "smallworld_horizontal", roles = "roles")
fig_stub   <- c(scalefree = "scalefree", smallworld = "smallworld", smallworld_pleio = "smallworld_pleio", roles = "roles")
## scen_title is displayed text only (tables, facet labels, titles): converted with display_text() here
scen_title <- display_text(c(scalefree = "scale-free", smallworld = "small-world", smallworld_pleio = "small-world, pleiotropy", roles = "exposures and outcomes"))
n_lab <- function(x) paste0("n = ", comma(x))
## factor of the method keys in the order of method_order, labelled with the displayed method names
meth_f <- function(m) { lv <- intersect(method_order, unique(m)); factor(m, levels = lv, labels = display_text(lv)) }
## Metrics reported in every summary table and figure, each for the methods that produce it:
##   network AUC and effect MAD: every method; confounding AUC: MR.RGM and MR.RGM+ (SSSL error covariance; in the roles
##   design also MR2 through its residual correlation); instrument-selection AUC: MR.RGM+ only (MR.RGM uses the known map).
METRICS <- list(
  auc  = list(metric = "AUC",        family = "graph",        label = "Network AUC",       what = "causal network recovery (AUC)"),
  mad  = list(metric = "MeanAbsDev", family = "effect",       label = "Effect MAD",        what = "causal effect estimation (mean absolute deviation over all off-diagonal entries; All-zero: the value of the estimate 0, i.e. of a method that reports no effect)"),
  mad_edges = list(metric = "MeanAbsDev", family = "effect_edges",    label = "Effect MAD, true edges", what = "causal effect estimation on the true edges (mean absolute deviation)"),
  mad_non   = list(metric = "MeanAbsDev", family = "effect_nonedges", label = "Mean |est.|, non-edges", what = "estimated effect where there is none (mean absolute estimate over the non-edges)"),
  conf = list(metric = "AUC",        family = "confounding",  label = "Confounding AUC",             what = "confounding structure recovery (AUC; MR.RGM and MR.RGM+)"),
  iv   = list(metric = "AUC",        family = "iv_selection", label = "Instrument selection AUC",  what = "instrument--trait selection by MR.RGM+ (AUC)"))
metric_methods <- function(s, mm) intersect(method_order, unique(s$method[s$family == mm$family & s$metric == mm$metric & !is.na(s$mean)]))

## ---------------------------------------------------------------------------
## Summary: mean, standard deviation across replicates, replicate count (per experiment and setting)
## ---------------------------------------------------------------------------
ok  <- res[!is.na(res$value), ]
key <- with(ok, paste(experiment, scenario, p, n, setting_label, method, family, metric, sep = "|"))
summ <- do.call(rbind, lapply(split(ok, key), function(d) {
  m <- nrow(d)
  data.frame(d[1, c("experiment", "scenario", "p", "n", "setting_label", "method", "family", "metric")],
             mean = mean(d$value), sd = if (m > 1) sd(d$value) else NA_real_, median = median(d$value), reps = m,
             seconds = if (all(is.na(d$seconds))) NA_real_ else mean(d$seconds, na.rm = TRUE))
}))
rownames(summ) <- NULL
write.csv(summ, file.path(TAB_DIR, "simulation_summary.csv"), row.names = FALSE)

fmt <- function(x) ifelse(is.na(x), "--",                       # one rule per column: three decimals, or "<0.001"
  ifelse(x > 0 & x < 0.0005, "$<0.001$", ifelse(x < 0 & x > -0.0005, "0.000",
  ifelse(x < 0, paste0("$-$", sprintf("%.3f", -x)), sprintf("%.3f", x)))))
cell <- function(s, exp, sc, p_, n_, meth, fam, met, lab = "base") {
  k <- n_failed(paste(exp, sc, p_, n_, lab, meth, sep = "|"))
  if (N_REP - k < MIN_FITS)                                                   # too few successful fits to summarise
    return(sprintf(if (n_time(paste(exp, sc, p_, n_, lab, meth, sep = "|")) == k) "not finished (%d/%d)" else "failed (%d/%d)", k, N_REP))
  r <- s[s$experiment == exp & s$scenario == sc & s$p == p_ & s$n == n_ & s$setting_label == lab & s$method == meth & s$family == fam & s$metric == met, ]
  if (!nrow(r)) return("--")
  paste0(fmt(r$mean), " (", fmt(r$sd), ")", if (r$reps < N_REP) "$^\\dagger$" else "")
}
## Cells over fewer than N_REP replicates (failed fits excluded): short_in() returns them, cap_tail() names the method,
## the reason and the number of successful fits in each such setting, e.g.
## "... MR2G stopped in some replicates (...); its entries marked $^\dagger$ are over its successful fits: 84 (scale-free, $p = 5$, $n = 500$)."
fail_reason <- c(MR2G = "it requires a genome-wide significant instrument ($p < 5\\times10^{-8}$, fixed in the package) for every trait",
                 GrIVET = "its peeling step could not order every trait, which a directed cycle can prevent; it stopped with a numerical error; or it did not finish within 48 hours")
setting_text <- function(sc, p_, n_, lab) paste0(if (sc %in% names(scen_title)) scen_title[[sc]] else sc, ", $p = ", p_, "$, $n = ", gsub(",", "{,}", format(n_, big.mark = ",", scientific = FALSE), fixed = TRUE), "$",
                                                 if (lab != "base") paste0(", ", gsub("_", " ", lab)) else "")
short_in <- function(s) s[s$reps < N_REP & s$method != "data" & s$family != "diagnostic", , drop = FALSE]
short_note <- function(sh, what = "entries marked $^\\dagger$ are") {
  if (!nrow(sh)) return("")
  ## the settings in which a method has fewer than N_REP replicates, with the smallest count over the metrics of the
  ## table or figure; a setting is attributed to failed fits when the method stopped with an error there (n_failed > 0),
  ## and otherwise to a metric that is undefined in some replicates (e.g. an AUC without both classes). Failed fits:
  ## one sentence per method; undefined metrics: one sentence per group of methods with the same counts.
  u <- aggregate(reps ~ experiment + method + scenario + p + n + setting_label, data = sh, FUN = min)
  u$failed <- vapply(seq_len(nrow(u)), function(i) n_failed(paste(u$experiment[i], u$scenario[i], u$p[i], u$n[i], u$setting_label[i], u$method[i], sep = "|")) > 0, logical(1))
  u <- u[order(match(u$method, method_order), u$scenario, u$p, u$n), ]
  txt <- function(x) paste(sprintf("%d (%s)", x$reps, mapply(setting_text, x$scenario, x$p, x$n, x$setting_label)), collapse = "; ")
  mtex <- function(m) if (m %in% names(method_tex)) method_tex[[m]] else m
  fail_txt <- vapply(split(u[u$failed, , drop = FALSE], factor(u$method[u$failed], levels = unique(u$method[u$failed]))), function(xf) {
    m <- xf$method[1]
    low_m <- any(endsWith(low_keys, paste0("|", m)) & sub("\\|.*", "", low_keys) %in% xf$experiment)
    paste0(mtex(m), " did not return a result in some replicates", if (m %in% names(fail_reason)) paste0(" (", fail_reason[[m]], ")") else "",
           "; its ", what, " over its successful fits only (", if (min(xf$reps) == max(xf$reps)) min(xf$reps) else paste0(min(xf$reps), "--", max(xf$reps)), " of ", N_REP, "; counts in Table~\\ref{tab:failures}).",
           if (low_m) sprintf(" Where fewer than %d fits succeeded, the entry gives the number of failed fits instead.", MIN_FITS) else "")
  }, character(1))
  uu <- u[!u$failed, , drop = FALSE]
  und_txt <- character(0)
  if (nrow(uu)) {
    per_m <- vapply(split(uu, factor(uu$method, levels = unique(uu$method))), txt, character(1))
    for (t in unique(per_m)) {
      ms <- names(per_m)[per_m == t]
      und_txt <- c(und_txt, paste0("For ", paste(vapply(ms, mtex, ""), collapse = ", "), " a metric is undefined in some replicates (the truth has a single class there, so the AUC is not defined); ",
                                   if (length(ms) > 1) "their " else "its ", what, " over the remaining replicates, out of ", N_REP,
                                   " (smallest count over the metrics shown): ", t, "."))
    }
  }
  paste(c(fail_txt, und_txt), collapse = " ")
}
cap_tail <- function(sh) paste0(" Entries are mean (SD) over ", N_REP, " replicates.", if (nrow(sh)) paste0(" ", short_note(sh)) else "")
## Every table is written at a readable font size (\\footnotesize, or \\scriptsize for tables with many numeric columns);
## adjustbox only guards against a table that is still wider than the text (needs \\usepackage{adjustbox}).
table_tex <- function(file, caption, label, colspec, header, body, size = "\\footnotesize") {
  L <- c("\\begin{table}[htbp]", "\\centering", size, "\\setlength{\\tabcolsep}{3.5pt}", sprintf("\\caption{%s}", caption), sprintf("\\label{%s}", label),
         "\\begin{adjustbox}{max width=\\textwidth, max totalheight=0.78\\textheight}",
         sprintf("\\begin{tabular}{%s}", colspec), "\\hline", paste(header, "\\\\"), "\\hline", "\\hline", body, "\\end{tabular}",
         "\\end{adjustbox}", "\\end{table}")
  L <- gsub("(?<=[0-9]),(?=[0-9]{3}(?![0-9]))", "{,}", L, perl = TRUE)   # thousands separator that also renders correctly in math mode
  writeLines(display_text(L), file.path(TAB_DIR, file))          # displayed method names and wording (sim_functions.R)
}
## Long table: one row per (setting, method), one column per sample size (or level) for one metric; the readable
## layout for a metric reported for every method. rows = data frame with columns group (first column, repeated
## values printed once), method; cols = named list of column keys; cellfun(row_index, col_key) returns the entry.
long_table_tex <- function(file, caption, label, rows, group_label, cols, cellfun, size = "\\footnotesize") {
  hdr  <- paste(c(sprintf("\\textbf{%s}", group_label), "\\textbf{Method}", sprintf("\\textbf{%s}", names(cols))), collapse = " & ")
  body <- c()
  for (i in seq_len(nrow(rows))) {
    first <- i == 1 || rows$group[i] != rows$group[i - 1]
    if (first && i > 1) body <- c(body, "\\hline")
    body <- c(body, paste(c(if (first) rows$group[i] else "", method_tex[[rows$method[i]]],
                            vapply(cols, function(cl) cellfun(i, cl), character(1))), collapse = " & "), "\\\\")
  }
  table_tex(file, caption, label, sprintf("|l|l|%s", paste(rep("c|", length(cols)), collapse = "")), hdr, body, size)
}
## Figure saving: paper size (large fonts, scaled down in LaTeX) unless a figure driver has defined THESIS_STYLE for
## print-size figures: then theme_bw(base_size = 9), sizes from the driver's CSV, cairo_pdf, no png.
STYLE  <- if (exists("THESIS_STYLE")) THESIS_STYLE else NULL
THESIS <- !is.null(STYLE)
save_fig <- function(g, file, width, height) {
  if (THESIS) return(STYLE$save(g, file.path(FIG_DIR, file), width, height))
  ggsave(file.path(FIG_DIR, file), g, width = width, height = height, units = "in")
  ggsave(file.path(FIG_DIR, "png", sub("\\.pdf$", ".png", file)), g, width = width, height = height, units = "in", dpi = 300)
}
## Caption notes for figures with fewer than N_REP replicates for some method (written to tables/figure_notes.txt):
## d = per-replicate rows (column rep) or summary rows (from_summary = TRUE, column reps)
FIG_NOTES <- character()
fig_note <- function(file, d, from_summary = FALSE) {
  if (is.null(d) || !nrow(d)) return(invisible(NULL))
  if (!from_summary) {
    d <- aggregate(rep ~ experiment + method + scenario + p + n + setting_label, data = d, FUN = function(x) length(unique(x)))
    names(d)[names(d) == "rep"] <- "reps"
  }
  sh <- d[d$reps < N_REP & d$method != "data", , drop = FALSE]
  if (nrow(sh)) FIG_NOTES <<- c(FIG_NOTES, paste0(file, ": ", short_note(sh, "results (boxes, means) are")))
  ## methods left out of this figure in some setting (fewer than MIN_FITS successful fits)
  ctx  <- unique(paste(d$experiment, d$scenario, d$p, d$n, d$setting_label, sep = "|"))
  gone <- low_keys[sub("\\|[^|]*$", "", low_keys) %in% ctx]
  if (length(gone)) FIG_NOTES <<- c(FIG_NOTES, paste0(file, ": not shown (fewer than ", MIN_FITS, " successful fits): ",
    paste(vapply(strsplit(gone, "|", fixed = TRUE), function(z) { key <- paste(z, collapse = "|")
      sprintf(if (n_time(key) == n_failed(key)) "%s, did not finish within 48 hours in %d of %d (%s)" else "%s, failed in %d of %d (%s)", z[6],
              n_failed(key), N_REP, setting_text(z[2], z[3], as.numeric(z[4]), z[5])) }, character(1)), collapse = "; "), "."))
  invisible(NULL)
}
theme_box <- theme_bw(base_size = 16) +
  theme(axis.text.x  = element_text(size = 24, angle = 45, hjust = 1),
        axis.text.y  = element_text(size = 20),
        axis.title.x = element_text(size = 26, face = "bold"),
        axis.title.y = element_text(size = 26, face = "bold"),
        strip.text   = element_text(size = 22, face = "bold"),
        plot.title   = element_text(size = 30, face = "bold", hjust = 0.5))
theme_line <- theme_bw(base_size = 16) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        plot.title   = element_text(size = 24, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 24, face = "bold", margin = margin(t = 8)),
        axis.title.y = element_text(size = 24, face = "bold", margin = margin(r = 8)),
        axis.text.x  = element_text(size = 24), axis.text.y = element_text(size = 24),
        strip.text   = element_text(size = 22, face = "bold"),
        legend.text  = element_text(size = 20), legend.title = element_blank(), legend.position = "top")
if (THESIS) {                                   # same layout, fonts from theme_bw(base_size = 9)
  theme_box <- theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold", hjust = 0.5))
  theme_line <- theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          legend.title = element_blank(), legend.position = "top", legend.margin = margin(0, 0, 0, 0))
}

## ===========================================================================
## 1. Main experiment: tables and figures per case and p (file names of the manuscript)
## ===========================================================================
if ("main" %in% names(exp_data)) {
  sm <- summ[summ$experiment == "main", ]; rm_ <- res[res$experiment == "main", ]
  n_levels <- sort(unique(rm_$n))

  write_method_table <- function(sc, p_, fam) {
    ## columns: metric (family) pairs; the effect tables also give the error on the true edges and the mean absolute
    ## estimate on the non-edges, so that the per-case tables are the complete record
    cols <- if (fam == "graph") list("AUC" = c("graph", "AUC"), "TPR" = c("graph", "TPR"), "FDR" = c("graph", "FDR"), "MCC" = c("graph", "MCC"))
            else list("Max.\\ abs.\\ dev." = c("effect", "MaxAbsDev"), "MAD" = c("effect", "MeanAbsDev"), "Mean sq.\\ dev." = c("effect", "MeanSqDev"),
                      "MAD, true edges" = c("effect_edges", "MeanAbsDev"), "Mean $|$est.$|$, non-edges" = c("effect_nonedges", "MeanAbsDev"))
    meths <- intersect(method_order, unique(sm$method[sm$scenario == sc & sm$family == fam]))
    what  <- if (fam == "graph") "Graph recovery" else "Causal effect estimation (mean absolute deviation, MAD, over all off-diagonal entries, on the true edges, and the mean absolute estimate on the true non-edges; All-zero: the estimate 0)"
    lab   <- sprintf("tab:feedback_%s_%s_90_p=%d", if (fam == "graph") "graph_structure" else "causal_effect", scen_label[[sc]], p_)
    body <- c()
    for (n_ in n_levels) {
      for (k in seq_along(meths)) {
        vals <- vapply(cols, function(cl) cell(sm, "main", sc, p_, n_, meths[k], cl[1], cl[2]), character(1))
        body <- c(body, paste0(if (k == 1) sprintf("$n = %s$", format(n_, big.mark = ",")) else "", " & ",
                               method_tex[[meths[k]]], " & ", paste(vals, collapse = " & "), " \\\\"))
      }
      body <- c(body, "\\hline")
    }
    fams <- unique(vapply(cols, `[`, "", 1))
    table_tex(paste0(sub("^tab:", "", lab), ".tex"),
              sprintf("%s in %s, with network size \\(p{=}%d\\).%s", what, scen_text[[sc]], p_,
                      cap_tail(short_in(sm[sm$scenario == sc & sm$family %in% fams & sm$p == p_, ]))),
              lab, sprintf("|c|c|%s", paste(rep("c|", length(cols)), collapse = "")),
              paste0("\\textbf{Setting} & \\textbf{Method} & ", paste0("\\textbf{", names(cols), "}", collapse = " & ")), body,
              if (fam == "graph") "\\footnotesize" else "\\scriptsize")
  }
  write_conf_table <- function(sc, meth, suffix) {
    mets <- c("AUC", "TPR", "FDR", "MCC")
    lab  <- sprintf("tab:feedback_conf_auc_%s_%s", scen_label[[sc]], suffix)
    body <- c()
    for (p_ in sort(unique(sm$p[sm$scenario == sc]))) {
      body <- c(body, sprintf("\\multirow{%d}{*}{$p = %d$}", length(n_levels), p_))
      for (n_ in n_levels) {
        vals <- vapply(mets, function(mt) cell(sm, "main", sc, p_, n_, meth, "confounding", mt), character(1))
        body <- c(body, paste0("    & ", format(n_, big.mark = ","), " & ", paste(vals, collapse = " & "), " \\\\"))
      }
      body <- c(body, "\\hline")
    }
    table_tex(paste0(sub("^tab:", "", lab), ".tex"),
              sprintf("Confounding structure recovery performance using %s in %s, across network sizes \\(p\\in\\{5,10\\}\\).%s",
                      method_tex[[meth]], scen_text[[sc]], cap_tail(short_in(sm[sm$scenario == sc & sm$family == "confounding" & sm$method == meth, ]))),
              lab, "|c|c|c|c|c|c|",
              "\\textbf{Setting} & \\textbf{Sample size} & \\textbf{AUC} & \\textbf{TPR} & \\textbf{FDR} & \\textbf{MCC}", body, "\\small")
  }
  for (sc in unique(rm_$scenario)) {
    for (p_ in sort(unique(rm_$p[rm_$scenario == sc]))) for (fam in c("graph", "effect")) write_method_table(sc, p_, fam)
    write_conf_table(sc, "MR.RGM", "A")
    if (any(sm$method == "MR.RGM+" & sm$scenario == sc)) write_conf_table(sc, "MR.RGM+", "B")
  }

  ## Compact overview: one table per metric (network AUC, effect MAD, confounding AUC, instrument-selection AUC),
  ## every method that produces the metric, all cases, p = 5 and 10, at each n
  for (key in c("conf", "iv")) {                    # the other metrics are in the per-case tables (one row per method)
    mm <- METRICS[[key]]; mt <- mm$metric; fam <- mm$family
    meths <- metric_methods(sm, mm); if (!length(meths)) next
    hdr <- paste(c("\\textbf{Case}", "\\textbf{$p$}", "\\textbf{$n$}", sprintf("\\textbf{%s}", method_tex[meths])), collapse = " & ")
    body <- c()
    for (sc in intersect(names(scen_text), unique(sm$scenario))) for (p_ in sort(unique(sm$p[sm$scenario == sc]))) {
      for (n_ in n_levels) body <- c(body, paste(c(if (n_ == n_levels[1]) sprintf("\\multirow{%d}{*}{%s}", length(n_levels), scen_title[[sc]]) else "",
                                                   if (n_ == n_levels[1]) p_ else "", format(n_, big.mark = ","),
                                                   vapply(meths, function(m) cell(sm, "main", sc, p_, n_, m, fam, mt), character(1))), collapse = " & "), "\\\\")
      body <- c(body, "\\hline")
    }
    table_tex(sprintf("overview_%s.tex", key),
              sprintf("Main simulation study, three cases: %s.%s", mm$what, cap_tail(short_in(sm[sm$family == fam, ]))),
              sprintf("tab:overview_%s", key), sprintf("|l|c|r|%s", paste(rep("c|", length(meths)), collapse = "")), hdr, body)
  }

  method_box <- function(sc, p_, fam, met, ylab, title, file) {
    d <- rm_[rm_$scenario == sc & rm_$p == p_ & rm_$family == fam & rm_$metric == met & !is.na(rm_$value), ]
    fig_note(file, d)
    d$Method <- meth_f(d$method)
    d$n <- factor(d$n, levels = n_levels, labels = n_lab(n_levels))
    g <- ggplot(d, aes(Method, value, fill = Method)) + geom_boxplot(colour = "grey15", outlier.size = 1.8) +
      scale_fill_manual(values = method_col_d) +
      facet_wrap(~ n, scales = if (fam == "graph") "free_x" else "free") + guides(fill = "none") + theme_box +   # effect errors: own y range per n
      labs(title = display_text(title), y = ylab, x = "Method")
    save_fig(g, file, 15, 12)
  }
  single_box <- function(d, title, file, ylim_low) {
    fig_note(file, d)
    d$n      <- factor(d$n, levels = n_levels, labels = n_lab(n_levels))
    d$Method <- factor(d$method, levels = c("MR.RGM", "MR.RGM+"), labels = display_text(c("MR.RGM", "MR.RGM+")))
    two <- length(unique(d$method)) > 1
    lo  <- min(ylim_low, min(d$value, na.rm = TRUE))
    g <- ggplot(d, aes(n, value)) + geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.4, colour = "grey60")
    if (two) {
      g <- g + geom_boxplot(aes(fill = Method), width = 0.6, colour = "grey15", alpha = 0.88,
                            outlier.shape = 16, outlier.size = 1.8, position = position_dodge(width = 0.75)) +
        scale_fill_manual(values = method_col_d[display_text(c("MR.RGM", "MR.RGM+"))])
    } else {
      g <- g + geom_boxplot(width = 0.6, fill = method_col[[d$method[1]]], colour = "grey15", alpha = 0.88,
                            outlier.shape = 16, outlier.size = 1.8) +
        geom_jitter(width = 0.08, height = 0, size = 1.4, alpha = 0.5, colour = "grey20")
    }
    g <- g + coord_cartesian(ylim = c(lo, 1)) + labs(title = display_text(title), x = "Sample size", y = "Area Under Curve (AUC)") + theme_line
    save_fig(g, file, 14, 8)
  }
  for (sc in unique(rm_$scenario)) for (p_ in sort(unique(rm_$p[rm_$scenario == sc]))) {
    st <- fig_stub[[sc]]
    method_box(sc, p_, "graph",  "AUC",        "Area Under Curve (AUC)",        sprintf("Causal network recovery (%s, p = %d)", scen_title[[sc]], p_), sprintf("%s_auc_p=%d.pdf", st, p_))
    method_box(sc, p_, "effect", "MeanAbsDev", "Mean Absolute Deviation (MAD)", sprintf("Causal effect estimation (%s, p = %d)", scen_title[[sc]], p_), sprintf("%s_mad_p=%d.pdf", st, p_))
    dc <- rm_[rm_$scenario == sc & rm_$p == p_ & rm_$family == "confounding" & rm_$metric == "AUC" & rm_$method %in% c("MR.RGM", "MR.RGM+") & !is.na(rm_$value), ]
    if (nrow(dc)) single_box(dc, sprintf("Confounding recovery (%s, p = %d)", scen_title[[sc]], p_), sprintf("%s_conf_p=%d.pdf", st, p_), 0.4)
    di <- rm_[rm_$scenario == sc & rm_$p == p_ & rm_$family == "iv_selection" & !is.na(rm_$value), ]
    if (nrow(di)) single_box(di, sprintf("Instrument selection by MR.RGM+ (%s, p = %d)", scen_title[[sc]], p_),
                             sprintf("%s_iv_p=%d.pdf", st, p_), 0.85)
  }

  ## Summary figure of the main experiment, one per p: mean (+- SD) of the network AUC and the effect MAD of every
  ## method, one column of points per method (so that no two methods overlap), the sample sizes side by side within
  ## a method (colour and symbol), one panel column per case; the MR.RGM family is shaded. All-zero (the MAD of the
  ## estimate 0) is drawn as a dashed line in the MAD panels.
  for (p_ in sort(unique(sm$p))) {
    dd <- do.call(rbind, lapply(METRICS[c("auc", "mad")], function(mm) {
      x <- sm[sm$p == p_ & sm$family == mm$family & sm$metric == mm$metric & !is.na(sm$mean), ]
      if (!nrow(x)) return(NULL); x$Metric <- mm$label; x }))
    if (is.null(dd) || !nrow(dd)) next
    fig_note(sprintf("main_summary_p=%d.pdf", p_), dd, from_summary = TRUE)
    dd$Metric <- factor(dd$Metric, levels = c("Network AUC", "Effect MAD"))
    dd$Case   <- factor(scen_title[dd$scenario], levels = scen_title[intersect(names(scen_title), unique(dd$scenario))])
    zero <- if (any(dd$method == "All-zero")) aggregate(mean ~ Case + Metric, data = dd[dd$method == "All-zero", ], FUN = mean) else NULL
    dd <- dd[dd$method != "All-zero", ]
    dd$Method <- meth_f(dd$method)
    dd$N      <- factor(dd$n, levels = n_levels, labels = n_lab(n_levels))
    n_fam <- sum(levels(dd$Method) %in% display_text(c("MR.RGM", "MR.RGM+", "MR.RGM_NoConf")))   # levels: displayed names
    g <- ggplot(dd, aes(Method, mean, colour = N, shape = N)) + scale_x_discrete() +     # discrete axis before the shading
      annotate("rect", xmin = 0.5, xmax = n_fam + 0.5, ymin = -Inf, ymax = Inf, fill = "grey92")
    if (!is.null(zero)) g <- g + geom_hline(data = zero, aes(yintercept = mean), inherit.aes = FALSE,
                                            linetype = "dashed", colour = "grey40", linewidth = 0.4)
    g <- g + geom_pointrange(aes(ymin = mean - sd, ymax = mean + sd), position = position_dodge(width = 0.75),
                      size = if (THESIS) 0.18 else 0.4, linewidth = if (THESIS) 0.35 else 0.6) +
      scale_colour_manual(values = setNames(n_col[seq_along(n_levels)], n_lab(n_levels)), name = NULL) +
      scale_shape_manual(values = setNames(n_shape[seq_along(n_levels)], n_lab(n_levels)), name = NULL) +
      facet_grid(Metric ~ Case, scales = "free_y") +
      labs(x = NULL, y = NULL, title = sprintf("Main simulation study (p = %d): mean and SD over replicates", p_)) +
      theme_line + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.x = element_blank())
    save_fig(g, sprintf("main_summary_p=%d.pdf", p_), 16, 10)
  }

  ## Confounding recovery in one figure per p: boxes against n, one column per case (the per-case figures are kept)
  for (p_ in sort(unique(rm_$p))) {
    dc <- rm_[rm_$p == p_ & rm_$family == "confounding" & rm_$metric == "AUC" & rm_$method %in% c("MR.RGM", "MR.RGM+") & !is.na(rm_$value), ]
    if (!nrow(dc)) next
    fig_note(sprintf("conf_recovery_p=%d.pdf", p_), dc)
    dc$n <- factor(dc$n, levels = n_levels, labels = n_lab(n_levels))
    dc$Method <- factor(dc$method, levels = c("MR.RGM", "MR.RGM+"), labels = display_text(c("MR.RGM", "MR.RGM+")))
    dc$Case <- factor(scen_title[dc$scenario], levels = scen_title[intersect(names(scen_title), unique(dc$scenario))])
    g <- ggplot(dc, aes(n, value)) + geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.4, colour = "grey60") +
      geom_boxplot(aes(fill = Method), width = 0.65, colour = "grey15", alpha = 0.88, outlier.shape = 16, outlier.size = 1.2,
                   position = position_dodge(width = 0.8)) +
      scale_fill_manual(values = method_col_d[display_text(c("MR.RGM", "MR.RGM+"))]) + facet_wrap(~ Case, nrow = 1) +
      labs(title = sprintf("Confounding structure recovery (p = %d)", p_), x = "Sample size", y = "Area Under Curve (AUC)") +
      theme_line + theme(axis.text.x = element_text(angle = 30, hjust = 1))
    save_fig(g, sprintf("conf_recovery_p=%d.pdf", p_), 16, 6)
  }
}

## ===========================================================================
## 2. One-factor experiments (noise, confounding): table + line figure of AUC and MAD against the factor
## ===========================================================================
factor_experiment <- function(exp, factor_name, factor_tex, factor_from_label, file_stub, caption_what) {
  if (!exp %in% names(exp_data)) return(invisible(NULL))
  se <- summ[summ$experiment == exp, ]; re <- res[res$experiment == exp, ]
  se$level <- factor_from_label(se$setting_label); re$level <- factor_from_label(re$setting_label)
  levels_ <- sort(unique(se$level)); ns <- sort(unique(se$n))
  ## tables: every method, one row per (level, method); columns = metric x n (two tables: network recovery and effect
  ## estimation; effect estimation on the true edges and non-edges), then the narrow confounding and
  ## instrument-selection tables (rows = level x n)
  lab_of <- function(lv) unique(se$setting_label[se$level == lv])[1]
  for (tb in list(list(key = "auc", cols = c("auc", "mad"), what = "causal network recovery (AUC) and causal effect estimation (MAD over all off-diagonal entries; All-zero: the estimate 0)"),
                  list(key = "mad", cols = c("mad_edges", "mad_non"), what = "causal effect estimation on the true edges (MAD) and the mean absolute estimate on the true non-edges"))) {
    meths <- unique(unlist(lapply(METRICS[tb$cols], function(mm) metric_methods(se, mm)))); meths <- intersect(method_order, meths)
    if (!length(meths)) next
    rows <- expand.grid(method = meths, level = levels_, stringsAsFactors = FALSE)[, c("level", "method")]
    rows$group <- as.character(rows$level)
    tex_lab <- c(auc = "Network AUC", mad = "Effect MAD", mad_edges = "MAD, true edges", mad_non = "Mean $|$est.$|$, non-edges")
    cols <- list(); for (ck in tb$cols) for (n_ in ns)                 # two-line column headers: metric, then n
      cols[[sprintf("\\shortstack{%s\\\\ $n = %s$}", tex_lab[[ck]], format(n_, big.mark = "{,}"))]] <- c(ck, n_)
    fams <- vapply(METRICS[tb$cols], `[[`, "", "family")
    long_table_tex(sprintf("%s_%s.tex", file_stub, tb$key),
                   sprintf("%s: %s by %s (small-world network with feedback loops and unmeasured confounding, $p = 10$).%s",
                           caption_what, tb$what, tolower(factor_name), cap_tail(short_in(se[se$family %in% fams, ]))),
                   sprintf("tab:%s_%s", file_stub, tb$key), rows, factor_tex, cols,
                   function(i, cl) { mm <- METRICS[[cl[1]]]; cell(se, exp, "smallworld", 10, as.numeric(cl[2]), rows$method[i], mm$family, mm$metric, lab_of(rows$level[i])) },
                   "\\scriptsize")
  }
  for (key in c("conf", "iv")) {
    mm <- METRICS[[key]]; mt <- mm$metric; fam <- mm$family
    meths <- metric_methods(se, mm); if (!length(meths)) next
    hdr <- paste(c(sprintf("\\textbf{%s}", factor_tex), "\\textbf{$n$}", sprintf("\\textbf{%s}", method_tex[meths])), collapse = " & ")
    body <- c()
    for (lv in levels_) { for (n_ in ns) {
      body <- c(body, paste(c(if (n_ == ns[1]) format(lv) else "", format(n_, big.mark = ","),
                              vapply(meths, function(m) cell(se, exp, "smallworld", 10, n_, m, fam, mt, lab_of(lv)), character(1))),
                            collapse = " & "), "\\\\") }
      body <- c(body, "\\hline") }
    table_tex(sprintf("%s_%s.tex", file_stub, key),
              sprintf("%s: %s by %s (small-world network with feedback loops and unmeasured confounding, $p = 10$).%s",
                      caption_what, mm$what, tolower(factor_name), cap_tail(short_in(se[se$family == fam, ]))),
              sprintf("tab:%s_%s", file_stub, key), sprintf("|c|r|%s", paste(rep("c|", length(meths)), collapse = "")), hdr, body)
  }
  ## figure: mean +- SD against the factor, one line per method, rows = the four metrics, panels = n
  ## the instrument-selection AUC of MR.RGM+ is at or near 1 in every setting of these experiments (tables); it is not drawn
  d <- do.call(rbind, lapply(METRICS[c("auc", "mad", "conf")], function(mm) { x <- se[se$family == mm$family & se$metric == mm$metric & !is.na(se$mean), ]
    if (!nrow(x)) return(NULL); x$Metric <- mm$label; x }))
  d$Metric <- factor(d$Metric, levels = vapply(METRICS, `[[`, "", "label"))
  fig_note(sprintf("%s.pdf", file_stub), d, from_summary = TRUE)
  d$Method <- meth_f(d$method); d$N <- factor(d$n, levels = ns, labels = n_lab(ns))
  d <- droplevels(d)
  g <- ggplot(d, aes(level, mean, colour = Method, shape = Method, linetype = Method, group = Method)) +
    geom_line(linewidth = if (THESIS) 0.45 else 0.8, position = position_dodge(width = 0.03 * diff(range(levels_)))) +
    geom_pointrange(aes(ymin = mean - sd, ymax = mean + sd), size = if (THESIS) 0.25 else 0.4, linewidth = if (THESIS) 0.3 else 0.5,
                    position = position_dodge(width = 0.03 * diff(range(levels_)))) +
    scale_methods() + facet_grid(Metric ~ N, scales = "free_y") +
    scale_x_continuous(breaks = levels_) + labs(x = factor_name, y = NULL, title = display_text(caption_what)) + theme_line +
    theme(legend.key.width = unit(1.6, "lines"))
  save_fig(g, sprintf("%s.pdf", file_stub), 14, 4 + 2.5 * nlevels(d$Metric))
}
factor_experiment("noise", "Error variance", "$\\sigma^2$", function(l) as.numeric(sub("^errvar_", "", l)), "noise", "Error variance")
factor_experiment("confounding", "Confounder loading scale", "Scale", function(l) as.numeric(sub("^conf_", "", l)), "confscale", "Strength of confounding")

## ---------------------------------------------------------------------------
## Selection metrics of the noise and confounding experiments: TPR and FDR of the selected causal graph (every method)
## and of the selected confounding links (MR.RGM family), mean (SD) over the replicates
## ---------------------------------------------------------------------------
selection_tables <- function(exp, factor_name, factor_tex, factor_from_label, file_stub, caption_what) {
  if (!exp %in% names(exp_data)) return(invisible(NULL))
  se <- summ[summ$experiment == exp, ]
  se$level <- factor_from_label(se$setting_label); levels_ <- sort(unique(se$level)); ns <- sort(unique(se$n))
  lab_of <- function(lv) unique(se$setting_label[se$level == lv])[1]
  meths <- intersect(method_order, unique(se$method[se$family == "graph" & se$metric == "TPR" & !is.na(se$mean)]))
  if (length(meths)) {
    rows <- expand.grid(method = meths, level = levels_, stringsAsFactors = FALSE)[, c("level", "method")]
    rows$group <- as.character(rows$level)
    cols <- list(); for (mt in c("TPR", "FDR")) for (n_ in ns)
      cols[[sprintf("\\shortstack{%s\\\\ $n = %s$}", mt, format(n_, big.mark = "{,}"))]] <- c(mt, n_)
    long_table_tex(sprintf("%s_sel.tex", file_stub),
                   sprintf("%s: true positive rate (TPR) and false discovery rate (FDR) of the selected causal graph (PIP $> 0.5$ for the MR.RGM family, MR2 and MrDAG; $p < 0.05$, or a posterior tail probability below 0.05 for mrbayes, for the others) by %s (small-world network with feedback loops and unmeasured confounding, $p = 10$).%s",
                           caption_what, tolower(factor_name), cap_tail(short_in(se[se$family == "graph", ]))),
                   sprintf("tab:%s_sel", file_stub), rows, factor_tex, cols,
                   function(i, cl) cell(se, exp, "smallworld", 10, as.numeric(cl[2]), rows$method[i], "graph", cl[1], lab_of(rows$level[i])),
                   "\\scriptsize")
  }
  cm <- intersect(method_order, unique(se$method[se$family == "confounding" & se$metric == "TPR" & !is.na(se$mean)]))
  if (length(cm)) {
    hdr  <- paste(c(sprintf("\\textbf{%s}", factor_tex), "\\textbf{$n$}",
                    unlist(lapply(cm, function(m) sprintf("\\textbf{%s %s}", c("TPR", "FDR"), method_tex[[m]])))), collapse = " & ")
    body <- c()
    for (lv in levels_) { for (n_ in ns) {
      body <- c(body, paste(c(if (n_ == ns[1]) format(lv) else "", format(n_, big.mark = ","),
                              unlist(lapply(cm, function(m) c(cell(se, exp, "smallworld", 10, n_, m, "confounding", "TPR", lab_of(lv)),
                                                               cell(se, exp, "smallworld", 10, n_, m, "confounding", "FDR", lab_of(lv)))))),
                            collapse = " & "), "\\\\") }
      body <- c(body, "\\hline") }
    table_tex(sprintf("%s_conf_sel.tex", file_stub),
              sprintf("%s: TPR and FDR of the confounding links selected at PIP $> 0.5$ by %s (small-world network with feedback loops and unmeasured confounding, $p = 10$; truth: the pairs whose $|\\sigma^*_{jh}|$ on the generating scale exceeds its mean over pairs).%s",
                      caption_what, tolower(factor_name), cap_tail(short_in(se[se$family == "confounding", ]))),
              sprintf("tab:%s_conf_sel", file_stub), sprintf("|c|r|%s", paste(rep("c|", 2 * length(cm)), collapse = "")), hdr, body, "\\scriptsize")
  }
}
selection_tables("noise", "Error variance", "$\\sigma^2$", function(l) as.numeric(sub("^errvar_", "", l)), "noise", "Error variance")
selection_tables("confounding", "Confounder loading scale", "Scale", function(l) as.numeric(sub("^conf_", "", l)), "confscale", "Strength of confounding")

## ---------------------------------------------------------------------------
## Main study: warnings of the package and the effect error of MR.RGM+ (mean and median), and the network AUC of
## MR.RGM on the replicates in which MR2G succeeded (the replicates over which MR2G is summarised)
## ---------------------------------------------------------------------------
if ("main" %in% names(exp_data)) local({
  rm_ <- res_all[res_all$experiment == "main", ]
  val <- function(sc, p_, n_, meth, fam, met) rm_$value[rm_$scenario == sc & rm_$p == p_ & rm_$n == n_ & rm_$method == meth &
                                                       rm_$family == fam & rm_$metric == met & !is.na(rm_$value)]
  pct <- function(x) if (!length(x)) "--" else sprintf("%.0f", 100 * mean(x))
  f3  <- function(x) if (!length(x)) "--" else sprintf("%.3f", x)
  mr2g_ok_diff <- function(sc, p_, n_) {
    a  <- rm_[rm_$scenario == sc & rm_$p == p_ & rm_$n == n_ & rm_$method == "MR.RGM" & rm_$family == "graph" & rm_$metric == "AUC" & !is.na(rm_$value), ]
    fr <- rm_$rep[rm_$scenario == sc & rm_$p == p_ & rm_$n == n_ & rm_$method == "MR2G" & rm_$family == "diagnostic" &
                  rm_$metric == "failed_fits" & !is.na(rm_$value) & rm_$value >= 1 & is.na(rm_$seconds)]
    if (!nrow(a) || !length(fr)) return("--")
    { d <- mean(a$value[!a$rep %in% fr]) - mean(a$value)          # math signs, not hyphens
      sprintf("$%s$%.3f (%d)", if (d < 0) "-" else "+", abs(d), length(unique(fr))) }
  }
  hdr <- paste(paste(c("", "", "", "\\multicolumn{2}{c|}{\\textbf{Any warning (\\%)}}", "\\multicolumn{3}{c|}{\\textbf{MR.RGM+ warnings (\\%)}}",
                       "\\multicolumn{2}{c|}{\\textbf{MR.RGM+ MAD}}", "\\textbf{MR.RGM AUC}"), collapse = " & "), "\\\\",
               paste(c("\\textbf{Case}", "\\textbf{$p$}", "\\textbf{$n$}", "\\textbf{MR.RGM}", "\\textbf{MR.RGM+}",
                       "\\textbf{Err.-var.}", "\\textbf{Spectral}", "\\textbf{Halves}", "\\textbf{Mean}", "\\textbf{Median}",
                       "\\textbf{diff.\\ (MR2G failed)}"), collapse = " & "))
  body <- c()
  for (sc in c("scalefree", "smallworld", "smallworld_pleio")) for (p_ in sort(unique(rm_$p))) {
    for (n_ in sort(unique(rm_$n))) {
      m <- val(sc, p_, n_, "MR.RGM+", "effect", "MeanAbsDev")
      body <- c(body, paste(c(if (n_ == min(rm_$n)) scen_title[[sc]] else "", if (n_ == min(rm_$n)) p_ else "", format(n_, big.mark = ","),
                              pct(val(sc, p_, n_, "MR.RGM", "diagnostic", "warn_any")), pct(val(sc, p_, n_, "MR.RGM+", "diagnostic", "warn_any")),
                              pct(val(sc, p_, n_, "MR.RGM+", "diagnostic", "warn_errvar")), pct(val(sc, p_, n_, "MR.RGM+", "diagnostic", "warn_spectral")),
                              pct(val(sc, p_, n_, "MR.RGM+", "diagnostic", "warn_halves")),
                              if (length(m)) f3(mean(m)) else "--", if (length(m)) f3(median(m)) else "--", mr2g_ok_diff(sc, p_, n_)),
                            collapse = " & "), "\\\\") }
    body <- c(body, "\\hline") }
  table_tex("main_diagnostics.tex",
            sprintf("Main study: percentage of the %d replicates in which \\texttt{RGM()} issued any warning (\\texttt{MR.RGM}, \\texttt{MR.RGM+}) and, for \\texttt{MR.RGM+}, the error-variance, spectral-radius and chain-halves warnings; mean and median over the replicates of the effect MAD of \\texttt{MR.RGM+}; difference between the mean network AUC of \\texttt{MR.RGM} over the replicates in which \\texttt{MR2G} succeeded and over all replicates (in parentheses: the number of replicates in which \\texttt{MR2G} stopped; --: none).", N_REP),
            "tab:main_diagnostics", "|l|r|r|c|c|c|c|c|c|c|c|", hdr, body, "\\scriptsize")
})

## ===========================================================================
## 3. Large n: table (n = 1e5, 1e6) and figure of AUC / MAD against n (main n plus large n, small-world p = 10)
## ===========================================================================
if ("largen" %in% names(exp_data)) {
  sl <- summ[summ$experiment == "largen", ]
  meths <- intersect(method_order, unique(sl$method[sl$family == "graph"]))
  ns <- sort(unique(sl$n))
  ## one row per (n, method): the network metrics and the effect errors; the confounding and instrument-selection
  ## AUCs of the MR.RGM family and the seconds per fit in a second table
  meths_e <- intersect(method_order, unique(sl$method[sl$family == "effect"]))
  rows <- expand.grid(method = union(meths, meths_e), n = ns, stringsAsFactors = FALSE)[, c("n", "method")]
  rows <- rows[order(match(rows$n, ns), match(rows$method, method_order)), ]; rows$group <- format(rows$n, big.mark = ",", scientific = FALSE)
  cols <- list("AUC" = c("graph", "AUC"), "TPR" = c("graph", "TPR"), "FDR" = c("graph", "FDR"), "MAD" = c("effect", "MeanAbsDev"),
               "MAD, true edges" = c("effect_edges", "MeanAbsDev"), "Mean $|$est.$|$, non-edges" = c("effect_nonedges", "MeanAbsDev"))
  long_table_tex("largen.tex", sprintf("Large samples ($n = 10^5$ and $10^6$; small-world network with feedback loops and unmeasured confounding, $p = 10$): causal network recovery (AUC; TPR and FDR of the selected graph) and causal effect estimation (MAD over all off-diagonal entries, on the true edges, and the mean absolute estimate on the non-edges; All-zero: the estimate 0).%s", cap_tail(short_in(sl[sl$family %in% c("graph", "effect", "effect_edges", "effect_nonedges"), ]))),
                 "tab:largen", rows, "$n$", cols, function(i, cl) cell(sl, "largen", "smallworld", 10, rows$n[i], rows$method[i], cl[1], cl[2]), "\\scriptsize")
  sec <- aggregate(seconds ~ n + method, data = sl[sl$family == "graph" & !is.na(sl$seconds), ], FUN = mean)
  secs <- function(n_, m) { v <- sec$seconds[sec$n == n_ & sec$method == m]; if (!length(v)) "--" else if (v < 0.5) "$<1$" else format(round(v), big.mark = ",", scientific = FALSE) }
  body <- c()
  for (n_ in ns) { for (k in seq_along(meths)) body <- c(body, paste(c(if (k == 1) format(n_, big.mark = ",", scientific = FALSE) else "", method_tex[[meths[k]]],
      cell(sl, "largen", "smallworld", 10, n_, meths[k], "confounding", "AUC"), cell(sl, "largen", "smallworld", 10, n_, meths[k], "iv_selection", "AUC"), secs(n_, meths[k])), collapse = " & "), "\\\\")
    body <- c(body, "\\hline") }
  table_tex("largen_seconds.tex", "Large samples (small-world network, $p = 10$): confounding structure recovery (AUC; MR.RGM and MR.RGM+), instrument--trait selection (AUC; MR.RGM+) and mean wall-clock seconds per fit (one core per replicate; MR.RGM family: 50,000 iterations). Entries are mean (SD) over 100 replicates.", "tab:largen_seconds",
            "|r|l|c|c|r|", paste(c("\\textbf{$n$}", "\\textbf{Method}", "\\textbf{Conf.\\ AUC}", "\\textbf{IV AUC}", "\\textbf{Seconds}"), collapse = " & "), body)
  ## figure over all n (main + largen), small-world p = 10
  sa <- summ[((summ$experiment == "main" & summ$scenario == "smallworld" & summ$p == 10 & summ$setting_label == "base") | summ$experiment == "largen") &
             summ$method %in% meths, ]
  dd <- do.call(rbind, lapply(METRICS, function(mm) { x <- sa[sa$family == mm$family & sa$metric == mm$metric & !is.na(sa$mean), ]
    if (!nrow(x)) return(NULL); x$Metric <- mm$label; x }))
  if (!is.null(dd) && nrow(dd)) {
    dd$Metric <- droplevels(factor(dd$Metric, levels = vapply(METRICS, `[[`, "", "label")))
    fig_note("largen.pdf", dd, from_summary = TRUE)
    dd$Method <- meth_f(dd$method)
    g <- ggplot(dd, aes(n, mean, colour = Method, shape = Method, linetype = Method, group = Method)) + geom_line(linewidth = if (THESIS) 0.45 else 0.8) +
      geom_pointrange(aes(ymin = mean - sd, ymax = mean + sd), size = if (THESIS) 0.25 else 0.4, linewidth = if (THESIS) 0.3 else 0.5) +
      scale_methods() + scale_x_log10(breaks = 10^(3:6), labels = c(expression(10^3), expression(10^4), expression(10^5), expression(10^6))) +
      facet_wrap(~ Metric, scales = "free_y", nrow = 2) +
      labs(x = "Sample size n (log scale)", y = NULL, title = "Performance against sample size (small-world, p = 10)") + theme_line +
      theme(legend.key.width = unit(1.6, "lines"))
    save_fig(g, "largen.pdf", 14, 10)
  }
}

## ===========================================================================
## 4. Sensitivity to the hyperparameters of MR.RGM and MR.RGM+
## ===========================================================================
if ("sensitivity" %in% names(exp_data)) {
  ss <- summ[summ$experiment == "sensitivity", ]
  labs_ <- c("base", setdiff(unique(ss$setting_label), "base")); ns <- sort(unique(ss$n))
  labs_ <- labs_[order(match(sub("_.*", "", labs_), c("base", "Beta", "nu1", "nu2")), -as.numeric(ifelse(grepl("^nu", labs_), sub("^nu[12]_", "", labs_), 0)))]   # base, Beta, nu1 (large to small), nu2
  sci <- function(x) { x <- as.numeric(x); e <- round(log10(x)); m <- x / 10^e
    if (abs(m - 1) < 1e-9) sprintf("10^{%d}", e) else sprintf("%g \\times 10^{%d}", m, e) }
  pretty_lab <- function(l) { if (l == "base") return("Default: Beta(3,1), $\\nu_1 = 10^{-3}$, $\\nu_2 = 10^{-4}$")
    if (grepl("^Beta_", l)) { ab <- strsplit(sub("^Beta_", "", l), "_")[[1]]; return(sprintf("Beta(%s,%s)", ab[1], ab[2])) }
    if (grepl("^nu1_", l)) return(sprintf("$\\nu_1 = %s$", sci(sub("^nu1_", "", l))))
    if (grepl("^nu2_", l)) return(sprintf("$\\nu_2 = %s$", sci(sub("^nu2_", "", l))))
    l }
  ## nu2 is the spike variance of the instrument effects, which only MR.RGM+ selects: the nu2 rows of MR.RGM and
  ## MR.RGM_NoConf would repeat their base rows exactly (same data and seeds), so they are shown as "--"
  rgm <- function(l, n_, fam, met, m = "MR.RGM") if (grepl("^nu2_", l)) "--" else cell(ss, "sensitivity", "smallworld", 10, n_, m, fam, met, l)
  cap0 <- "Sensitivity of %s to the prior inclusion probability (Beta($a_\\rho, b_\\rho$)) and the spike variances $\\nu_1$ (causal effects) and $\\nu_2$ (instrument effects), one changed at a time (small-world network with feedback loops and unmeasured confounding, $p = 10$); FDR of the edges selected at PIP $> 0.5$. %sThe competing methods do not use these hyperparameters; the data are those of the main study, so their results are those of the main study for the same case and $n$.%s"
  ## table 1: MR.RGM and MR.RGM_NoConf (known instrument map); table 2: MR.RGM+
  hdr <- paste(c("\\textbf{Hyperparameters}", "\\textbf{$n$}", "\\textbf{AUC}", "\\textbf{FDR}", "\\textbf{MAD}", "\\textbf{Conf.\\ AUC}", "\\textbf{AUC, NoConf}", "\\textbf{MAD, NoConf}"), collapse = " & ")
  body <- c()
  for (l in labs_) { for (n_ in ns) body <- c(body, paste(c(if (n_ == ns[1]) pretty_lab(l) else "", format(n_, big.mark = ","),
      rgm(l, n_, "graph", "AUC"), rgm(l, n_, "graph", "FDR"), rgm(l, n_, "effect", "MeanAbsDev"), rgm(l, n_, "confounding", "AUC"),
      rgm(l, n_, "graph", "AUC", "MR.RGM_NoConf"), rgm(l, n_, "effect", "MeanAbsDev", "MR.RGM_NoConf")), collapse = " & "), "\\\\")
    body <- c(body, "\\hline") }
  table_tex("sensitivity.tex", sprintf(cap0, "MR.RGM and MR.RGM\\_NoConf", "$\\nu_2$ applies only to MR.RGM+ (Table~\\ref{tab:sensitivity_plus}); MR.RGM and MR.RGM\\_NoConf (known instrument map, no selection) do not use it, so their entries in the $\\nu_2$ rows are marked ``--''. ", cap_tail(short_in(ss[ss$method != "MR.RGM+", ]))),
            "tab:sensitivity", "|l|r|c|c|c|c|c|c|", hdr, body, "\\scriptsize")
  hdr <- paste(c("\\textbf{Hyperparameters}", "\\textbf{$n$}", "\\textbf{AUC}", "\\textbf{FDR}", "\\textbf{MAD}", "\\textbf{Conf.\\ AUC}", "\\textbf{IV AUC}"), collapse = " & ")
  body <- c()
  for (l in labs_) { for (n_ in ns) body <- c(body, paste(c(if (n_ == ns[1]) pretty_lab(l) else "", format(n_, big.mark = ","),
      cell(ss, "sensitivity", "smallworld", 10, n_, "MR.RGM+", "graph", "AUC", l), cell(ss, "sensitivity", "smallworld", 10, n_, "MR.RGM+", "graph", "FDR", l),
      cell(ss, "sensitivity", "smallworld", 10, n_, "MR.RGM+", "effect", "MeanAbsDev", l), cell(ss, "sensitivity", "smallworld", 10, n_, "MR.RGM+", "confounding", "AUC", l),
      cell(ss, "sensitivity", "smallworld", 10, n_, "MR.RGM+", "iv_selection", "AUC", l)), collapse = " & "), "\\\\")
    body <- c(body, "\\hline") }
  table_tex("sensitivity_plus.tex", sprintf(cap0, "MR.RGM+", "", cap_tail(short_in(ss[ss$method == "MR.RGM+", ]))),
            "tab:sensitivity_plus", "|l|r|c|c|c|c|c|", hdr, body, "\\scriptsize")
}

## ===========================================================================
## 5. Exposure -> outcome design (roles): exposure->outcome block, outcome->outcome block, confounding among outcomes
## ===========================================================================
if ("roles" %in% names(exp_data)) {
  sr <- summ[summ$experiment == "roles", ]; rr <- res[res$experiment == "roles", ]
  ns <- sort(unique(sr$n))
  meths <- intersect(method_order, unique(sr$method[sr$family == "graph_xy"]))
  rows <- expand.grid(method = meths, n = ns, stringsAsFactors = FALSE)[, c("n", "method")]; rows$group <- format(rows$n, big.mark = ",")
  cap_common <- "MR2 and MrDAG are given the true split into exposures and outcomes; the other methods are not. The pairwise methods (MendelianRandomization, OneSampleMR, mrbayes) estimate the total effect of each trait on each other trait."
  cols <- list("AUC" = c("graph", "AUC"), "FDR" = c("graph", "FDR"), "MAD" = c("effect", "MeanAbsDev"),
               "AUC, out.\\ $\\to$ out." = c("graph_yy", "AUC"), "Conf.\\ AUC, outcomes" = c("confounding", "AUC"))
  long_table_tex("roles.tex", sprintf("Exposure--outcome design (4 exposures, 4 outcomes with a causal DAG among the outcomes, unmeasured confounding): recovery of the whole graph (AUC; FDR at the decision rule PIP $>0.5$ for the Bayesian methods and $p<0.05$ for the others; MAD over all off-diagonal entries), of the outcome $\\to$ outcome effects (AUC; the network methods, MVMR and the pairwise methods estimate them, MR2 does not, MrDAG through its DAG among the outcomes), and of the confounding structure among the outcomes (AUC; MR2 through its residual correlation). %s%s", cap_common, cap_tail(short_in(sr[sr$family %in% c("graph", "effect", "graph_yy", "confounding"), ]))),
                 "tab:roles", rows, "$n$", cols, function(i, cl) cell(sr, "roles", "roles", unique(sr$p)[1], rows$n[i], rows$method[i], cl[1], cl[2]), "\\scriptsize")
  cols <- list("AUC" = c("graph_xy", "AUC"), "TPR" = c("graph_xy", "TPR"), "FDR" = c("graph_xy", "FDR"), "MAD" = c("effect_xy", "MeanAbsDev"))
  long_table_tex("roles_xy.tex", sprintf("Exposure--outcome design: recovery of the exposure $\\to$ outcome effects (AUC; TPR and FDR at the decision rule PIP $>0.5$ for the Bayesian methods and $p<0.05$ for the others; MAD). %s%s", cap_common, cap_tail(short_in(sr[sr$family %in% c("graph_xy", "effect_xy"), ]))),
                 "tab:roles_xy", rows, "$n$", cols, function(i, cl) cell(sr, "roles", "roles", unique(sr$p)[1], rows$n[i], rows$method[i], cl[1], cl[2]))
  ## figure: boxplots of the three AUCs per method
  d <- rr[rr$family %in% c("graph_xy", "graph_yy", "confounding") & rr$metric == "AUC" & !is.na(rr$value) & rr$method %in% meths, ]
  if (nrow(d)) {
    for (fm in unique(d$family)) fig_note(sprintf("roles_auc.pdf [%s]", fm), d[d$family == fm, ])
    d$Block <- factor(c(graph_xy = "Exposure to outcome", graph_yy = "Outcome to outcome", confounding = "Confounding (outcomes)")[d$family],
                      levels = c("Exposure to outcome", "Outcome to outcome", "Confounding (outcomes)"))
    d$Method <- meth_f(d$method); d$N <- factor(d$n, levels = ns, labels = n_lab(ns))
    g <- ggplot(d, aes(Method, value, fill = Method)) + geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey60") +
      geom_boxplot(colour = "grey15", outlier.size = 1.5) + scale_fill_manual(values = method_col_d) + guides(fill = "none") +
      facet_grid(N ~ Block, scales = "free_x", space = "free_x", labeller = labeller(Block = label_wrap_gen(12))) + labs(title = "Exposures and outcomes with known roles", x = "Method", y = "Area Under Curve (AUC)") + theme_box
    save_fig(g, "roles_auc.pdf", 15, 10)
  }
}

## ===========================================================================
## 5b. Failed fits: number of failed fits per method and setting (every table and figure above summarises the
## successful fits only; their number is given in the captions). Failed fits are flagged in
## tables/simulation_all_replicates.csv (column failed_fit) and their error messages are in results_<experiment>/errors.log.
## ===========================================================================
{
  fd <- res_all[res_all$family == "diagnostic" & res_all$metric == "failed_fits" & res_all$method != "data", ]
  if (nrow(fd)) {
    fd$failed <- fd$failed_fit
    fd$cause  <- ifelse(fd$failed, fit_cause(fd), "")
    cause_lab <- c(significance = "no significant instrument", peeling = "peeling", numerical = "numerical", time = "time", error = "error")
    cause_text <- function(r) {                                   # e.g. "peeling 3, numerical 2"
      x <- fd$cause[fd$failed & fd$experiment == r$experiment & fd$scenario == r$scenario & fd$p == r$p & fd$n == r$n &
                    fd$setting_label == r$setting_label & fd$method == r$method]
      tb <- table(factor(x, levels = names(cause_lab))); tb <- tb[tb > 0]
      if (length(tb) == 1 && names(tb) == "significance") return(cause_lab[["significance"]])
      paste(sprintf("%s %d", cause_lab[names(tb)], as.integer(tb)), collapse = ", ")
    }
    fr <- aggregate(failed ~ experiment + scenario + p + n + setting_label + method, data = fd, FUN = function(x) c(k = sum(x), m = length(x)))
    fr <- do.call(data.frame, fr)
    fr <- fr[fr$failed.k > 0, ]
    if (nrow(fr)) {
      exp_lab <- c(main = "Main study", noise = "Error variance", confounding = "Confounding strength", largen = "Large samples")
      set_lab <- function(l) if (l == "base") "default" else if (grepl("^errvar_", l)) sprintf("$\\sigma^2 = %s$", sub("^errvar_", "", l)) else
        if (grepl("^conf_", l)) sprintf("scale %s", sub("^conf_", "", l)) else gsub("_", "\\\\_", l)
      fr <- fr[order(match(fr$experiment, names(exp_lab)), match(fr$method, method_order), fr$scenario, fr$p, fr$n, fr$setting_label), ]
      body <- unlist(lapply(seq_len(nrow(fr)), function(i) { r <- fr[i, ]
        first <- i == 1 || fr$experiment[i] != fr$experiment[i - 1]
        c(if (first && i > 1) "\\hline",
          paste(c(if (first) (if (r$experiment %in% names(exp_lab)) exp_lab[[r$experiment]] else r$experiment) else "",
                  if (r$scenario %in% names(scen_title)) scen_title[[r$scenario]] else r$scenario, r$p,
                  format(r$n, big.mark = ",", scientific = FALSE), set_lab(r$setting_label),
                  if (r$method %in% names(method_tex)) method_tex[[r$method]] else r$method,
                  sprintf("%d / %d", r$failed.k, r$failed.m), sprintf("%d", r$failed.m - r$failed.k), cause_text(r)), collapse = " & "), "\\\\") }))
      body <- c(body, "\\hline")
      table_tex("failures.tex",
                paste0("Fits without a result, by experiment and setting: number of such fits out of the replicates, number of successful fits over which the method is summarized in every table and figure, and cause (MR2G: ", fail_reason[["MR2G"]],
                       "; GrIVET: peeling, its peeling step could not order every trait, which a directed cycle can prevent; numerical, it stopped with a numerical error; time, it did not finish within 48 hours)."),
                "tab:failures", "|l|l|r|r|l|l|c|c|l|",
                paste(c("\\textbf{Experiment}", "\\textbf{Case}", "\\textbf{$p$}", "\\textbf{$n$}", "\\textbf{Setting}", "\\textbf{Method}",
                        "\\textbf{Failed}", "\\textbf{Successful}", "\\textbf{Cause}"), collapse = " & "),
                body)
    }
  }
}
## Caption notes for the figures (to be written under the figures in the paper)
writeLines(display_text(if (length(FIG_NOTES)) FIG_NOTES else "Every method has all replicates in every figure."), file.path(TAB_DIR, "figure_notes.txt"))

## ===========================================================================
## 6. Timing: mean seconds per fit by method, experiment main (p = 5 and 10, each n)
## ===========================================================================
if ("main" %in% names(exp_data)) {
  st <- summ[summ$experiment == "main" & summ$family == "graph" & summ$metric == "AUC" & !is.na(summ$seconds), ]
  if (nrow(st)) {
    meths <- intersect(method_order, unique(st$method)); ns <- sort(unique(st$n))
    body <- c()
    for (p_ in sort(unique(st$p))) { for (n_ in ns) body <- c(body, paste(c(if (n_ == ns[1]) p_ else "", format(n_, big.mark = ","),
        vapply(meths, function(m) { v <- st$seconds[st$p == p_ & st$n == n_ & st$method == m & st$scenario == "smallworld"]; if (!length(v)) "--" else if (v < 0.05) "$<0.1$" else format(round(v, 1), nsmall = 1, big.mark = ",", scientific = FALSE) }, character(1))), collapse = " & "), "\\\\")
      body <- c(body, "\\hline") }
    table_tex("timing.tex", "Mean wall-clock seconds per fit (small-world case, one core per replicate; MCMC methods: 50,000 iterations; MendelianRandomization, OneSampleMR and mrbayes: all ordered pairs; MrDAG: one fit per outcome). The three MendelianRandomization estimators are computed in one call and share its time.",
              "tab:timing", sprintf("|c|r|%s", paste(rep("c|", length(meths)), collapse = "")),
              paste(c("\\textbf{$p$}", "\\textbf{$n$}", sprintf("\\textbf{%s}", method_tex[meths])), collapse = " & "), body, "\\scriptsize")
  }
}
message("Wrote ", length(list.files(FIG_DIR, pattern = "\\.pdf$")), " figures and ", length(list.files(TAB_DIR, pattern = "\\.tex$")), " tables.")
TRUE
}, error = function(e) {
  message("make_figures_tables.R stopped with an error: ", conditionMessage(e))
  stop(e)
})
