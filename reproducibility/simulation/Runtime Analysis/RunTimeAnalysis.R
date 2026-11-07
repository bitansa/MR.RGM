###############################################################################
## RunTimeAnalysis.R
##
## Runtime benchmark ("Scalability"): scale-free network, n = 30,000, p in {2, 5, 10, 20}; every method run
## RUNTIME_TIMES (100) times on the same data; median wall-clock time reported.
##
## Uses the same data generator, standardization and method settings as the simulation study (sim_functions.R),
## so the timings correspond exactly to what was run:
##   MR.RGM family: 50,000 iterations, 10,000 burn-in, thinning 10 (sufficient statistics: cost does not grow with n)
##   MR2G:          package defaults (100 perturbations, ScreenAug)
##   MVMR:          MendelianRandomization::mr_mvivw, one outcome at a time
##   mrbayes:       3 chains (package default), 10,000 burn-in + 40,000 iterations per chain, per ordered pair
##   MendelianRandomization: 10,000 bootstrap draws (median-based SEs) per ordered pair
##   OneSampleMR:   TSRI per ordered pair
##   MrDAG:         50,000 iterations, 10,000 burn-in, thinning 10, one outcome at a time with the other traits as
##                  exposures (by_outcome() in sim_functions.R, as in the simulation designs without roles)
##   GrIVET:        data-driven pipeline with the tuning grid of the simulation study, on the unstandardized data
##   MR2 needs an exposure/outcome split with at least two responses and is not part of this benchmark (its seconds
##   per fit in the roles design are stored in the results of that experiment)
## Every other method of the simulation study is timed. Methods whose package is not installed are skipped, and a method that
## cannot be fitted to a data set is reported with median_sec = NA.
##
## Usage (from the simulation folder):
##   source("Runtime Analysis/RunTimeAnalysis.R")
## Outputs: Runtime Analysis/runtime_vs_p_n30000.csv (medians), runtime_vs_p_n30000_all_runs.csv (every run),
##          figures/runtime_vs_p_n30000.pdf (+ png). Redraw only: RUNTIME_PLOT_ONLY <- TRUE before sourcing.
###############################################################################

source("sim_functions.R")
suppressPackageStartupMessages({ library(microbenchmark); library(ggplot2) })

if (!exists("RUNTIME_TIMES")) RUNTIME_TIMES <- 100     # runs per method and p (same number as the simulation replicates)
set    <- SIM_SETTINGS
n      <- 30000
p_vals <- c(2, 5, 10, 20)
FIG_DIR <- "figures"; stub <- ""
## RUNTIME_SMOKE <- TRUE: 1 run per method, p in {2, 5}, short chains and few perturbations, separate output files
if (exists("RUNTIME_SMOKE") && isTRUE(RUNTIME_SMOKE)) {
  RUNTIME_TIMES <- 1; p_vals <- c(2, 5); FIG_DIR <- "figures_smoke"; stub <- "_smoke"
  set$mcmc <- list(nIter = 2000, nBurnin = 500, Thin = 2); set$mr_boot <- 200
  set$mrbayes <- list(n.chains = 3, n.burn = 200, n.iter = 800)
  set$mr2g <- modifyList(set$mr2g, list(num_pert = 5, num_pert_M = 5))
  set$mr2  <- modifyList(set$mr2,  list(niter = 600, burnin = 200, thin = 2))
  set$mrdag <- modifyList(set$mrdag, list(niter = 600, burnin = 200, thin = 2))
}
dir.create(file.path(FIG_DIR, "png"), recursive = TRUE, showWarnings = FALSE)

csv <- file.path("Runtime Analysis", sprintf("runtime_vs_p_n%d%s.csv", n, stub))
## RUNTIME_PLOT_ONLY <- TRUE: redraw the figure from the saved CSV without rerunning the benchmark
if (exists("RUNTIME_PLOT_ONLY") && isTRUE(RUNTIME_PLOT_ONLY)) {
  if (!file.exists(csv)) stop("No saved timings: ", csv, ". Run the benchmark first.")
  rt <- read.csv(csv)
  ## A method whose fit stopped with an error in some blocks of runs of the cluster benchmark has median_sec = NA in
  ## the summary, and the runs of those blocks are missing from the per-run file; its median is then taken over the
  ## runs that completed, and the counts are written to runtime_notes.txt (for the caption)
  ar_f <- sub("\\.csv$", "_all_runs.csv", csv); RT_NOTES <- character(0)
  if (file.exists(ar_f)) {
    ar <- read.csv(ar_f)
    for (i in which(is.na(rt$median_sec))) {
      s <- ar$sec[ar$Method == rt$Method[i] & ar$p == rt$p[i] & !is.na(ar$sec)]
      if (length(s)) {
        rt$median_sec[i] <- median(s)
        RT_NOTES <- c(RT_NOTES, sprintf("%s, p = %d: median over the %d of %d runs that completed (the other runs were not timed because a fit stopped with an error).",
                                        rt$Method[i], rt$p[i], length(s), RUNTIME_TIMES))
      }
    }
  }
  writeLines(display_text(if (length(RT_NOTES)) RT_NOTES else "Every method has all runs at every p."), file.path("Runtime Analysis", "runtime_notes.txt"))
} else {
## Restartable: the timings of every finished p are saved in runtime_partial*.rds and skipped on a restart
partial <- file.path("Runtime Analysis", sprintf("runtime_partial%s.rds", stub))
rows <- list(); raw <- list()
if (file.exists(partial)) {
  pr <- readRDS(partial)
  if (identical(pr$sampler, SAMPLER_TAG)) {
    rows <- pr$rows; raw <- pr$raw
    message("Restart: timings for p = ", paste(vapply(rows, function(r) r$p[1], numeric(1)), collapse = ", "), " read from ", partial)
  } else message("Ignoring ", partial, ": no current sampler tag; the benchmark starts from the beginning.")
}
for (p in setdiff(p_vals, vapply(rows, function(r) r$p[1], numeric(1)))) {
  rd <- runtime_data(p, n, set); X <- rd$X; Y <- rd$Y; D <- rd$D; Dfull <- rd$Dfull; ms <- rd$ms   # sim_functions.R
  X_raw <- rd$X_raw; Y_raw <- rd$Y_raw; sY <- rd$sY                                                  # GrIVET
  message(sprintf("[%s] p = %d: %d runs per method", format(Sys.time(), "%H:%M:%S"), p, RUNTIME_TIMES))
  exprs <- RUNTIME_EXPRS                                              # the call of every method (sim_functions.R)
  for (m in intersect(names(exprs), names(METHOD_PACKAGE))) if (!has_pkg(METHOD_PACKAGE[[m]])) exprs[[m]] <- NULL
  ## A competitor that fails on this data set is not timed (otherwise the time to its error would be recorded):
  ## one trial fit each, before the benchmark; a failing method is reported with median_sec = NA
  failed_m <- character(0)
  for (m in intersect(names(exprs), c("MVMR", "MR2G", "MendelianRandomization", "OneSampleMR", "mrbayes", "MR2", "MrDAG", "GrIVET")))
    if (is.null(eval(exprs[[m]]))) failed_m <- c(failed_m, m)
  if (length(failed_m)) { message("   not timed (fit failed): ", paste(failed_m, collapse = ", ")); exprs[failed_m] <- NULL }
  mb <- microbenchmark(list = exprs, times = RUNTIME_TIMES)
  med <- tapply(mb$time, mb$expr, median) / 1e9                     # seconds
  rows[[length(rows) + 1]] <- rbind(data.frame(Method = names(med), p = p, median_sec = as.numeric(med)),
                                    if (length(failed_m)) data.frame(Method = failed_m, p = p, median_sec = NA_real_))
  raw[[length(raw) + 1]]   <- data.frame(Method = as.character(mb$expr), p = p, sec = mb$time / 1e9)
  saveRDS(list(rows = rows, raw = raw, sampler = SAMPLER_TAG,
               MR.RGM_build = unlist(utils::packageDescription("MR.RGM")[c("Version", "Packaged", "RemoteSha", "Built")])), partial)                    # progress: a restart skips this p
}
rt <- do.call(rbind, rows)
rt <- rt[order(rt$p), ]
attr(rt, "package_versions") <- package_versions()
write.csv(rt, csv, row.names = FALSE)
write.csv(do.call(rbind, raw), sub("\\.csv$", "_all_runs.csv", csv), row.names = FALSE)   # every timed run
writeLines(c(paste("Run on", format(Sys.time())), paste(names(attr(rt, "package_versions")), attr(rt, "package_versions"), sep = " ")),
           file.path("Runtime Analysis", sprintf("runtime_versions%s.txt", stub)))
}
print(reshape(rt, idvar = "Method", timevar = "p", direction = "wide"))

## Plot (same style as the manuscript figure)
method_order <- c("MR.RGM", "MR.RGM+", "MR.RGM_NoConf", "MR2G", "MVMR", "MR2", "MrDAG", "GrIVET", "MendelianRandomization", "OneSampleMR", "mrbayes")
rt_lv     <- intersect(method_order, unique(rt$Method))
rt$Method <- factor(rt$Method, levels = rt_lv, labels = display_text(rt_lv))   # labels: displayed method names (sim_functions.R)
## same colours as make_figures_tables.R, keyed by the displayed method names
rt_col <- c("MR.RGM" = "#0072B2", "MR.RGM+" = "#E69F00", "MR.RGM_NoConf" = "#56B4E9",
            "MR2G" = "#D55E00", "MVMR" = "#009E73", "MR2" = "#CC79A7", "MrDAG" = "#F0E442", "GrIVET" = "#882255",
            "MendelianRandomization" = "#999999", "OneSampleMR" = "#000000", "mrbayes" = "#8C510A")
rt_col <- setNames(rt_col, display_text(names(rt_col)))
rt_plot   <- rt[!is.na(rt$median_sec), ]                                   # methods not timed at some p (NA) are left out there
end_pts   <- do.call(rbind, lapply(split(rt_plot, rt_plot$Method), function(d) d[which.max(d$p), ]))
g <- ggplot(rt_plot, aes(p, median_sec, colour = Method)) +
  geom_line(aes(group = Method), linewidth = 1.2) + geom_point(size = 2.6, stroke = 0.5) +
  geom_text(data = end_pts, aes(label = Method), hjust = -0.05, vjust = 0.5, size = if (exists("THESIS_STYLE")) 2.5 else 5.3,
            fontface = "bold", show.legend = FALSE) +
  scale_colour_manual(values = rt_col, guide = "none") +                         # same colours as make_figures_tables.R
  scale_x_continuous(breaks = p_vals, expand = expansion(mult = c(0.02, 0.22))) +
  scale_y_log10(breaks = 10^(-2:4), labels = c("0.01", "0.1", "1", "10", "100", "1,000", "10,000")) +
  labs(title = display_text("Median runtime versus number of traits (p)"),
       subtitle = sprintf("Fixed sample size n = %s · medians across %d runs · log scale", format(n, big.mark = ","), RUNTIME_TIMES),
       x = display_text("Number of Traits (p)"), y = "Median time per fit (seconds, log scale)") +
  theme_classic(base_size = 12) +
  theme(plot.title.position = "plot", plot.title = element_text(face = "bold", size = 24),
        plot.subtitle = element_text(margin = margin(b = 6), size = 18),
        axis.title = element_text(face = "bold", size = 24), axis.text = element_text(size = 24),
        axis.ticks.length = unit(3, "pt"), panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.4),
        panel.grid.minor = element_blank(), plot.margin = margin(8, 42, 8, 8)) +
  coord_cartesian(clip = "off")
## Print-size style: off unless a figure driver has defined
## THESIS_STYLE (then: theme_bw(base_size = 9), sizes from the driver's CSV, cairo_pdf, no png). The default figures are unaffected.
STYLE <- if (exists("THESIS_STYLE")) THESIS_STYLE else NULL
if (!is.null(STYLE)) {
  ## small figure: a legend instead of the labels at the line ends, which overlap where the curves are close
  g$layers <- Filter(function(l) !inherits(l$geom, "GeomText"), g$layers)
  g <- g + scale_colour_manual(values = rt_col) +
    scale_x_continuous(breaks = p_vals, expand = expansion(mult = c(0.03, 0.03))) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    theme_bw(base_size = 9) +
    theme(plot.title.position = "plot", plot.title = element_text(face = "bold"), axis.title = element_text(face = "bold"),
          panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), plot.margin = margin(4, 6, 4, 4),
          legend.position = "top", legend.title = element_blank(), legend.margin = margin(0, 0, 0, 0))
  if (exists("DROP_TITLES") && isTRUE(DROP_TITLES)) g <- g + labs(subtitle = NULL)   # the LaTeX caption gives n and the number of runs
  STYLE$save(g, file.path(FIG_DIR, sprintf("runtime_vs_p_n%d.pdf", n)), 14, 8)
} else {
ggsave(file.path(FIG_DIR, sprintf("runtime_vs_p_n%d.pdf", n)), g, width = 14, height = 8, units = "in")
ggsave(file.path(FIG_DIR, "png", sprintf("runtime_vs_p_n%d.png", n)), g, width = 14, height = 8, units = "in", dpi = 300)
}
