###############################################################################
## RunTimeAnalysis.R
##
## Runtime benchmark (Main text, "Scalability Test"): scale-free network,
## n = 30,000, p in {2, 5, 10, 20}; every method run 20 times; median
## wall-clock time reported.
##
## Uses the SAME data generator and the SAME method settings as the simulation
## study (sim_functions.R), so the timings correspond exactly to what was run:
##   MR.RGM family: 50,000 iterations, 10,000 burn-in, thinning 10
##   mrbayes:       2 chains, 10,000 burn-in, 50,000 iterations per pair
##   MendelianRandomization: 10,000 bootstrap draws (median-based SEs)
##
## Usage (from the simulation folder):
##   setwd(".../reproducibility/simulation"); source("Runtime Analysis/RunTimeAnalysis.R")
## Outputs: Runtime Analysis/runtime_vs_p_n30000.csv and figures/runtime_vs_p_n30000.pdf
###############################################################################

source("sim_functions.R")
suppressPackageStartupMessages({ library(microbenchmark); library(ggplot2) })

if (!exists("RUNTIME_TIMES")) RUNTIME_TIMES <- 20
set    <- SIM_SETTINGS
n      <- 30000
p_vals <- c(2, 5, 10, 20)
dir.create("figures", showWarnings = FALSE)

rows <- list()
for (p in p_vals) {
  set.seed(truth_seed(1));                       truth <- gen_truth("scalefree", p, set)
  set.seed(as.integer(9e6 + p));                 dat   <- gen_data(truth, n, set)
  X <- dat$X; Y <- dat$Y; D <- truth$D; Dfull <- matrix(1, p, ncol(X))
  message(sprintf("[%s] p = %d: %d runs per method", format(Sys.time(), "%H:%M:%S"), p, RUNTIME_TIMES))
  mb <- microbenchmark(
    `MR.RGM`                 = fit_rgm(X, Y, D,     "SSSL",     set),
    `MR.RGM+`                = fit_rgm(X, Y, Dfull, "SSSL",     set),
    `MR.RGM_NoConf`          = fit_rgm(X, Y, D,     "diagonal", set),
    `MendelianRandomization` = fit_mr_package(X, Y, D, set),
    `OneSampleMR`            = fit_onesamplemr(X, Y, D),
    `mrbayes`                = fit_mrbayes(X, Y, D, seed = 1L, set = set),
    times = RUNTIME_TIMES)
  med <- tapply(mb$time, mb$expr, median) / 1e9                     # seconds
  rows[[length(rows) + 1]] <- data.frame(Method = names(med), p = p, median_sec = as.numeric(med))
}
rt <- do.call(rbind, rows)
write.csv(rt, file.path("Runtime Analysis", sprintf("runtime_vs_p_n%d.csv", n)), row.names = FALSE)
print(reshape(rt, idvar = "Method", timevar = "p", direction = "wide"))

## Plot (same style as the manuscript figure)
method_order <- c("MR.RGM", "MR.RGM+", "MR.RGM_NoConf", "MendelianRandomization", "OneSampleMR", "mrbayes")
rt$Method <- factor(rt$Method, levels = method_order)
end_pts   <- do.call(rbind, lapply(split(rt, rt$Method), function(d) d[which.max(d$p), ]))
g <- ggplot(rt, aes(p, median_sec, colour = Method)) +
  geom_line(aes(group = Method), linewidth = 1.2) + geom_point(size = 2.6, stroke = 0.5) +
  geom_text(data = end_pts, aes(label = Method), hjust = -0.05, vjust = 0.5, size = 5.3,
            fontface = "bold", show.legend = FALSE) +
  scale_colour_manual(values = c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#000000", "#CC79A7"), guide = "none") +
  scale_x_continuous(breaks = p_vals, expand = expansion(mult = c(0.02, 0.22))) +
  scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  labs(title = "Median runtime versus number of responses (p)",
       subtitle = sprintf("Fixed sample size n = %s · medians across %d runs", format(n, big.mark = ","), RUNTIME_TIMES),
       x = "Number of Responses (p)", y = "Median Runtime (seconds)") +
  theme_classic(base_size = 12) +
  theme(plot.title.position = "plot", plot.title = element_text(face = "bold", size = 24),
        plot.subtitle = element_text(margin = margin(b = 6), size = 18),
        axis.title = element_text(face = "bold", size = 24), axis.text = element_text(size = 24),
        axis.ticks.length = unit(3, "pt"), panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.4),
        panel.grid.minor = element_blank(), plot.margin = margin(8, 42, 8, 8)) +
  coord_cartesian(clip = "off")
ggsave(file.path("figures", sprintf("runtime_vs_p_n%d.pdf", n)), g, width = 14, height = 8, units = "in")
