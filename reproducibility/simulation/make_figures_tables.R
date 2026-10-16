###############################################################################
## make_figures_tables.R
##
## Reads results/*.rds (from run_simulations.R) and writes
##   figures/  every simulation figure of the paper (same file names as in the
##             manuscript's \includegraphics)
##   tables/   every simulation table of the Supplement as a .tex file (same
##             labels and captions), to be \input{} into Supplementary.tex
##   tables/simulation_summary.csv   all means, SDs and replicate counts
##
## No number is typed by hand: everything is computed from the result files.
##
## Usage (from this folder):  source("make_figures_tables.R")
###############################################################################

suppressPackageStartupMessages({ library(ggplot2); library(scales) })

if (!exists("RESULTS_DIR")) RESULTS_DIR <- "results"
dir.create("figures", showWarnings = FALSE); dir.create("tables", showWarnings = FALSE)

files <- list.files(RESULTS_DIR, pattern = "^sim_.*\\.rds$", full.names = TRUE)
if (!length(files)) stop("No result files found in ", RESULTS_DIR)
res <- do.call(rbind, lapply(files, readRDS))
N_REP <- 20

## ---------------------------------------------------------------------------
## Labels and order (identical in figures, tables and text)
## ---------------------------------------------------------------------------
method_order <- c("MR.RGM", "MR.RGM+", "MR.RGM_NoConf", "SimpleMedian", "WeightedMedian", "IVW", "OneSampleMR", "mrbayes")
method_label <- c("MR.RGM" = "MR.RGM", "MR.RGM+" = "MR.RGM+", "MR.RGM_NoConf" = "MR.RGM_NoConf",
                  "SimpleMedian" = "MR-SimpleMedian", "WeightedMedian" = "MR-WeightedMedian", "IVW" = "MR-IVW",
                  "OneSampleMR" = "OneSampleMR", "mrbayes" = "mrbayes")
method_tex <- c("MR.RGM" = "\\texttt{MR.RGM}", "MR.RGM+" = "\\texttt{MR.RGM+}", "MR.RGM_NoConf" = "\\texttt{MR.RGM\\_NoConf}",
                "SimpleMedian" = "\\texttt{MR-SimpleMedian}", "WeightedMedian" = "\\texttt{MR-WeightedMedian}",
                "IVW" = "\\texttt{MR-IVW}", "OneSampleMR" = "\\texttt{OneSampleMR}", "mrbayes" = "\\texttt{mrbayes}")
scen_text <- c(scalefree        = "a scale-free network with feedback loops and unmeasured confounding",
               smallworld       = "a small-world network with feedback loops and unmeasured confounding",
               smallworld_pleio = "a small-world network with feedback loops, unmeasured confounding, and horizontal pleiotropy")
scen_label <- c(scalefree = "scalefree", smallworld = "smallworld", smallworld_pleio = "smallworld_horizontal")
fig_stub   <- c(scalefree = "scalefree", smallworld = "smallworld", smallworld_pleio = "smallworld_pleio")
n_levels <- sort(unique(res$n))

## ---------------------------------------------------------------------------
## Summary: mean, population SD (as in the original tables), replicate count
## ---------------------------------------------------------------------------
ok  <- res[!is.na(res$value) & res$family != "diagnostic", ]
key <- with(ok, interaction(scenario, p, n, method, family, metric, drop = TRUE, sep = "|"))
summ <- do.call(rbind, lapply(split(ok, key), function(d) {
  m <- nrow(d)
  data.frame(d[1, c("scenario", "p", "n", "method", "family", "metric")],
             mean = mean(d$value), sd = if (m > 1) sd(d$value) * sqrt((m - 1) / m) else 0, reps = m)
}))
rownames(summ) <- NULL
write.csv(summ, "tables/simulation_summary.csv", row.names = FALSE)

fmt <- function(x) ifelse(abs(x) >= 0.001 | x == 0, sprintf("%.3f", x), format(signif(x, 1), scientific = FALSE))
cell <- function(s, sc, p_, n_, meth, fam, met) {
  r <- s[s$scenario == sc & s$p == p_ & s$n == n_ & s$method == meth & s$family == fam & s$metric == met, ]
  if (!nrow(r)) return("--")
  paste0(fmt(r$mean), " (", fmt(r$sd), ")", if (r$reps < N_REP) "$^\\dagger$" else "")
}
any_short <- function(s, sc, fam) any(s$reps[s$scenario == sc & s$family == fam] < N_REP)
cap_tail <- function(short) paste0(" Entries are mean (SD) over ", N_REP, " replicates",
                                   if (short) "; $^\\dagger$ fewer valid replicates (see \\texttt{simulation\\_summary.csv})" else "", ".")

## ---------------------------------------------------------------------------
## Tables
## ---------------------------------------------------------------------------
write_method_table <- function(sc, p_, fam) {
  mets <- if (fam == "graph") c("AUC", "TPR", "FDR", "MCC") else c("MaxAbsDev", "MeanAbsDev", "MeanSqDev")
  meths <- intersect(method_order, unique(summ$method[summ$scenario == sc & summ$family == fam]))
  what  <- if (fam == "graph") "Graph recovery" else "Causal effect estimation"
  lab   <- sprintf("tab:feedback_%s_%s_90_p=%d", if (fam == "graph") "graph_structure" else "causal_effect", scen_label[[sc]], p_)
  L <- c("\\begin{table}[htb]", "\\centering", "\\tiny",
         sprintf("\\caption{%s performance in %s, with network size \\(p{=}%d\\).%s}", what, scen_text[[sc]], p_,
                 cap_tail(any_short(summ, sc, fam))),
         sprintf("\\label{%s}", lab),
         sprintf("\\begin{tabular}{|c|c|%s}", paste(rep("c|", length(mets)), collapse = "")),
         "\\hline",
         paste0("\\textbf{Setting} & \\textbf{Method} & ", paste0("\\textbf{", mets, "}", collapse = " & "), " \\\\"),
         "\\hline", "\\hline")
  for (n_ in n_levels) {
    for (k in seq_along(meths)) {
      vals <- vapply(mets, function(mt) cell(summ, sc, p_, n_, meths[k], fam, mt), character(1))
      L <- c(L, paste0(if (k == 1) sprintf("$n = %s$", format(n_, big.mark = ",")) else "", " & ",
                       method_tex[[meths[k]]], " & ", paste(vals, collapse = " & "), " \\\\"))
    }
    L <- c(L, "\\hline")
  }
  L <- c(L, "\\end{tabular}", "\\end{table}")
  writeLines(L, file.path("tables", paste0(sub("^tab:", "", lab), ".tex")))
}

write_conf_table <- function(sc, meth, suffix) {
  mets <- c("AUC", "TPR", "FDR", "MCC")
  lab  <- sprintf("tab:feedback_conf_auc_%s_%s", scen_label[[sc]], suffix)
  L <- c("\\begin{table}[htb]", "\\centering", "\\small",
         sprintf("\\caption{Confounding structure recovery performance using %s in %s, across network sizes \\(p\\in\\{5,10\\}\\).%s}",
                 method_tex[[meth]], scen_text[[sc]], cap_tail(any(summ$reps[summ$scenario == sc & summ$family == "confounding" & summ$method == meth] < N_REP))),
         sprintf("\\label{%s}", lab),
         "\\begin{tabular}{|c|c|c|c|c|c|}", "\\hline",
         "\\textbf{Setting} & \\textbf{Sample Size} & \\textbf{AUC} & \\textbf{TPR} & \\textbf{FDR} & \\textbf{MCC} \\\\",
         "\\hline", "\\hline")
  for (p_ in sort(unique(summ$p))) {
    L <- c(L, sprintf("\\multirow{%d}{*}{$p = %d$}", length(n_levels), p_))
    for (n_ in n_levels) {
      vals <- vapply(mets, function(mt) cell(summ, sc, p_, n_, meth, "confounding", mt), character(1))
      L <- c(L, paste0("    & ", format(n_, big.mark = ","), " & ", paste(vals, collapse = " & "), " \\\\"))
    }
    L <- c(L, "\\hline")
  }
  L <- c(L, "\\end{tabular}", "\\end{table}")
  writeLines(L, file.path("tables", paste0(sub("^tab:", "", lab), ".tex")))
}

for (sc in unique(res$scenario)) {
  for (p_ in sort(unique(res$p[res$scenario == sc]))) for (fam in c("graph", "effect")) write_method_table(sc, p_, fam)
  write_conf_table(sc, "MR.RGM", "A")
  if (sc == "smallworld_pleio") write_conf_table(sc, "MR.RGM+", "B")
}

## ---------------------------------------------------------------------------
## Figures
## ---------------------------------------------------------------------------
n_lab <- function(x) paste0("Sample Size: ", x)
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
        legend.text  = element_text(size = 20), legend.title = element_blank(), legend.position = "top")

method_box <- function(sc, p_, fam, met, ylab, title, file) {
  d <- res[res$scenario == sc & res$p == p_ & res$family == fam & res$metric == met & !is.na(res$value), ]
  d$Method <- factor(unname(method_label[d$method]), levels = unname(method_label[intersect(method_order, unique(d$method))]))
  d$n <- factor(d$n, levels = n_levels, labels = n_lab(n_levels))
  g <- ggplot(d, aes(Method, value, fill = Method)) + geom_boxplot() +
    facet_wrap(~ n, scales = "free_x") + guides(fill = "none") + theme_box +
    labs(title = title, y = ylab, x = "Methods")
  ggsave(file.path("figures", file), g, width = 15, height = 12, units = "in")
}

single_box <- function(d, title, file, ylim_low) {
  d$n      <- factor(d$n, levels = n_levels, labels = paste0("n = ", comma(n_levels)))
  d$Method <- factor(unname(method_label[d$method]), levels = c("MR.RGM", "MR.RGM+"))
  two <- length(unique(d$method)) > 1
  lo  <- min(ylim_low, min(d$value, na.rm = TRUE))
  g <- ggplot(d, aes(n, value)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.4, colour = "grey60")
  if (two) {
    g <- g + geom_boxplot(aes(fill = Method), width = 0.6, colour = "grey15", alpha = 0.88,
                          outlier.shape = 16, outlier.size = 1.8, position = position_dodge(width = 0.75)) +
      scale_fill_manual(values = c("MR.RGM" = "#0072B2", "MR.RGM+" = "#E69F00"))
  } else {
    g <- g + geom_boxplot(width = 0.6, fill = "#0072B2", colour = "grey15", alpha = 0.88,
                          outlier.shape = 16, outlier.size = 1.8) +
      geom_jitter(width = 0.08, height = 0, size = 1.4, alpha = 0.5, colour = "grey20")
  }
  g <- g + coord_cartesian(ylim = c(lo, 1)) +
    labs(title = title, x = "Sample size", y = "Area Under Curve (AUC)") + theme_line
  ggsave(file.path("figures", file), g, width = 14, height = 8, units = "in")
}

scen_title <- c(scalefree = "scale-free", smallworld = "small-world", smallworld_pleio = "small-world + horizontal pleiotropy")
for (sc in unique(res$scenario)) for (p_ in sort(unique(res$p[res$scenario == sc]))) {
  st <- fig_stub[[sc]]
  method_box(sc, p_, "graph",  "AUC",        "Area Under Curve (AUC)",        sprintf("AUC comparison (%s, network size = %d)", scen_title[[sc]], p_), sprintf("%s_auc_p=%d.pdf", st, p_))
  method_box(sc, p_, "effect", "MeanAbsDev", "Mean Absolute Deviation (MAD)", sprintf("MAD comparison (%s, network size = %d)", scen_title[[sc]], p_), sprintf("%s_mad_p=%d.pdf", st, p_))
  cm <- if (sc == "smallworld_pleio") c("MR.RGM", "MR.RGM+") else "MR.RGM"
  dc <- res[res$scenario == sc & res$p == p_ & res$family == "confounding" & res$metric == "AUC" & res$method %in% cm & !is.na(res$value), ]
  single_box(dc, sprintf("Confounding recovery AUC (%s, p = %d)", scen_title[[sc]], p_), sprintf("%s_conf_p=%d.pdf", st, p_), 0.4)
  if (sc == "smallworld_pleio") {
    di <- res[res$scenario == sc & res$p == p_ & res$family == "iv_selection" & !is.na(res$value), ]
    single_box(di, sprintf("Instrument-trait selection AUC (MR.RGM+, %s, p = %d)", scen_title[[sc]], p_),
               sprintf("AUC_SW_MR.RGMPlus_B_p_%d.pdf", p_), 0.85)
  }
}
message("Wrote ", length(list.files("figures")), " figures and ", length(list.files("tables", pattern = "\\.tex$")), " tables.")
