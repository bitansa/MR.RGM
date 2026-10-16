###############################################################################
## OneK1K_Results.R  —  tables, network figures and motifs from output/OneK1K_fit.rds
##
## Writes (all generated from the fitted model; nothing typed by hand):
##   tables/onek_mcmc_diagnostics.csv       split-half PIP agreement, acceptance rates; figures/onek_loglik_trace.pdf
##   tables/onek_causal_PIPs_all.csv, tables/onek_confounding_PIPs_all.csv
##   tables/onek_threshold_sweep.csv      edges passing each PIP threshold (+ Bayes factors)
##   tables/onek_causal_edges.tex         reported causal edges, explicit direction
##   tables/onek_confounding_links.tex    reported confounding links
##   tables/onek_label_map.tex            node label -> gene symbol (Supplement)
##   tables/onek_modules.csv              genes in each displayed module
##   tables/onek_motifs_all.csv           every motif among reported edges + probability
##   figures/onek1k_network.pdf           full network
##   figures/onek1k_cluster_XX.pdf        one panel per module (Louvain communities of the
##                                        displayed network, >= 3 genes; data-driven, not hand-picked)
##   figures/motif_onek_{feedback,feedforward,cascade}.pdf   most probable motif of each type
##
## Display thresholds: ONEK_THR_DIR = 0.50 (causal; median probability model, the rule used in the
## simulations) and ONEK_THR_CONF = 0.50 (confounding; never below its prior of 0.5). Causal edges have
## prior inclusion probability 0.75, so PIPs between 0.5 and 0.75 are retained by the median
## probability model but are not strongly supported; see tables/onek_threshold_sweep.csv.
##
## Usage (this folder): source("OneK1K_Results.R")
###############################################################################

source("../real_data_functions.R")
if (!exists("ONEK_THR_DIR"))  ONEK_THR_DIR  <- 0.50   # median probability model (same rule as the simulations)
if (!exists("ONEK_THR_CONF")) ONEK_THR_CONF <- 0.50   # not below the prior inclusion probability of a confounding link

res <- readRDS(file.path("output", "OneK1K_fit.rds"))
fit <- res$fit; labels <- res$labels
dir.create("tables", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)
message(sprintf("OneK1K fit: n = %d, p = %d, k = %d, l = %d (MR.RGM %s)", res$n, res$p, res$k, res$l, res$package_version))

mcmc_diagnostics(fit, "tables", "figures", "onek")   # convergence: split-half PIP agreement + log-likelihood trace

edge_tables(fit, labels, ONEK_THR_DIR, ONEK_THR_CONF, "tables", "onek")

label_table(labels, res$symbols, file.path("tables", "onek_label_map.tex"),
            caption = "OneK1K analysis: node labels used in the figures and tables, and the gene symbols whose expression was analysed.",
            label = "tab:onek_label_map")

plot_network(fit$GammaEst, fit$ZEst, labels, ONEK_THR_DIR, ONEK_THR_CONF,
             title = sprintf("OneK1K B-cell network (causal inc. prob. ≥ %.2f, confounding inc. prob. ≥ %.2f)", ONEK_THR_DIR, ONEK_THR_CONF),
             out_file = file.path("figures", "onek1k_network.pdf"), label_size = 3.5)

mods <- network_modules(fit$GammaEst, fit$ZEst, labels, ONEK_THR_DIR, ONEK_THR_CONF, min_size = 3)
write.csv(data.frame(module = rep(seq_along(mods), lengths(mods)), gene = unlist(mods)),
          file.path("tables", "onek_modules.csv"), row.names = FALSE)
old <- list.files("figures", pattern = "^onek1k_cluster_[0-9]+\\.pdf$", full.names = TRUE)
if (length(old)) file.rename(old, sub("\\.pdf$", ".previous.pdf", old))   # never overwrite silently
for (i in seq_along(mods))
  plot_network(fit$GammaEst, fit$ZEst, labels, ONEK_THR_DIR, ONEK_THR_CONF, genes_keep = mods[[i]],
               title = sprintf("OneK1K B cells — module %d (%d genes)", i, length(mods[[i]])),
               out_file = file.path("figures", sprintf("onek1k_cluster_%02d.pdf", i)),
               width = 8, height = 6, label_size = 4.8, conf_curvature = 0.25)
message(sprintf("%d modules with >= 3 genes drawn.", length(mods)))

report_motifs(fit, labels, ONEK_THR_DIR, "tables", "figures", "onek")
message("OneK1K results written to tables/ and figures/.")
