###############################################################################
## GTEx_Results.R  —  tables, network figure and motifs from output/GTEx_fit.rds
##
## Writes (all generated from the fitted model; nothing typed by hand):
##   tables/gtex_mcmc_diagnostics.csv       split-half PIP agreement, acceptance rates; figures/gtex_loglik_trace.pdf
##   tables/gtex_causal_PIPs_all.csv, tables/gtex_confounding_PIPs_all.csv
##   tables/gtex_threshold_sweep.csv      edges passing each PIP threshold (+ Bayes factors)
##   tables/gtex_causal_edges.tex         reported causal edges, explicit direction
##   tables/gtex_confounding_links.tex    reported confounding links
##   tables/gtex_label_map.tex            node label -> gene symbol (Supplement)
##   tables/gtex_motifs_all.csv           every motif among reported edges + probability
##   figures/gtex_mtor_network.pdf        network figure (Main text)
##   figures/motif_gtex_{feedback,feedforward,cascade}.pdf   most probable motif of each type
##
## Display thresholds: GTEX_THR_DIR (causal) and GTEX_THR_CONF (confounding). Read them against the
## priors (causal 0.75, confounding 0.5); see tables/gtex_threshold_sweep.csv.
##
## Usage (this folder): source("GTEx_Results.R")
###############################################################################

source("../real_data_functions.R")
if (!exists("GTEX_THR_DIR"))  GTEX_THR_DIR  <- 0.85
if (!exists("GTEX_THR_CONF")) GTEX_THR_CONF <- 0.50

res <- readRDS(file.path("output", "GTEx_fit.rds"))
fit <- res$fit; labels <- res$labels
dir.create("tables", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)
message(sprintf("GTEx fit: n = %d, p = %d, k = %d, l = %d (MR.RGM %s)", res$n, res$p, res$k, res$l, res$package_version))

mcmc_diagnostics(fit, "tables", "figures", "gtex")   # convergence: split-half PIP agreement + log-likelihood trace

edge_tables(fit, labels, GTEX_THR_DIR, GTEX_THR_CONF, "tables", "gtex")

label_table(labels, res$symbols, file.path("tables", "gtex_label_map.tex"),
            caption = "GTEx analysis: node labels used in the figures and tables, and the gene symbols whose expression was analysed.",
            label = "tab:gtex_label_map")

plot_network(fit$GammaEst, fit$ZEst, labels, GTEX_THR_DIR, GTEX_THR_CONF,
             title = sprintf("GTEx skeletal muscle network (causal inc. prob. ≥ %.2f, confounding inc. prob. ≥ %.2f)", GTEX_THR_DIR, GTEX_THR_CONF),
             out_file = file.path("figures", "gtex_mtor_network.pdf"))

report_motifs(fit, labels, GTEX_THR_DIR, "tables", "figures", "gtex")
message("GTEx results written to tables/ and figures/.")
