###############################################################################
## OneK1K_Results.R  —  tables, network figures and motifs from output/OneK1K_fit_MRRGM.rds (or _MRRGMplus.rds)
##
## Writes (all generated programmatically from the fitted model):
##   tables/onek_mcmc_diagnostics.csv       convergence of the chains (split-R-hat, ESS, between-chain PIPs), chain halves,
##                                        acceptance; figures/onek_loglik_trace.pdf (chains overlaid), onek_chains_PIPs.pdf
##   tables/onek_causal_PIPs_all.csv, tables/onek_confounding_PIPs_all.csv
##   tables/onek_threshold_sweep.csv      edges passing each PIP threshold (+ Bayes factors)
##   tables/onek_causal_edges.tex         reported causal edges, explicit direction
##   tables/onek_confounding_links.tex    reported confounding links
##   tables/onek_label_map.tex            pathway node -> gene symbol (Supplement)
##   tables/onek_edge_identification.csv  reported edges with instrument strength and weak-identification flags
##   tables/onek_instruments.tex, onek_instrument_strength.csv, onek_instruments.csv   instruments per gene
##   tables/onek_modules.csv              genes in each displayed module
##   tables/onek_motifs_all.csv           every motif among reported edges + probability
##   figures/onek1k_network.pdf           full network
##   figures/onek1k_cluster_XX.pdf        one panel per module (Louvain communities of the
##                                        displayed network, >= 3 genes; data-driven modules)
##   figures/motif_onek_all.pdf (+ one file per type)   most probable motif of each type, with its probability
##   tables/onek_candidates.tex, onek_instruments_list.tex, onek[_plus]_all_edges.tex, realdata_diagnostics.tex,
##   onek_model_comparison.tex   paper tables (LaTeX), generated from the fitted objects
##   tables/onek_methods_comparison.tex, onek_edge_support.tex, onek_methods_edges.csv   every method on the same data
##                                        (written when output/OneK1K_competitors.rds exists, OneK1K_Competitors.R)
##
## Display thresholds: ONEK_THR_DIR = 0.50 (median probability model, the rule used in the simulations) and
## ONEK_THR_CONF = 0.50. Prior inclusion probabilities: causal 0.75, confounding 0.5 (edge width and shade in
## the figures increase with the PIP); see tables/onek_threshold_sweep.csv.
##
## Usage (this folder): source("OneK1K_Results.R")
###############################################################################

source("../real_data_functions.R")
if (!exists("ONEK_THR_DIR"))  ONEK_THR_DIR  <- 0.50   # median probability model (rule used in the simulations)
if (!exists("ONEK_THR_CONF")) ONEK_THR_CONF <- 0.50   # not below the prior inclusion probability of a confounding link

MODEL <- if (exists("ONEK_MODEL")) ONEK_MODEL else "MR.RGM"    # "MR.RGM" (main) or "MR.RGM+" (sensitivity)
stopifnot(MODEL %in% c("MR.RGM", "MR.RGM+"))
pre <- if (MODEL == "MR.RGM") "onek" else "onek_plus"           # file prefix for tables
fig <- if (MODEL == "MR.RGM") "onek1k" else "onek1k_plus"       # file prefix for network figures
OUT_DIR <- if (exists("ONEK_OUT_DIR")) ONEK_OUT_DIR else "."    # folder holding output/ (fits); tables/ and figures/ are written there
res <- readRDS(file.path(OUT_DIR, "output", if (MODEL == "MR.RGM") "OneK1K_fit_MRRGM.rds" else "OneK1K_fit_MRRGMplus.rds"))
fit <- res$fit; labels <- res$labels
old_wd <- setwd(normalizePath(OUT_DIR))                         # every path below is relative to OUT_DIR; restored at the end
dir.create("tables", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)
message(sprintf("OneK1K fit: n = %d, p = %d, k = %d, l = %d (MR.RGM %s)", res$n, res$p, res$k, res$l, res$package_version))

mcmc_diagnostics(fit, "tables", "figures", pre, convergence = res$convergence)   # convergence of the chains, package diagnostics
## The pooled posterior is interpreted only if the chains converged. Otherwise the network, modules and motifs are not
## drawn; the PIPs of every chain are written instead (tables/<prefix>_chains_PIPs.csv), to show the disagreement.
converged <- is.null(res$convergence) || isTRUE(res$convergence$converged)
if (!converged) {
  cp <- res$convergence$pip_chains; off_ <- which(row(cp[[1]]) != col(cp[[1]]), arr.ind = TRUE)
  chains_tab <- data.frame(from = labels[off_[, 2]], to = labels[off_[, 1]], sapply(cp, function(g) round(g[off_], 3)))
  names(chains_tab)[-(1:2)] <- sprintf("PIP_chain_%d", seq_along(cp))
  write.csv(chains_tab, file.path("tables", paste0(pre, "_chains_PIPs.csv")), row.names = FALSE)
  message(sprintf("%s: the %d chains did not converge (largest split-R-hat %.2f, smallest between-chain PIP correlation %.2f): network, modules and motifs are not drawn; per-chain PIPs in tables/%s_chains_PIPs.csv.",
                  MODEL, res$convergence$chains, res$convergence$rhat_edges_max, res$convergence$pip_pairwise_min_cor, pre))
}

edge_tables(fit, labels, ONEK_THR_DIR, ONEK_THR_CONF, "tables", pre)
edge_identification(fit, labels, res$instrument_strength, ONEK_THR_DIR, "tables", pre)

label_table(res$pathway_node, res$symbols, file.path("tables", "onek_label_map.tex"),
            caption = "OneK1K analysis: genes analyzed (labels in the figures and tables) and the KEGG hsa04662 node (orthology group) each belongs to.",
            label = "tab:onek_label_map")

snp_df <- data.frame(gene = rep(names(res$top_snps), lengths(res$top_snps)), variant = unlist(res$top_snps), row.names = NULL)
instrument_tables(res$instrument_strength, snp_df, res$symbols, file.path("tables", "onek"),
                  file.path("tables", "onek_instruments.tex"),
                  caption = sprintf("OneK$1$K B cells: instruments of the analyzed genes ($n = %d$; cis window $\\pm$%s). Partial $R^2$ of a gene's own SNPs given the covariates; concentration parameter $\\mu^2 = nR^2/(1-R^2)$.", res$n, sub("^\\+/- *", "", res$settings$cis_window)),
                  label = "tab:onek_instruments")

if (converged) {
plot_network(fit$GammaEst, fit$ZEst, labels, ONEK_THR_DIR, ONEK_THR_CONF,
             title = sprintf("OneK1K B cells, %s (causal inc. prob. > %.2f, confounding inc. prob. > %.2f)", MODEL, ONEK_THR_DIR, ONEK_THR_CONF),
             out_file = file.path("figures", paste0(fig, "_network.pdf")))

mods <- network_modules(fit$GammaEst, fit$ZEst, labels, ONEK_THR_DIR, ONEK_THR_CONF, min_size = 3)
write.csv(data.frame(module = rep(seq_along(mods), lengths(mods)), gene = unlist(mods)),
          file.path("tables", paste0(pre, "_modules.csv")), row.names = FALSE)
old <- list.files("figures", pattern = paste0("^", fig, "_cluster_[0-9]+\\.pdf$"), full.names = TRUE)
if (length(old)) unlink(old)                                                # existing panels, replaced below
for (i in seq_along(mods))
  plot_network(fit$GammaEst, fit$ZEst, labels, ONEK_THR_DIR, ONEK_THR_CONF, genes_keep = mods[[i]],
               title = sprintf("OneK1K B cells, %s — module %d (%d genes)", MODEL, i, length(mods[[i]])),
               out_file = file.path("figures", sprintf("%s_cluster_%02d.pdf", fig, i)),
               width = 8, height = 6.5, conf_curvature = 0.25)
message(sprintf("%d modules with >= 3 genes drawn.", length(mods)))

report_motifs(fit, labels, ONEK_THR_DIR, "tables", "figures", pre)
} else {                                                         # outputs of an earlier converged run are removed
  unlink(c(file.path("figures", paste0(fig, "_network.pdf")), list.files("figures", paste0("^", fig, "_cluster_[0-9]+\\.pdf$"), full.names = TRUE),
           list.files("figures", paste0("^motif_", pre, "_(all|cascade|feedback|feedforward)\\.pdf$"), full.names = TRUE), file.path("tables", paste0(pre, c("_modules.csv", "_motifs_all.csv")))))
}

## Paper tables written from the fitted object(s)
DATA_NAME <- "OneK$1$K B cells"
candidates_tex(res, file.path("tables", "onek_candidates.tex"), DATA_NAME)
instruments_list_tex(res, file.path("tables", "onek_instruments_list.tex"), DATA_NAME)
all_edges_tex(fit, labels, ONEK_THR_DIR, file.path("tables", paste0(pre, "_all_edges.tex")), DATA_NAME, MODEL)

## MR.RGM vs MR.RGM+ side by side (written once both fits exist); a model whose chains did not converge is marked
plus_note <- ""
f_main <- file.path("output", "OneK1K_fit_MRRGM.rds"); f_plus <- file.path("output", "OneK1K_fit_MRRGMplus.rds")
r_main <- if (file.exists(f_main)) readRDS(f_main); r_plus <- if (file.exists(f_plus)) readRDS(f_plus)
if (!is.null(r_plus) && !is.null(r_plus$convergence) && !isTRUE(r_plus$convergence$converged))
  plus_note <- " The chains of \\texttt{MR.RGM+} did not converge (Table~\\ref{tab:realdata_diagnostics}); its pooled PIPs are shown for reference only."
if (!is.null(r_main) && !is.null(r_plus) && identical(r_main$symbols, r_plus$symbols)) {
  compare_models(r_main$fit, r_plus$fit, labels, ONEK_THR_DIR, ONEK_THR_CONF, "tables", "onek")
  diagnostics_tex(list(r_main$fit, r_plus$fit), c("MR.RGM", "MR.RGM+"), DATA_NAME, file.path("tables", "realdata_diagnostics.tex"), r_main$settings$mcmc, ONEK_THR_DIR,
                  conv = list(r_main$convergence, r_plus$convergence))
  model_comparison_tex(r_main$fit, r_plus$fit, labels, ONEK_THR_DIR, file.path("tables", "onek_model_comparison.tex"), DATA_NAME, note = plus_note)
} else {
  diagnostics_tex(list(fit), MODEL, DATA_NAME, file.path("tables", paste0(pre, "_diagnostics.tex")), res$settings$mcmc, ONEK_THR_DIR,
                  conv = list(res$convergence))
}
## Every method on the same data (OneK1K_Competitors.R): agreement with MR.RGM
f_comp <- file.path("output", "OneK1K_competitors.rds")
if (!is.null(r_main) && file.exists(f_comp)) {
  comp <- readRDS(f_comp)
  if (identical(comp$labels, r_main$labels))
    methods_comparison_tables(r_main$fit, if (!is.null(r_plus) && identical(r_main$symbols, r_plus$symbols)) r_plus$fit,
                              comp, r_main$labels, ONEK_THR_DIR, "tables", "onek", DATA_NAME, note = plus_note)
  else warning("OneK1K_competitors.rds does not match the genes of the MR.RGM fit; rerun OneK1K_Competitors.R")
}
message("OneK1K results written to ", file.path(normalizePath("."), "tables"), " and figures/.")
setwd(old_wd)
