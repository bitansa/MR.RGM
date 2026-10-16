###############################################################################
## real_data_functions.R
##
## Shared code for the GTEx and OneK1K analyses (Version 4):
##   - preprocessing helpers (centering/scaling, summary statistics)
##   - result tables written straight from the fitted model (no hand copying)
##   - network figure, data-driven modules, network motifs
##
## Conventions (MR.RGM output): matrices are [to, from], i.e. GammaEst[j, i] is the
## posterior inclusion probability (PIP) of the causal edge i -> j. ZEst is the
## symmetric matrix of confounding PIPs.
##
## Priors (package defaults, stated in the manuscript):
##   causal edge  : rho_jh ~ Beta(3, 1)  -> prior inclusion probability 0.75
##   confounding  : pi = 0.5             -> prior inclusion probability 0.50
## A PIP should therefore be read against its prior; bayes_factor() converts it.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggraph)
  library(igraph); library(tidygraph); library(scales); library(grid); library(MR.RGM)
})

## Version 4 results require MR.RGM >= 0.1.2 (corrected samplers; see NEWS.md)
if (packageVersion("MR.RGM") < "0.1.2")
  stop("MR.RGM ", packageVersion("MR.RGM"), " is installed; these scripts need >= 0.1.2: remotes::install_github(\"bitansa/MR.RGM\")")


PRIOR_CAUSAL <- 0.75
PRIOR_CONF   <- 0.50
bayes_factor <- function(pip, prior) (pip / (1 - pip)) / (prior / (1 - prior))

## ---------------------------------------------------------------------------
## Preprocessing
## ---------------------------------------------------------------------------
## Center and scale every column; drop columns with zero variance (reported).
standardize_cols <- function(M, what) {
  M <- as.matrix(M)
  s <- apply(M, 2, sd, na.rm = TRUE)
  keep <- is.finite(s) & s > 0
  if (any(!keep))
    message(sprintf("%s: dropped %d constant column(s): %s", what, sum(!keep),
                    paste(head(colnames(M)[!keep], 10), collapse = ", ")))
  M <- M[, keep, drop = FALSE]
  if (anyNA(M)) stop(what, " contains missing values; impute or remove them before standardizing.")
  Z <- scale(M, center = TRUE, scale = TRUE)
  attr(Z, "scaled:center") <- NULL; attr(Z, "scaled:scale") <- NULL
  Z
}

## Summary statistics from centered data (covariance-type matrices)
summary_stats <- function(Y, X, U = NULL) {
  n <- nrow(Y); stopifnot(nrow(X) == n, is.null(U) || nrow(U) == n)
  out <- list(n = n, Syy = crossprod(Y) / n, Syx = crossprod(Y, X) / n, Sxx = crossprod(X) / n)
  if (!is.null(U)) out <- c(out, list(Syu = crossprod(Y, U) / n, Sxu = crossprod(X, U) / n, Suu = crossprod(U) / n))
  out
}

## ---------------------------------------------------------------------------
## Result tables (B5: every reciprocal pair written as "X -> Y / Y -> X")
## ---------------------------------------------------------------------------
edge_tables <- function(fit, labels, thr_dir, thr_conf, out_dir, prefix) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  G <- fit$GammaEst; Z <- fit$ZEst; A <- fit$AEst
  p <- length(labels); dimnames(G) <- dimnames(Z) <- dimnames(A) <- list(labels, labels)

  ## All ordered pairs (from -> to)
  causal <- expand.grid(to = seq_len(p), from = seq_len(p)) %>% filter(to != from) %>%
    mutate(edge = paste(labels[from], "->", labels[to]),
           from_label = labels[from], to_label = labels[to],
           PIP = G[cbind(to, from)], PIP_reverse = G[cbind(from, to)],
           effect = A[cbind(to, from)], BF = bayes_factor(PIP, PRIOR_CAUSAL)) %>%
    arrange(desc(PIP))
  conf <- expand.grid(i = seq_len(p), j = seq_len(p)) %>% filter(i < j) %>%
    mutate(pair = paste(labels[i], "--", labels[j]), PIP = Z[cbind(i, j)],
           BF = bayes_factor(PIP, PRIOR_CONF)) %>% arrange(desc(PIP))
  write.csv(causal %>% select(edge, from = from_label, to = to_label, PIP, PIP_reverse, effect, BF),
            file.path(out_dir, paste0(prefix, "_causal_PIPs_all.csv")), row.names = FALSE)
  write.csv(conf %>% select(pair, PIP, BF), file.path(out_dir, paste0(prefix, "_confounding_PIPs_all.csv")), row.names = FALSE)

  ## How many edges pass at several thresholds (to choose display thresholds)
  thr <- c(0.50, 0.75, 0.85, 0.90, 0.95)
  sweep <- data.frame(threshold = thr,
                      causal_edges = sapply(thr, function(t) sum(causal$PIP >= t)),
                      causal_BF    = round(bayes_factor(thr, PRIOR_CAUSAL), 2),
                      confounding_links = sapply(thr, function(t) sum(conf$PIP >= t)),
                      confounding_BF    = round(bayes_factor(thr, PRIOR_CONF), 2))
  write.csv(sweep, file.path(out_dir, paste0(prefix, "_threshold_sweep.csv")), row.names = FALSE)
  message("Edges passing each threshold (BF = Bayes factor implied by that PIP given the prior):")
  print(sweep, row.names = FALSE)

  ## Reported edges: unidirectional and reciprocal pairs, with explicit direction
  sel <- causal %>% filter(PIP >= thr_dir)
  recip <- sel %>% filter(PIP_reverse >= thr_dir) %>% filter(from < to)
  uni   <- sel %>% filter(PIP_reverse <  thr_dir)
  tex_row <- function(x) gsub("_", "\\\\_", x)
  L <- c(sprintf("%% Generated by real_data_functions.R from the fitted model (%s). Causal PIP >= %.2f.", prefix, thr_dir),
         "\\begin{tabular}{|l|c|}", "\\hline", "\\textbf{Causal edge} & \\textbf{PIP} \\\\", "\\hline")
  for (r in seq_len(nrow(uni)))
    L <- c(L, sprintf("\\textit{%s} $\\rightarrow$ \\textit{%s} & %.3f \\\\", tex_row(uni$from_label[r]), tex_row(uni$to_label[r]), uni$PIP[r]))
  for (r in seq_len(nrow(recip)))
    L <- c(L, sprintf("\\textit{%s} $\\rightleftarrows$ \\textit{%s} & %.3f\\,/\\,%.3f \\\\",
                      tex_row(recip$from_label[r]), tex_row(recip$to_label[r]), recip$PIP[r], recip$PIP_reverse[r]))
  L <- c(L, "\\hline", "\\end{tabular}",
         sprintf("%% Reciprocal pairs: first PIP = %s, second = reverse direction.", "left -> right"))
  writeLines(L, file.path(out_dir, paste0(prefix, "_causal_edges.tex")))

  cs <- conf %>% filter(PIP >= thr_conf)
  L <- c(sprintf("%% Generated from the fitted model (%s). Confounding PIP >= %.2f (prior inclusion probability %.2f).", prefix, thr_conf, PRIOR_CONF),
         "\\begin{tabular}{|l|c|}", "\\hline", "\\textbf{Confounded pair} & \\textbf{PIP} \\\\", "\\hline")
  for (r in seq_len(nrow(cs))) {
    g <- strsplit(cs$pair[r], " -- ")[[1]]
    L <- c(L, sprintf("\\textit{%s} -- \\textit{%s} & %.3f \\\\", tex_row(g[1]), tex_row(g[2]), cs$PIP[r]))
  }
  L <- c(L, "\\hline", "\\end{tabular}")
  writeLines(L, file.path(out_dir, paste0(prefix, "_confounding_links.tex")))
  invisible(list(causal = causal, conf = conf, sweep = sweep))
}

## C2: node label -> gene symbol used in the data
label_table <- function(labels, symbols, out_file, caption, label) {
  esc <- function(x) gsub("_", "\\\\_", x)
  rows <- sprintf("%s & \\textit{%s} \\\\", esc(labels), esc(symbols))
  half <- ceiling(length(rows) / 2)
  left <- rows[seq_len(half)]; right <- c(rows[-seq_len(half)], rep(" & \\\\", 2 * half - length(rows)))
  L <- c("\\begin{table}[htb]", "\\centering", "\\small", sprintf("\\caption{%s}", caption), sprintf("\\label{%s}", label),
         "\\begin{tabular}{|l|l||l|l|}", "\\hline",
         "\\textbf{Node label} & \\textbf{Gene symbol} & \\textbf{Node label} & \\textbf{Gene symbol} \\\\", "\\hline",
         paste(sub(" \\\\\\\\$", "", left), "&", right), "\\hline", "\\end{tabular}", "\\end{table}")
  writeLines(L, out_file)
}

## ---------------------------------------------------------------------------
## Network figure (style of the manuscript; legend colour fixed: orange)
## ---------------------------------------------------------------------------
plot_network <- function(G, Z, labels, thr_dir, thr_conf, title, genes_keep = NULL,
                         layout_alg = "kk", seed = 1, out_file = NULL, width = 14, height = 14,
                         label_size = 5, conf_curvature = 0.45, title_size = 12, caption_size = 12) {
  set.seed(seed)
  dimnames(G) <- dimnames(Z) <- list(labels, labels)
  if (!is.null(genes_keep)) { k <- intersect(genes_keep, labels); G <- G[k, k, drop = FALSE]; Z <- Z[k, k, drop = FALSE]; labels <- k }
  di <- which(G >= thr_dir & row(G) != col(G), arr.ind = TRUE)
  ec <- tibble(from = colnames(G)[di[, 2]], to = rownames(G)[di[, 1]], prob = as.numeric(G[di]))
  if (nrow(ec)) ec <- ec %>% mutate(type = "causal", edge_w = rescale(prob, to = c(0.55, 1.55)), edge_a = rescale(prob, to = c(0.65, 1)))
  ci <- which(Z >= thr_conf & upper.tri(Z), arr.ind = TRUE)
  ef <- tibble(from = colnames(Z)[ci[, 2]], to = rownames(Z)[ci[, 1]], prob = as.numeric(Z[ci]))
  if (nrow(ef)) ef <- ef %>% mutate(type = "conf", edge_w = rescale(prob, to = c(0.35, 1.05)), edge_a = rescale(prob, to = c(0.45, 0.9)))
  edges <- bind_rows(ec, ef)
  if (!nrow(edges)) { message("No edges above the thresholds: ", title); return(invisible(NULL)) }
  lay <- create_layout(tbl_graph(nodes = tibble(name = labels), edges = edges, directed = TRUE), layout = layout_alg)
  p <- ggraph(lay) +
    geom_node_label(aes(label = name), size = label_size, label.size = 0.25,
                    label.padding = unit(0.18, "lines"), fill = "white", colour = "grey10") +
    geom_edge_arc(aes(edge_width = edge_w, edge_alpha = edge_a, filter = type == "conf"),
                  curvature = conf_curvature, colour = "#D55E00", show.legend = FALSE,
                  start_cap = circle(4, "mm"), end_cap = circle(4, "mm"), lineend = "round") +
    geom_edge_link(aes(edge_width = edge_w, edge_alpha = edge_a, filter = type == "causal"),
                   colour = "#0072B2", show.legend = FALSE,
                   arrow = arrow(length = unit(3.2, "mm"), type = "closed"),
                   start_cap = circle(4, "mm"), end_cap = circle(4, "mm"), lineend = "round", linejoin = "mitre") +
    scale_edge_width(range = c(0.25, 1.6)) + scale_edge_alpha(range = c(0.8, 1)) +
    theme_void(base_size = 12) +
    theme(plot.margin = margin(3, 3, 3, 3), plot.title = element_text(face = "bold", size = title_size),
          plot.caption = element_text(size = caption_size),
          panel.border = element_rect(colour = "grey85", fill = NA, linewidth = 0.4)) +
    labs(title = title, caption = "Blue: causal (double-headed if bidirectional); orange: confounding (curved).")
  if (!is.null(out_file)) ggsave(out_file, p, width = width, height = height, units = "in", device = grDevices::cairo_pdf)
  invisible(p)
}

## Data-driven modules (replaces hand-picked panels): Louvain communities of the
## displayed network (causal and confounding edges, undirected), size >= min_size.
network_modules <- function(G, Z, labels, thr_dir, thr_conf, min_size = 3, seed = 1) {
  set.seed(seed)
  W <- pmax((G >= thr_dir) * G, t((G >= thr_dir) * G), (Z >= thr_conf) * Z); diag(W) <- 0
  dimnames(W) <- list(labels, labels)
  g <- graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE, diag = FALSE)
  g <- delete_vertices(g, which(degree(g) == 0))
  if (vcount(g) == 0) return(list())
  cm <- membership(cluster_louvain(g))
  mods <- split(names(cm), cm)
  mods <- mods[vapply(mods, length, 1L) >= min_size]
  mods[order(-vapply(mods, length, 1L))]
}

## ---------------------------------------------------------------------------
## Network motifs: enumerate among edges with PIP >= thr_dir, score each with
## NetworkMotif() (posterior probability that all motif edges are present)
## ---------------------------------------------------------------------------
make_gamma0 <- function(labels, from, to) {
  p <- length(labels); M <- matrix(0, p, p, dimnames = list(labels, labels))
  for (e in seq_along(from)) M[to[e], from[e]] <- 1
  M
}

enumerate_motifs <- function(fit, labels, thr_dir) {
  G <- fit$GammaEst; dimnames(G) <- list(labels, labels)
  E <- (G >= thr_dir) & (row(G) != col(G))                     # E[to, from]
  has <- function(a, b) E[b, a]                                  # edge a -> b
  out <- list()
  add <- function(type, nodes, from, to) {
    prob <- NetworkMotif(make_gamma0(labels, from, to), fit$GammaPst)
    out[[length(out) + 1]] <<- data.frame(type = type, nodes = paste(nodes, collapse = ", "),
                                           edges = paste(paste(from, "->", to), collapse = "; "),
                                           probability = as.numeric(prob), stringsAsFactors = FALSE)
  }
  idx <- seq_along(labels)
  for (a in idx) for (b in idx) for (c in idx) {
    if (length(unique(c(a, b, c))) < 3) next
    A <- labels[a]; B <- labels[b]; C <- labels[c]
    if (has(A, B) && has(B, C) && has(C, A) && a < b && a < c)            # feedback loop, one rotation
      add("feedback", c(A, B, C), c(A, B, C), c(B, C, A))
    if (has(A, B) && has(A, C) && has(B, C))                              # feed-forward loop
      add("feedforward", c(A, B, C), c(A, A, B), c(B, C, C))
    if (has(A, B) && has(B, C) && !has(A, C))                             # cascade (no shortcut)
      add("cascade", c(A, B, C), c(A, B), c(B, C))
  }
  if (!length(out)) return(data.frame())
  do.call(rbind, out) %>% arrange(type, desc(probability))
}

## Motif drawing (style of the manuscript motif panels)
save_motif_pdf <- function(type, nodes, file) {
  coords <- switch(type,
    feedback    = matrix(c(0, 1, -0.87, -0.5, 0.87, -0.5), ncol = 2, byrow = TRUE),
    feedforward = matrix(c(-1, 0, 0, 1, 1, 0), ncol = 2, byrow = TRUE),
    cascade     = matrix(c(0, 1, 0, -0.2, 1, -0.2), ncol = 2, byrow = TRUE))
  from <- switch(type, feedback = nodes, feedforward = nodes[c(1, 1, 2)], cascade = nodes[1:2])
  to   <- switch(type, feedback = nodes[c(2, 3, 1)], feedforward = nodes[c(2, 3, 3)], cascade = nodes[2:3])
  g <- graph_from_data_frame(data.frame(from = from, to = to), directed = TRUE,
                             vertices = data.frame(name = nodes))
  rownames(coords) <- nodes; lay <- coords[V(g)$name, , drop = FALSE]
  pdf(file, width = 4.2, height = 4.2, useDingbats = FALSE); on.exit(dev.off(), add = TRUE)
  par(mar = c(1.6, 1.6, 1.6, 1.6), xpd = NA)
  plot(g, layout = lay, rescale = FALSE, asp = 0,
       xlim = range(lay[, 1]) + c(-0.38, 0.38), ylim = range(lay[, 2]) + c(-0.38, 0.38),
       vertex.size = 30, vertex.color = "white", vertex.frame.color = "black",
       vertex.label = V(g)$name, vertex.label.cex = 0.95, vertex.label.color = "black",
       vertex.label.family = "serif", edge.width = 2, edge.color = "black", edge.arrow.size = 0.6)
  box()
}

## Write motif table and draw the most probable motif of each type
report_motifs <- function(fit, labels, thr_dir, out_dir, fig_dir, prefix) {
  m <- enumerate_motifs(fit, labels, thr_dir)
  if (!nrow(m)) { message("No motifs among edges with PIP >= ", thr_dir); return(invisible(m)) }
  write.csv(m, file.path(out_dir, paste0(prefix, "_motifs_all.csv")), row.names = FALSE)
  top <- m %>% group_by(type) %>% slice_max(probability, n = 1, with_ties = FALSE) %>% ungroup()
  print(top)
  for (r in seq_len(nrow(top)))
    save_motif_pdf(top$type[r], strsplit(top$nodes[r], ", ")[[1]],
                   file.path(fig_dir, sprintf("motif_%s_%s.pdf", prefix, top$type[r])))
  invisible(m)
}

## ---------------------------------------------------------------------------
## MCMC diagnostics
## ---------------------------------------------------------------------------
## Causal PIPs from the first and the second half of the kept draws should agree if the
## chain has converged; also draws the log-likelihood trace.
mcmc_diagnostics <- function(fit, out_dir, fig_dir, prefix) {
  G <- fit$GammaPst; m <- dim(G)[3]; h <- floor(m / 2)
  off <- row(G[, , 1]) != col(G[, , 1])
  p1 <- apply(G[, , 1:h, drop = FALSE], c(1, 2), mean)[off]
  p2 <- apply(G[, , (h + 1):m, drop = FALSE], c(1, 2), mean)[off]
  res <- data.frame(kept_draws = m,
                    pip_correlation_halves = round(cor(p1, p2), 4),
                    pip_max_abs_diff_halves = round(max(abs(p1 - p2)), 4),
                    edges_crossing_0.5 = sum((p1 > 0.5) != (p2 > 0.5)),
                    accept_rate_A_pct = round(fit$AccptA, 1),
                    accept_rate_B_pct = if (!is.null(fit$AccptB)) round(fit$AccptB, 1) else NA)
  write.csv(res, file.path(out_dir, paste0(prefix, "_mcmc_diagnostics.csv")), row.names = FALSE)
  print(res)
  pdf(file.path(fig_dir, paste0(prefix, "_loglik_trace.pdf")), width = 8, height = 4)
  plot(fit$LLPst, type = "l", xlab = "Kept draw", ylab = "Log-likelihood", main = paste(prefix, "log-likelihood trace"))
  dev.off()
  invisible(res)
}
