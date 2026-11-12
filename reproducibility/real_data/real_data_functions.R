###############################################################################
## real_data_functions.R
##
## Shared code for the real-data analyses (OneK1K; the helpers are generic), MR.RGM >= 0.1.4:
##   - preprocessing helpers (centering/scaling, summary statistics)
##   - result tables generated programmatically from the fitted model
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

## These scripts require MR.RGM >= 0.1.4 (MR.RGM = known map without selection of the instrument effects; MR.RGM+ =
## selection with a point mass at zero; see NEWS.md)
if (packageVersion("MR.RGM") < "0.1.4")
  stop("MR.RGM ", packageVersion("MR.RGM"), " is installed; these scripts need MR.RGM >= 0.1.4: remotes::install_github(\"bitansa/MR.RGM\")")


PRIOR_CAUSAL <- 0.75
PRIOR_CONF   <- 0.50
bayes_factor <- function(pip, prior) (pip / (1 - pip)) / (prior / (1 - prior))

## ---------------------------------------------------------------------------
## Displayed text (tables, captions, figure labels). The result files, data frames and lists keep the method keys
## "MR.RGM", "MR.RGM+", "MR.RGM_NoConf" and the scenario keys (e.g. "smallworld_pleio"); only the text that is shown is
## converted: the methods are displayed as IV-RGM, IV-RGM+ and IV-RGM_NoConf, (horizontal) pleiotropy as invalid
## instruments, and trait(s) as response(s). References to the R package MR.RGM ("package MR.RGM", "MR.RGM >= 0.1.3",
## "MR.RGM::", ...) are kept. Character vector in, character vector out (names kept); a second application changes nothing.
## sim_functions.R defines the same function, which also renames SNPs as instruments for the simulated data.
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
  word_map <- c("Traits" = "Responses", "traits" = "responses", "Trait" = "Response", "trait" = "response")
  for (i in seq_along(word_map)) x <- gsub(paste0("\\b", names(word_map)[i], "\\b"), word_map[[i]], x, perl = TRUE)
  for (i in seq_along(keep)) x <- gsub(token[i], keep[i], x, fixed = TRUE)
  names(x) <- nm
  x
}

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
## Result tables (every reciprocal pair written as "X -> Y / Y -> X")
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
  sel <- causal %>% filter(PIP > thr_dir)
  recip <- sel %>% filter(PIP_reverse > thr_dir) %>% filter(from < to)
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
  writeLines(display_text(L), file.path(out_dir, paste0(prefix, "_causal_edges.tex")))

  cs <- conf %>% filter(PIP > thr_conf)
  L <- c(sprintf("%% Generated from the fitted model (%s). Confounding PIP >= %.2f (prior inclusion probability %.2f).", prefix, thr_conf, PRIOR_CONF),
         "\\begin{tabular}{|l|c|}", "\\hline", "\\textbf{Confounded pair} & \\textbf{PIP} \\\\", "\\hline")
  for (r in seq_len(nrow(cs))) {
    g <- strsplit(cs$pair[r], " -- ")[[1]]
    L <- c(L, sprintf("\\textit{%s} -- \\textit{%s} & %.3f \\\\", tex_row(g[1]), tex_row(g[2]), cs$PIP[r]))
  }
  L <- c(L, "\\hline", "\\end{tabular}")
  writeLines(display_text(L), file.path(out_dir, paste0(prefix, "_confounding_links.tex")))
  invisible(list(causal = causal, conf = conf, sweep = sweep))
}

## Node label -> gene symbol used in the data
label_table <- function(labels, symbols, out_file, caption, label) {
  esc <- function(x) gsub("_", "\\\\_", x)
  rows <- sprintf("%s & \\textit{%s} \\\\", esc(labels), esc(symbols))
  half <- ceiling(length(rows) / 2)
  left <- rows[seq_len(half)]; right <- c(rows[-seq_len(half)], rep(" & \\\\", 2 * half - length(rows)))
  L <- c("\\begin{table}[htb]", "\\centering", "\\small", sprintf("\\caption{%s}", caption), sprintf("\\label{%s}", label),
         "\\begin{tabular}{|l|l||l|l|}", "\\hline",
         "\\textbf{Pathway node} & \\textbf{Gene symbol} & \\textbf{Pathway node} & \\textbf{Gene symbol} \\\\", "\\hline",
         paste(sub(" \\\\\\\\$", "", left), "&", right), "\\hline", "\\end{tabular}", "\\end{table}")
  writeLines(display_text(L), out_file)
}

## How well each reported edge is identified. An edge i -> j is identified through the instruments of
## the exposure i, so every edge with PIP > thr_dir is listed with the instrument strength (mu2) of
## both genes, their error variances, and a flag for opposite-sign 2-cycles (i -> j and j -> i both
## reported with effects of opposite sign), the pattern that weak instruments produce in simulations
## (reproducibility/package_check/check_design_matched.R).
edge_identification <- function(fit, labels, strength, thr_dir, out_dir, prefix, strong_mu2 = 50) {
  G <- fit$GammaEst; A <- fit$AEst; s2 <- diag(fit$SigmaEst)
  mu <- strength$mu2[match(labels, strength$gene)]
  idx <- which(G > thr_dir & row(G) != col(G), arr.ind = TRUE)            # [to, from]
  if (!nrow(idx)) { message("No causal edge with PIP > ", thr_dir, "."); return(invisible(NULL)) }
  to <- idx[, 1]; from <- idx[, 2]
  recip <- G[cbind(from, to)] > thr_dir
  out <- data.frame(edge = paste(labels[from], "->", labels[to]), PIP = round(G[idx], 3), effect = round(A[idx], 3),
                    mu2_from = mu[from], mu2_to = mu[to], error_var_from = round(s2[from], 2), error_var_to = round(s2[to], 2),
                    exposure_strong = mu[from] >= strong_mu2,
                    opposite_sign_2cycle = recip & sign(A[idx]) != sign(A[cbind(from, to)]))
  out <- out[order(-out$PIP), ]
  write.csv(out, file.path(out_dir, paste0(prefix, "_edge_identification.csv")), row.names = FALSE)
  message(sprintf("Reported edges: %d; out of a strongly instrumented gene (mu2 >= %g): %d; in an opposite-sign 2-cycle: %d",
                  nrow(out), strong_mu2, sum(out$exposure_strong), sum(out$opposite_sign_2cycle)))
  invisible(out)
}

## Instruments of the analysed genes (data description in the manuscript / Supplement)
##   strength: data frame with gene, node, instruments, partial_R2, mu2 (+ optional extra columns)
##   snps:     data frame with one row per instrument (gene, variant, ...)
instrument_tables <- function(strength, snps, genes, out_prefix, tex_file, caption, label) {
  st <- strength[match(genes, strength$gene), , drop = FALSE]
  write.csv(st, paste0(out_prefix, "_instrument_strength.csv"), row.names = FALSE)
  write.csv(snps[snps$gene %in% genes, , drop = FALSE], paste0(out_prefix, "_instruments.csv"), row.names = FALSE)
  esc <- function(x) gsub("_", "\\\\_", x)
  rows <- sprintf("\\textit{%s} & %s & %d & %.3f & %.1f \\\\", esc(st$gene), esc(st$node), as.integer(st$instruments), st$partial_R2, st$mu2)
  L <- c("\\begin{table}[htb]", "\\centering", "\\small", sprintf("\\caption{%s}", caption), sprintf("\\label{%s}", label),
         "\\begin{tabular}{|l|l|r|r|r|}", "\\hline",
         "\\textbf{Gene} & \\textbf{Pathway node} & \\textbf{Instruments} & \\textbf{Partial $R^2$} & \\textbf{$\\mu^2$} \\\\", "\\hline",
         rows, "\\hline", "\\end{tabular}", "\\end{table}")
  writeLines(display_text(L), tex_file)
}

## ---------------------------------------------------------------------------
## Network figure (style of the manuscript; confounding links drawn in orange)
## ---------------------------------------------------------------------------
plot_network <- function(G, Z, labels, thr_dir, thr_conf, title, genes_keep = NULL,
                         layout_alg = "kk", seed = 1, out_file = NULL, width = 9, height = 8,
                         label_size = 6.5, conf_curvature = 0.45, title_size = 15, caption_size = 11,
                         solid_from = 0, recip_curvature = 0.22) {   # causal edges below solid_from are dashed (0: all solid);
                                                               # the two directions of a 2-cycle are drawn as arcs of this curvature
  thesis <- exists("THESIS_STYLE")                              # print-size style, set by a figure driver: 9 pt text, labels 2.5
  if (thesis) { label_size <- 2.5; title_size <- 9; caption_size <- 7 }
  arrow_mm <- if (thesis) 2 else 3.2
  set.seed(seed)
  dimnames(G) <- dimnames(Z) <- list(labels, labels)
  if (!is.null(genes_keep)) { k <- intersect(genes_keep, labels); G <- G[k, k, drop = FALSE]; Z <- Z[k, k, drop = FALSE]; labels <- k }
  di <- which(G > thr_dir & row(G) != col(G), arr.ind = TRUE)
  ec <- tibble(from = colnames(G)[di[, 2]], to = rownames(G)[di[, 1]], prob = as.numeric(G[di]))
  ## widths and shades on a fixed scale (PIP from the threshold to 1), so that they compare across figures and do not
  ## exaggerate small differences between moderate PIPs
  if (nrow(ec)) ec <- ec %>% mutate(type = "causal", edge_w = rescale(prob, to = c(0.55, 1.55), from = c(thr_dir, 1)),
                                  edge_a = rescale(prob, to = c(0.65, 1), from = c(thr_dir, 1)),
                                  lt = ifelse(prob >= solid_from, "solid", "dashed"),
                                  recip = paste(to, from) %in% paste(from, to))     # both directions selected: a 2-cycle
  ci <- which(Z > thr_conf & upper.tri(Z), arr.ind = TRUE)
  ef <- tibble(from = colnames(Z)[ci[, 2]], to = rownames(Z)[ci[, 1]], prob = as.numeric(Z[ci]))
  if (nrow(ef)) ef <- ef %>% mutate(type = "conf", edge_w = rescale(prob, to = c(0.35, 1.05), from = c(thr_conf, 1)),
                                  edge_a = rescale(prob, to = c(0.45, 0.9), from = c(thr_conf, 1)), lt = "solid",
                                  recip = FALSE)
  edges <- bind_rows(ec, ef)
  if (!nrow(edges)) { message("No edges above the thresholds: ", title); return(invisible(NULL)) }
  lay <- create_layout(tbl_graph(nodes = tibble(name = labels), edges = edges, directed = TRUE), layout = layout_alg)
  fs_pt <- label_size * .pt; pad <- margin(2.5, 2.5, 2.5, 2.5, "mm")   # arrows stop at the label boxes
  if (thesis) pad <- margin(1.2, 1.2, 1.2, 1.2, "mm")
  p <- ggraph(lay) +
    geom_node_label(aes(label = name), size = label_size,
                    label.padding = unit(0.18, "lines"), fill = "white", colour = "grey10") +
    geom_edge_arc(aes(edge_width = edge_w, edge_alpha = edge_a, filter = type == "conf",
                      start_cap = label_rect(node1.name, fontsize = fs_pt, padding = pad),
                      end_cap = label_rect(node2.name, fontsize = fs_pt, padding = pad)),
                  strength = conf_curvature, colour = "#D55E00", show.legend = FALSE,
                  lineend = "round") +
    geom_edge_link(aes(edge_width = edge_w, edge_alpha = edge_a, edge_linetype = lt, filter = type == "causal" & !recip,
                       start_cap = label_rect(node1.name, fontsize = fs_pt, padding = pad),
                       end_cap = label_rect(node2.name, fontsize = fs_pt, padding = pad)),
                   colour = "#0072B2", show.legend = FALSE,
                   arrow = arrow(length = unit(arrow_mm, "mm"), type = "closed"),
                   lineend = "round", linejoin = "mitre") +
    ## 2-cycles: each direction as an arc; an arc bends to the same side of its direction of travel, so the two
    ## directions of a pair separate and both arrowheads are visible
    geom_edge_arc(aes(edge_width = edge_w, edge_alpha = edge_a, edge_linetype = lt, filter = type == "causal" & recip,
                      start_cap = label_rect(node1.name, fontsize = fs_pt, padding = pad),
                      end_cap = label_rect(node2.name, fontsize = fs_pt, padding = pad)),
                  strength = recip_curvature, colour = "#0072B2", show.legend = FALSE,
                  arrow = arrow(length = unit(arrow_mm, "mm"), type = "closed"),
                  lineend = "round", linejoin = "mitre") +
    scale_edge_width(range = c(0.25, 1.6), limits = c(0.35, 1.55)) + scale_edge_alpha(range = c(0.8, 1), limits = c(0.45, 1)) +
    scale_x_continuous(expand = expansion(mult = 0.14)) + scale_y_continuous(expand = expansion(mult = 0.08)) +   # keep labels inside
    scale_edge_linetype_manual(values = c(solid = "solid", dashed = "dashed"), guide = "none") +
    theme_void(base_size = if (thesis) 9 else 12) +
    theme(plot.margin = margin(3, 3, 3, 3), plot.title = element_text(face = "bold", size = title_size),
          plot.caption = element_text(size = caption_size),
          panel.border = element_rect(colour = "grey85", fill = NA, linewidth = 0.4)) +
    labs(title = display_text(title), caption = sprintf("Blue: causal (a 2-cycle as two curved arrows); orange: confounding (curved). Width and shade increase with the\nposterior inclusion probability, on a fixed scale from the threshold to 1 (prior inclusion probability: causal %.2f, confounding %.2f).", PRIOR_CAUSAL, PRIOR_CONF))
  if (thesis && exists("DROP_TITLES") && isTRUE(DROP_TITLES)) p <- p + labs(title = NULL, caption = NULL)   # the LaTeX caption carries both
  if (!is.null(out_file)) {
    if (thesis) THESIS_STYLE$save(p, out_file, width, height)
    else ggsave(out_file, p, width = width, height = height, units = "in", device = grDevices::cairo_pdf)
  }
  invisible(p)
}

## Data-driven modules: Louvain communities of the
## displayed network (causal and confounding edges, undirected), size >= min_size.
network_modules <- function(G, Z, labels, thr_dir, thr_conf, min_size = 3, seed = 1) {
  set.seed(seed)
  W <- pmax((G > thr_dir) * G, t((G > thr_dir) * G), (Z > thr_conf) * Z); diag(W) <- 0
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
## Network motifs: enumerate among edges with PIP > thr_dir, score each with
## NetworkMotif() (posterior probability that all motif edges are present)
## ---------------------------------------------------------------------------
make_gamma0 <- function(labels, from, to) {
  p <- length(labels); M <- matrix(0, p, p, dimnames = list(labels, labels))
  for (e in seq_along(from)) M[to[e], from[e]] <- 1
  M
}

enumerate_motifs <- function(fit, labels, thr_dir) {
  G <- fit$GammaEst; dimnames(G) <- list(labels, labels)
  E <- (G > thr_dir) & (row(G) != col(G))                     # E[to, from]
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
## Most probable motif of each type, drawn in the style of the network figures, with its posterior probability
motif_layout <- function(type, nodes) {
  xy <- switch(type,
    feedback    = rbind(c(0, 0.85), c(0.85, -0.55), c(-0.85, -0.55)),
    feedforward = rbind(c(-0.85, -0.55), c(0, 0.85), c(0.85, -0.55)),
    cascade     = rbind(c(-1.1, 0.15), c(0, 0.15), c(1.1, 0.15)))
  from <- switch(type, feedback = nodes, feedforward = nodes[c(1, 1, 2)], cascade = nodes[1:2])
  to   <- switch(type, feedback = nodes[c(2, 3, 1)], feedforward = nodes[c(2, 3, 3)], cascade = nodes[2:3])
  list(nodes = data.frame(name = nodes, x = xy[, 1], y = xy[, 2]), edges = data.frame(from = from, to = to))
}
plot_motifs <- function(top, file) {
  thesis <- exists("THESIS_STYLE")                              # print-size style: 9 pt text, labels 2.5
  type_lab <- c(cascade = "Cascade", feedforward = "Feed-forward loop", feedback = "Feedback loop")
  top <- as.data.frame(top)[order(match(top$type, names(type_lab))), , drop = FALSE]
  N <- E <- list()
  for (r in seq_len(nrow(top))) {
    L <- motif_layout(top$type[r], strsplit(top$nodes[r], ", ")[[1]])
    panel <- sprintf(if (thesis) "%s\n(probability %.2f)" else "%s (probability %.2f)", type_lab[[top$type[r]]], top$probability[r])
    p1 <- L$nodes[match(L$edges$from, L$nodes$name), c("x", "y")]
    p2 <- L$nodes[match(L$edges$to,   L$nodes$name), c("x", "y")]
    d <- p2 - p1; len <- sqrt(rowSums(d^2)); sh <- 0.4                   # stop arrows short of the labels
    N[[r]] <- cbind(L$nodes, panel = panel)
    E[[r]] <- data.frame(x = p1$x + d$x * sh / len, y = p1$y + d$y * sh / len,
                         xend = p2$x - d$x * sh / len, yend = p2$y - d$y * sh / len, panel = panel)
  }
  N <- do.call(rbind, N); E <- do.call(rbind, E)
  lev <- unique(N$panel); N$panel <- factor(N$panel, lev); E$panel <- factor(E$panel, lev)
  g <- ggplot() +
    geom_segment(data = E, aes(x = x, y = y, xend = xend, yend = yend), colour = "#0072B2", linewidth = if (thesis) 0.6 else 1,
                 arrow = arrow(length = unit(if (thesis) 2 else 3.2, "mm"), type = "closed")) +
    geom_label(data = N, aes(x = x, y = y, label = name), size = if (thesis) 2.5 else 5.5, fill = "white", colour = "grey10",
               label.padding = unit(if (thesis) 0.2 else 0.3, "lines")) +
    facet_wrap(~ panel, nrow = 1) + coord_equal(xlim = c(-1.55, 1.55), ylim = c(-0.95, 1.2), clip = "off") +
    theme_void(base_size = if (thesis) 9 else 14) +
    theme(strip.text = element_text(face = "bold", size = if (thesis) 8 else 15, margin = margin(b = if (thesis) 3 else 8)),
          panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.4),
          plot.margin = margin(6, 10, 6, 10))
  if (thesis) return(THESIS_STYLE$save(g, file, 4.4 * nrow(top) + 0.4, 4.2))
  ggsave(file, g, width = 4.4 * nrow(top) + 0.4, height = 4.2, units = "in", device = grDevices::cairo_pdf)
}

## Write motif table and draw the most probable motif of each type
report_motifs <- function(fit, labels, thr_dir, out_dir, fig_dir, prefix) {
  ## delete the existing motif figures first, so that the folder holds only the motifs of this fit
  stale <- file.path(fig_dir, sprintf("motif_%s_%s.pdf", prefix, c("cascade", "feedforward", "feedback", "all")))
  unlink(stale[file.exists(stale)])
  m <- enumerate_motifs(fit, labels, thr_dir)
  csv <- file.path(out_dir, paste0(prefix, "_motifs_all.csv"))
  if (!nrow(m)) {
    message("No motifs among edges with PIP > ", thr_dir)
    write.csv(data.frame(type = character(), nodes = character(), edges = character(), probability = numeric()), csv, row.names = FALSE)
    return(invisible(m))
  }
  write.csv(m, csv, row.names = FALSE)
  top <- m %>% group_by(type) %>% slice_max(probability, n = 1, with_ties = FALSE) %>% ungroup()
  print(top)
  plot_motifs(top, file.path(fig_dir, sprintf("motif_%s_all.pdf", prefix)))          # all types in one figure
  for (r in seq_len(nrow(top)))                                                        # and one file per type
    plot_motifs(top[r, ], file.path(fig_dir, sprintf("motif_%s_%s.pdf", prefix, top$type[r])))
  invisible(m)
}

## ---------------------------------------------------------------------------
## MCMC diagnostics
## ---------------------------------------------------------------------------
## Causal PIPs from the first and the second half of the kept draws should agree if the
## chain has converged; also draws the log-likelihood trace.
## Convergence and identification diagnostics: chain halves (correlation and largest difference of the causal PIPs),
## the convergence of the chains (fit_chains()), the package diagnostics (fit$diagnostics) and the package warnings
## recorded with the fit
## Convergence and identification diagnostics of a (pooled) fit; `convergence` = fit_chains()$convergence
mcmc_diagnostics <- function(fit, out_dir, fig_dir, prefix, convergence = NULL) {
  G <- fit$GammaPst; m <- dim(G)[3]; h <- floor(m / 2)
  off <- row(G[, , 1]) != col(G[, , 1])
  p1 <- apply(G[, , 1:h, drop = FALSE], c(1, 2), mean)[off]
  p2 <- apply(G[, , (h + 1):m, drop = FALSE], c(1, 2), mean)[off]
  dg <- fit$diagnostics
  res <- data.frame(kept_draws = m,
                    pip_correlation_halves = round(cor(p1, p2), 4),
                    pip_max_abs_diff_halves = round(max(abs(p1 - p2)), 4),
                    edges_crossing_0.5 = sum((p1 > 0.5) != (p2 > 0.5)),
                    chains = if (is.null(convergence)) 1 else convergence$chains,
                    rhat_loglik = if (is.null(convergence)) NA else round(convergence$rhat_loglik, 4),
                    rhat_edges_max = if (is.null(convergence)) NA else round(convergence$rhat_edges_max, 4),
                    ess_bulk_loglik = if (is.null(convergence)) NA else round(convergence$ess_bulk_loglik),
                    ess_bulk_edges_min = if (is.null(convergence)) NA else round(convergence$ess_bulk_edges_min),
                    chains_pip_max_diff = if (is.null(convergence)) NA else round(convergence$pip_pairwise_max_diff, 4),
                    chains_pip_min_cor = if (is.null(convergence)) NA else round(convergence$pip_pairwise_min_cor, 4),
                    converged = if (is.null(convergence)) NA else convergence$converged,
                    accept_rate_A_pct = round(fit$AccptA, 1),
                    accept_rate_B_pct = if (!is.null(fit$AccptB)) round(fit$AccptB, 1) else NA,
                    min_first_stage_F = if (is.null(dg$instrumentF)) NA else round(min(dg$instrumentF, na.rm = TRUE), 2),
                    spectral_radius_A = round(max(Mod(eigen(fit$AEst, only.values = TRUE)$values)), 3),   # of the (pooled) posterior mean
                    package_warnings = if (is.null(fit$warnings) || !length(fit$warnings)) "none" else paste(sub("[:(].*", "", fit$warnings), collapse = " | "))
  write.csv(res, file.path(out_dir, paste0(prefix, "_mcmc_diagnostics.csv")), row.names = FALSE)
  if (!is.null(fit$warnings) && length(fit$warnings)) writeLines(fit$warnings, file.path(out_dir, paste0(prefix, "_package_warnings.txt")))
  print(res)
  nice <- c(onek = "OneK1K, MR.RGM", onek_plus = "OneK1K, MR.RGM+")
  name <- display_text(if (prefix %in% names(nice)) nice[[prefix]] else prefix)   # displayed in the figure titles
  f_trace <- file.path(fig_dir, paste0(prefix, "_loglik_trace.pdf"))
  if (exists("THESIS_STYLE")) {                                  # print-size style: 9 pt base, tick labels 7.5 pt, fonts embedded
    sz <- THESIS_STYLE$size(f_trace, 9, 4.5)
    grDevices::cairo_pdf(f_trace, width = sz[1], height = sz[2], pointsize = 9)
    par(mar = c(3, 3.6, if (no_title()) 0.6 else 2.6, 0.6), mgp = c(1.9, 0.5, 0), tcl = -0.3, cex.main = 1, cex.lab = 1, cex.axis = 7.5 / 9)
  } else {
  pdf(f_trace, width = 9, height = 4.5)
  par(mar = c(4.5, 5.5, 3, 1), cex.main = 1.4, cex.lab = 1.3, cex.axis = 1.15)
  }
  main <- sprintf("MCMC Log-Likelihood Trace (%s)", name)
  if (exists("THESIS_STYLE")) main <- THESIS_STYLE$wrap(main, sz[1], 9)   # two lines if it does not fit the width
  if (no_title()) main <- NULL
  if (!is.null(convergence) && convergence$chains > 1) {        # one trace per chain, overlaid
    LL <- convergence$loglik_chains
    cols <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#555555")[seq_len(ncol(LL))]
    matplot(LL, type = "l", lty = 1, col = grDevices::adjustcolor(cols, 0.7),
            xlab = "Kept Draw (per chain)", ylab = "Log-Likelihood", main = main,
            ylim = range(LL) + c(0, 0.30 * diff(range(LL))))             # headroom for the legend
    legend("top", legend = sprintf("chain %d", seq_len(ncol(LL))), lty = 1, lwd = 2, bty = "n", ncol = min(ncol(LL), 2),
           cex = 0.8, col = cols, seg.len = 1.2, x.intersp = 0.5)
  } else plot(fit$LLPst, type = "l", col = "grey25", xlab = "Kept Draw", ylab = "Log-Likelihood", main = main)
  dev.off()
  if (!is.null(convergence) && convergence$chains > 1) chain_pip_figure(convergence, file.path(fig_dir, paste0(prefix, "_chains_PIPs.pdf")), name)
  invisible(res)
}

## TRUE when a figure driver asks for figures without titles (the LaTeX captions carry them)
no_title <- function() exists("THESIS_STYLE") && exists("DROP_TITLES") && isTRUE(DROP_TITLES)

## Agreement of the chains: causal PIP of every ordered pair in chains 2, 3, ... against chain 1 (points on the diagonal:
## the chains agree; scattered points: they do not)
chain_pip_figure <- function(convergence, file, name) {
  P <- convergence$pip_chains; off <- row(P[[1]]) != col(P[[1]]); K <- length(P)
  cols <- c("#E69F00", "#009E73", "#CC79A7", "#555555")[seq_len(K - 1)]
  if (exists("THESIS_STYLE")) {
    sz <- THESIS_STYLE$size(file, 4.5, 4.5)
    grDevices::cairo_pdf(file, width = sz[1], height = sz[2], pointsize = 9)
    par(mar = c(3, 3.2, if (no_title()) 0.6 else 2.6, 0.6), mgp = c(1.9, 0.5, 0), tcl = -0.3, cex.main = 1, cex.axis = 7.5 / 9)
  } else {
    pdf(file, width = 5.5, height = 5.5); par(mar = c(4.5, 4.8, 3, 1), cex.main = 1.2, cex.lab = 1.2)
  }
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), asp = 1, xlab = "Causal PIP, chain 1", ylab = "Causal PIP, other chains",
       main = if (no_title()) NULL else sprintf("Agreement of the chains (%s)", name))
  abline(0, 1, col = "grey60", lty = 2)
  for (k in 2:K) points(P[[1]][off], P[[k]][off], pch = 16, cex = 0.8, col = grDevices::adjustcolor(cols[k - 1], 0.75))
  legend("topleft", legend = sprintf("chain %d", 2:K), pch = 16, col = cols, bty = "n", cex = 0.8)
  legend("bottomright", legend = sprintf("smallest correlation %.2f", convergence$pip_pairwise_min_cor), bty = "n", cex = 0.8)
  dev.off()
}

## ---------------------------------------------------------------------------
## MR.RGM vs MR.RGM+ on the same genes and instruments
## ---------------------------------------------------------------------------
## Every causal edge reported by either model (PIP > thr_dir) and every confounding link reported by
## either model (PIP > thr_conf), with both PIPs side by side, plus the two models' diagnostics.
compare_models <- function(fit_main, fit_plus, labels, thr_dir, thr_conf, out_dir, prefix) {
  G1 <- fit_main$GammaEst; G2 <- fit_plus$GammaEst; dimnames(G1) <- dimnames(G2) <- list(labels, labels)
  Z1 <- fit_main$ZEst;     Z2 <- fit_plus$ZEst;     dimnames(Z1) <- dimnames(Z2) <- list(labels, labels)
  idx <- which((G1 > thr_dir | G2 > thr_dir) & row(G1) != col(G1), arr.ind = TRUE)
  edges <- data.frame(from = labels[idx[, 2]], to = labels[idx[, 1]],
                      PIP_MRRGM = round(G1[idx], 3), PIP_MRRGMplus = round(G2[idx], 3))
  edges$status <- ifelse(edges$PIP_MRRGM > thr_dir & edges$PIP_MRRGMplus > thr_dir, "both",
                         ifelse(edges$PIP_MRRGM > thr_dir, "MR.RGM only", "MR.RGM+ only"))
  edges <- edges[order(edges$status, -edges$PIP_MRRGM), ]
  cidx <- which((Z1 > thr_conf | Z2 > thr_conf) & upper.tri(Z1), arr.ind = TRUE)
  links <- data.frame(gene1 = labels[cidx[, 1]], gene2 = labels[cidx[, 2]],
                      PIP_MRRGM = round(Z1[cidx], 3), PIP_MRRGMplus = round(Z2[cidx], 3))
  v1 <- diag(fit_main$SigmaEst); v2 <- diag(fit_plus$SigmaEst)
  diag_tab <- data.frame(model = c("MR.RGM", "MR.RGM+"),
                         error_variance_median = round(c(median(v1), median(v2)), 2),
                         error_variance_max = round(c(max(v1), max(v2)), 2),
                         causal_PIP_median = round(c(median(G1[row(G1) != col(G1)]), median(G2[row(G2) != col(G2)])), 3),
                         causal_edges_reported = c(sum(G1 > thr_dir & row(G1) != col(G1)), sum(G2 > thr_dir & row(G2) != col(G2))),
                         confounding_PIP_range = c(paste(round(range(Z1[upper.tri(Z1)]), 2), collapse = "-"),
                                                   paste(round(range(Z2[upper.tri(Z2)]), 2), collapse = "-")))
  write.csv(edges, file.path(out_dir, paste0(prefix, "_model_comparison_edges.csv")), row.names = FALSE)
  write.csv(links, file.path(out_dir, paste0(prefix, "_model_comparison_confounding.csv")), row.names = FALSE)
  write.csv(diag_tab, file.path(out_dir, paste0(prefix, "_model_comparison_summary.csv")), row.names = FALSE)
  message(sprintf("Model comparison: %d causal edges reported by either model (%d both, %d MR.RGM only, %d MR.RGM+ only)",
                  nrow(edges), sum(edges$status == "both"), sum(edges$status == "MR.RGM only"), sum(edges$status == "MR.RGM+ only")))
  print(diag_tab)
  invisible(list(edges = edges, links = links, summary = diag_tab))
}


## ---------------------------------------------------------------------------
## Paper tables generated programmatically from the fitted objects
## ---------------------------------------------------------------------------
tex_esc <- function(x) gsub("_", "\\\\_", x)

## Diagnostics of one fit (the numbers quoted in the paper)
fit_diagnostics <- function(fit, thr_dir = 0.5) {
  G <- fit$GammaEst; A <- fit$AEst; p <- nrow(G); off <- row(G) != col(G)
  GP <- fit$GammaPst; m <- dim(GP)[3]; h <- floor(m / 2)
  p1 <- apply(GP[, , 1:h, drop = FALSE], c(1, 2), mean)[off]
  p2 <- apply(GP[, , (h + 1):m, drop = FALSE], c(1, 2), mean)[off]
  v  <- if (is.matrix(fit$SigmaEst)) diag(fit$SigmaEst) else as.numeric(fit$SigmaEst)
  rep <- off & G > thr_dir
  opp <- rep & t(rep) & sign(A) != sign(t(A))                    # edges in opposite-sign 2-cycles
  Z <- fit$ZEst
  list(halves_cor = cor(p1, p2), halves_maxdiff = max(abs(p1 - p2)),
       errvar_min = min(v), errvar_median = median(v), errvar_max = max(v),
       spectral_radius = max(Mod(eigen(A)$values)), max_abs_a = max(abs(A[off])),
       median_pip = median(G[off]), edges = sum(rep), opposite_sign_edges = sum(opp),
       conf_min = if (!is.null(Z)) min(Z[upper.tri(Z)]) else NA, conf_max = if (!is.null(Z)) max(Z[upper.tri(Z)]) else NA)
}

## Diagnostics table for one or two fits (MR.RGM and MR.RGM+), rows = diagnostics
diagnostics_tex <- function(fits, labels_models, data_name, file, mcmc, thr_dir = 0.5, conv = NULL) {
  d <- lapply(fits, fit_diagnostics, thr_dir = thr_dir)
  cell <- function(f, x, fmt = "%.2f") sprintf(fmt, d[[f]][[x]])
  rows <- list(
    c("Halves: correlation", "halves_cor", "%.2f"), c("Halves: largest difference", "halves_maxdiff", "%.2f"),
    c("Error variance: minimum", "errvar_min", "%.2f"), c("Error variance: median", "errvar_median", "%.2f"),
    c("Error variance: maximum", "errvar_max", "%.2f"), c("Spectral radius of $\\widehat{\\mathbf{A}}$", "spectral_radius", "%.2f"),
    c("$\\max|\\widehat a_{jh}|$", "max_abs_a", "%.2f"), c("Median causal PIP", "median_pip", "%.3f"),
    c(sprintf("Edges with PIP $> %.1f$", thr_dir), "edges", "%d"), c("Edges in opposite-sign 2-cycles", "opposite_sign_edges", "%d"))
  L <- c("\\begin{table}[htb]", "\\centering", "\\small",
         sprintf("\\caption{%s: diagnostics of the fits (iterations, burn-in and thinning per chain in the first row; chains with different random seeds from the same starting point, pooled when every rank-normalized split-$\\widehat{R}$ is at most 1.01; \\citealp{vehtari2021rank}). For a fit that fails this rule, the pooled values are shown for reference only. Bulk ESS = bulk effective sample size; halves = correlation (and largest absolute difference) of the causal PIPs computed from the two halves of the chain; error variances relative to the (unit) response variances; spectral radius and $\\max|\\widehat a_{jh}|$ = largest absolute entry of the posterior mean of $\\mathbf{A}$; median causal PIP (prior $%.2f$); edges = causal edges with PIP $> %.1f$; opposite-sign 2-cycles = number of these edges that belong to a pair $j \\rightleftarrows h$ with both PIPs $> %.1f$ and effects of opposite signs; confounding PIPs = range over gene pairs.}",
                 data_name, PRIOR_CAUSAL, thr_dir, thr_dir),
         sprintf("\\label{tab:%s}", sub("\\.tex$", "", basename(file))),
         sprintf("\\begin{tabular}{|l|%s|}", paste(rep("c", length(fits)), collapse = "|")), "\\hline",
         paste0("\\textbf{Diagnostic} & ", paste(sprintf("\\texttt{%s}", labels_models), collapse = " & "), " \\\\"), "\\hline")
  ## chain length of every model, then the convergence of the chains (fit_chains()), one row each
  L <- c(L, paste0("Iterations / burn-in / thinning per chain & ", paste(vapply(fits, function(f)
    sprintf("%s / %s / %d", format(f$settings$nIter, big.mark = ",", scientific = FALSE), format(f$settings$nBurnin, big.mark = ",", scientific = FALSE),
            as.integer(f$settings$Thin)), ""), collapse = " & "), " \\\\"))
  if (!is.null(conv)) {
    cv <- function(name, fmt) paste(vapply(conv, function(x) if (is.null(x)) "--" else sprintf(fmt, x[[name]]), ""), collapse = " & ")
    ## counts are printed with a thousands separator, as elsewhere in the table
    cvn <- function(name) paste(vapply(conv, function(x) if (is.null(x)) "--" else
                       format(round(x[[name]]), big.mark = ",", scientific = FALSE), ""), collapse = " & ")
    L <- c(L, paste0("Chains (pooled draws) & ", paste(vapply(conv, function(x) if (is.null(x)) "1" else sprintf("%d (%s)", x$chains, format(x$chains * x$draws_per_chain, big.mark = ",")), ""), collapse = " & "), " \\\\"),
           paste0("Split-$\\widehat{R}$: log-likelihood & ", cv("rhat_loglik", "%.3f"), " \\\\"),
           paste0("Split-$\\widehat{R}$: largest over edge indicators & ", cv("rhat_edges_max", "%.3f"), " \\\\"),
           paste0("Bulk ESS: log-likelihood & ", cvn("ess_bulk_loglik"), " \\\\"),
           paste0("Bulk ESS: smallest over edge indicators & ", cvn("ess_bulk_edges_min"), " \\\\"),
           paste0("Between chains: largest PIP difference & ", cv("pip_pairwise_max_diff", "%.3f"), " \\\\"),
           paste0("Between chains: smallest PIP correlation & ", cv("pip_pairwise_min_cor", "%.2f"), " \\\\"),
           paste0(sprintf("Converged (split-$\\widehat{R} \\le %.2f$) & ", RHAT_MAX), paste(vapply(conv, function(x) if (is.null(x)) "--" else if (isTRUE(x$converged)) "yes" else "no", ""), collapse = " & "), " \\\\"),
           "\\hline")
  }
  for (r in rows) L <- c(L, paste0(r[1], " & ", paste(vapply(seq_along(fits), function(f) sprintf(r[3], d[[f]][[r[2]]]), ""), collapse = " & "), " \\\\"))
  L <- c(L, paste0("Confounding PIPs (range) & ", paste(vapply(seq_along(fits), function(f) sprintf("%.2f--%.2f", d[[f]]$conf_min, d[[f]]$conf_max), ""), collapse = " & "), " \\\\"),
         "\\hline", "\\end{tabular}", "\\end{table}")
  writeLines(display_text(L), file)
  invisible(d)
}

## All causal edges with PIP > thr_dir (plus optional extra rows), with effect, reverse PIP and Bayes factor
all_edges_tex <- function(fit, labels, thr_dir, file, data_name, model_name, extra_edges = NULL) {
  G <- fit$GammaEst; A <- fit$AEst; dimnames(G) <- dimnames(A) <- list(labels, labels)
  idx <- which(G > thr_dir & row(G) != col(G), arr.ind = TRUE)
  e <- data.frame(from = labels[idx[, 2]], to = labels[idx[, 1]], PIP = G[idx], effect = A[idx], rev = G[cbind(idx[, 2], idx[, 1])])
  e <- e[order(-e$PIP), ]
  L <- c("\\begin{table}[htb]", "\\centering", "\\small",
         sprintf("\\caption{%s, \\texttt{%s}: all causal edges with posterior inclusion probability $> %.2f$, with the posterior mean effect (standard-deviation units), the inclusion probability of the reverse edge, and the Bayes factor implied by the PIP under the prior inclusion probability $%.2f$.%s}",
                 data_name, model_name, thr_dir, PRIOR_CAUSAL, if (is.null(extra_edges)) "" else " The rows below the line are edges listed for reference although their PIP is below the threshold."),
         sprintf("\\label{tab:%s}", sub("\\.tex$", "", basename(file))),
         "\\begin{tabular}{|l|c|c|c|c|}", "\\hline",
         "\\textbf{Edge} & \\textbf{PIP} & \\textbf{Effect} & \\textbf{PIP reverse} & \\textbf{Bayes factor} \\\\", "\\hline")
  row_of <- function(fr, to) sprintf("%s $\\rightarrow$ %s & %.3f & %+.3f & %.3f & %.2f \\\\", tex_esc(fr), tex_esc(to), G[to, fr], A[to, fr], G[fr, to], bayes_factor(G[to, fr], PRIOR_CAUSAL))
  for (r in seq_len(nrow(e))) L <- c(L, row_of(e$from[r], e$to[r]))
  if (!is.null(extra_edges)) { L <- c(L, "\\hline"); for (r in seq_len(nrow(extra_edges))) L <- c(L, row_of(extra_edges$from[r], extra_edges$to[r])) }
  L <- c(L, "\\hline", "\\end{tabular}", "\\end{table}")
  writeLines(display_text(L), file)
  invisible(e)
}

## cis window for captions: "+/- 1000 kb" is written as 1 Mb (whole megabases), any other value without "+/-"
cis_window_tex <- function(w) {
  kb <- suppressWarnings(as.numeric(sub("^\\+/- *([0-9.]+) *kb$", "\\1", w)))
  if (length(kb) == 1 && is.finite(kb) && kb %% 1000 == 0) sprintf("%g\\,Mb", kb / 1000) else sub("^\\+/- *", "", w)
}

## Candidate genes: every gene with at least one instrument, with the selection statistics and whether it was analysed
candidates_tex <- function(res, file, data_name) {
  st <- res$instrument_strength; it <- res$instrument_table
  st <- st[order(-st$mu2), ]
  kept <- st$gene %in% res$symbols
  L <- c("\\begin{table}[htb]", "\\centering", "\\scriptsize",
         sprintf("\\caption{%s: the %d pathway genes (KEGG %s) with at least one cis-SNP with $F > %g$, ordered by the concentration parameter $\\mu^2 = nR^2/(1-R^2)$ of their instruments. Cis-SNPs = SNPs tested in the $\\pm$%s window; strong = SNPs with $F > %g$; instruments = SNPs kept after LD clumping ($r^2 \\le %.1f$, at most %d); partial $R^2$ of the instruments given the covariates, $1-\\mathrm{RSS}_1/\\mathrm{RSS}_0$ with $\\mathrm{RSS}_0$ and $\\mathrm{RSS}_1$ the residual sums of squares of the expression on the covariates without and with the instruments; $F$ = joint first-stage $F$ of the $m$ instruments given the $l$ covariates, $(R^2/m)/\\{(1-R^2)/(n-m-l-1)\\}$, which can be below $10$ although every instrument has $F=t^2>10$; $p_{\\mathrm{Bonf}}$ = Bonferroni-adjusted $p$-value of the strongest cis-SNP over the cis-SNPs tested (candidate genes: $p_{\\mathrm{Bonf}}<0.05$). The %d analyzed genes are marked.}",
                 data_name, nrow(st), res$kegg_pathway, res$settings$f_min, cis_window_tex(res$settings$cis_window), res$settings$f_min, res$settings$ld_r2, res$settings$top_per_gene, length(res$symbols)),
         sprintf("\\label{tab:%s}", sub("\\.tex$", "", basename(file))),
         "\\begin{adjustbox}{max width=\\textwidth}", "\\begin{tabular}{|l|l|r|r|r|r|r|r|r|c|}", "\\hline",
         "\\textbf{Gene} & \\textbf{KEGG node} & \\textbf{Cis-SNPs} & \\textbf{Strong} & \\textbf{Instruments} & \\textbf{Partial $R^2$} & \\textbf{$F$} & \\textbf{$\\mu^2$} & \\textbf{$p_{\\mathrm{Bonf}}$} & \\textbf{Analyzed} \\\\", "\\hline")
  for (r in seq_len(nrow(st))) {
    g <- st$gene[r]; i <- match(g, it$gene)
    L <- c(L, sprintf("\\textit{%s} & %s & %d & %d & %d & %.3f & %.1f & %.1f & %s & %s \\\\", tex_esc(g), tex_esc(st$node[r]), it$cis_snps[i], it$strong[i],
                      st$instruments[r], st$partial_R2[r], st$F[r], st$mu2[r], format(st$best_p_bonf[r], digits = 2), if (kept[r]) "$\\checkmark$" else ""))
  }
  L <- c(L, "\\hline", "\\end{tabular}", "\\end{adjustbox}", "\\end{table}")
  writeLines(display_text(L), file)
  invisible(st)
}

## Instruments of the analysed genes (SNP identifiers) and their cross-gene associations
instruments_list_tex <- function(res, file, data_name) {
  ct <- res$cross_gene_t
  L <- c("\\begin{table}[htb]", "\\centering", "\\scriptsize",
         sprintf("\\caption{%s: instruments of the %d analyzed genes (%d SNPs). For each SNP, $t$ = $t$-statistic on its own gene given the covariates, and $\\max|t|$ other = largest absolute $t$-statistic on any other analyzed gene (a value above $%.1f$, i.e.\\ $F > %g$, indicates an association with another gene).}",
                 data_name, length(res$symbols), length(res$snps), sqrt(res$settings$f_min), res$settings$f_min),
         sprintf("\\label{tab:%s}", sub("\\.tex$", "", basename(file))),
         "\\begin{tabular}{|l|l|r|r|}", "\\hline", "\\textbf{Gene} & \\textbf{SNP} & \\textbf{$t$} & \\textbf{$\\max|t|$ other} \\\\", "\\hline")
  for (g in res$symbols) for (s in res$top_snps[[g]]) {
    other <- setdiff(res$symbols, g)
    L <- c(L, sprintf("\\textit{%s} & %s & %.1f & %.1f \\\\", tex_esc(g), tex_esc(s), ct[g, s], if (length(other)) max(abs(ct[other, s])) else NA))
  }
  L <- c(L, "\\hline", "\\end{tabular}", "\\end{table}")
  writeLines(display_text(L), file)
}

## MR.RGM vs MR.RGM+: edges with PIP > thr under either model
model_comparison_tex <- function(fit_main, fit_plus, labels, thr, file, data_name, note = "") {
  G1 <- fit_main$GammaEst; G2 <- fit_plus$GammaEst; dimnames(G1) <- dimnames(G2) <- list(labels, labels)
  idx <- which((G1 > thr | G2 > thr) & row(G1) != col(G1), arr.ind = TRUE)
  e <- data.frame(from = labels[idx[, 2]], to = labels[idx[, 1]], p1 = G1[idx], p2 = G2[idx])
  e$status <- ifelse(e$p1 > thr & e$p2 > thr, "both", ifelse(e$p1 > thr, "MR.RGM only", "MR.RGM+ only"))
  e <- e[order(match(e$status, c("both", "MR.RGM only", "MR.RGM+ only")), -e$p1), ]
  L <- c("\\begin{table}[htb]", "\\centering", "\\small",
         sprintf("\\caption{%s: causal edges with posterior inclusion probability $> %.2f$ under \\texttt{MR.RGM} or \\texttt{MR.RGM+} (same genes, SNPs and covariates).%s}", data_name, thr, note),
         sprintf("\\label{tab:%s}", sub("\\.tex$", "", basename(file))),
         "\\begin{tabular}{|l|c|c|l|}", "\\hline", "\\textbf{Edge} & \\textbf{PIP \\texttt{MR.RGM}} & \\textbf{PIP \\texttt{MR.RGM+}} & \\textbf{Reported by} \\\\", "\\hline")
  for (r in seq_len(nrow(e))) L <- c(L, sprintf("\\textit{%s} $\\rightarrow$ \\textit{%s} & %.3f & %.3f & %s \\\\", tex_esc(e$from[r]), tex_esc(e$to[r]), e$p1[r], e$p2[r], e$status[r]))
  L <- c(L, "\\hline", "\\end{tabular}", "\\end{table}")
  writeLines(display_text(L), file)
  invisible(e)
}

## ---------------------------------------------------------------------------
## Every method on the same data (OneK1K_Competitors.R): agreement with MR.RGM
##   <prefix>_methods_comparison.tex  per method: edge rule, edges reported, edges shared with MR.RGM, Jaccard index of
##                                    the two edge sets, share of the shared edges with the same sign of the effect
##   <prefix>_edge_support.tex        every edge reported by MR.RGM: its PIP (and the PIP under MR.RGM+) and, for every other
##                                    method, whether it reports the edge (+ / - = sign of its estimate; blank = not reported)
##   <prefix>_methods_edges.csv       the same for every ordered pair, every method
## Edges [to, from]; MR.RGM / MR.RGM+: PIP > thr (the reporting rule of the analysis); the other methods: their rule.
## ---------------------------------------------------------------------------
methods_comparison_tables <- function(fit_main, fit_plus, comp, labels, thr, out_dir, prefix, data_name, note = "") {
  p <- length(labels); off <- row(diag(p)) != col(diag(p))
  dn <- list(labels, labels)
  base <- list(`MR.RGM` = list(edge = fit_main$GammaEst > thr, est = fit_main$AEst, score = fit_main$GammaEst, rule = sprintf("PIP $> %.2f$", thr)))
  if (!is.null(fit_plus)) base$`MR.RGM+` <- list(edge = fit_plus$GammaEst > thr, est = fit_plus$AEst, score = fit_plus$GammaEst, rule = sprintf("PIP $> %.2f$", thr))
  all_m <- c(base, comp$methods)
  for (m in names(all_m)) if (!inherits(all_m[[m]], "failed")) for (k in c("edge", "est", "score")) dimnames(all_m[[m]][[k]]) <- dn
  ref <- all_m$`MR.RGM`$edge & off
  rows <- lapply(names(all_m), function(m) {
    x <- all_m[[m]]
    if (inherits(x, "failed")) return(data.frame(method = m, rule = "--", estimated = NA, edges = NA, shared = NA, jaccard = NA, sign_agree = NA, note = "failed"))
    e <- x$edge & off; e[is.na(e)] <- FALSE; sh <- e & ref
    sa <- if (sum(sh)) mean(sign(x$est[sh]) == sign(all_m$`MR.RGM`$est[sh]), na.rm = TRUE) else NA
    data.frame(method = m, rule = x$rule, estimated = sum(!is.na(x$score[off])), edges = sum(e), shared = sum(sh),
               jaccard = if (sum(e | ref)) sum(sh) / sum(e | ref) else NA, sign_agree = sa, note = "")
  })
  tab <- do.call(rbind, rows)
  mtex <- function(m) sprintf("\\texttt{%s}", gsub("_", "\\\\_", m))
  fmt  <- function(v, d = 2) ifelse(is.na(v), "--", formatC(v, format = "f", digits = d))
  L <- c("\\begin{table}[htb]", "\\centering", "\\small",
         sprintf("\\caption{%s: every method on the same genes, SNPs and donors (summary-level and pairwise methods on covariate-adjusted data). Ordered pairs the method could estimate (of %d), edges reported under its rule, edges shared with \\texttt{MR.RGM} (PIP $> %.2f$), Jaccard index of the two edge sets (shared edges divided by the edges in either set), and the share of the shared edges whose estimated effects have the same sign. The pairwise methods estimate total effects of one gene on another; MR2 needs an exposure/outcome split and does not apply.}", data_name, p * (p - 1), thr),
         sprintf("\\label{tab:%s_methods_comparison}", prefix), "\\begin{adjustbox}{max width=\\textwidth}",
         "\\begin{tabular}{|l|l|c|c|c|c|c|}", "\\hline",
         "\\textbf{Method} & \\textbf{Edge rule} & \\textbf{Pairs estimated} & \\textbf{Edges} & \\textbf{Shared with \\texttt{MR.RGM}} & \\textbf{Jaccard} & \\textbf{Same sign} \\\\", "\\hline")
  for (i in seq_len(nrow(tab))) L <- c(L, sprintf("%s & %s & %s & %s & %s & %s & %s \\\\", mtex(tab$method[i]), tab$rule[i],
    fmt(tab$estimated[i], 0), if (tab$note[i] == "failed") "failed" else tab$edges[i], fmt(tab$shared[i], 0), fmt(tab$jaccard[i]), fmt(tab$sign_agree[i])))
  L <- c(L, "\\hline", "\\end{tabular}", "\\end{adjustbox}", "\\end{table}")
  writeLines(display_text(L), file.path(out_dir, paste0(prefix, "_methods_comparison.tex")))

  ## every ordered pair, every method (csv) and the edges of MR.RGM / MR.RGM+ with the support of the other methods (tex)
  idx <- which(off, arr.ind = TRUE)
  long <- data.frame(from = labels[idx[, 2]], to = labels[idx[, 1]])
  for (m in names(all_m)) { x <- all_m[[m]]
    long[[paste0(m, "_edge")]]  <- if (inherits(x, "failed")) NA else x$edge[idx]
    long[[paste0(m, "_est")]]   <- if (inherits(x, "failed")) NA else x$est[idx]
    long[[paste0(m, "_score")]] <- if (inherits(x, "failed")) NA else x$score[idx] }
  write.csv(long, file.path(out_dir, paste0(prefix, "_methods_edges.csv")), row.names = FALSE)
  rep_e <- long[long$`MR.RGM_edge`, , drop = FALSE]                  # the edges of the main analysis (every pair: the csv)
  rep_e <- rep_e[order(-rep_e$`MR.RGM_score`), , drop = FALSE]
  others <- setdiff(names(all_m), c("MR.RGM", "MR.RGM+"))
  ab <- c(MR.RGM_NoConf = "NoConf", MVMR = "MVMR", MR2G = "MR2G", MrDAG = "MrDAG", GrIVET = "GrIVET", `MR-IVW` = "IVW", `MR-SimpleMedian` = "SMed",
          `MR-WeightedMedian` = "WMed", OneSampleMR = "TSRI", mrbayes = "mrbayes")
  hdr <- c("\\textbf{Edge}", "\\textbf{PIP \\texttt{MR.RGM}}", if (!is.null(fit_plus)) "\\textbf{PIP \\texttt{MR.RGM+}}",
           sprintf("\\textbf{%s}", ifelse(others %in% names(ab), ab[others], others)))
  cellm <- function(r, m) { e <- rep_e[r, paste0(m, "_edge")]; v <- rep_e[r, paste0(m, "_est")]
    if (is.na(e) || is.na(rep_e[r, paste0(m, "_score")])) "--" else if (isTRUE(e)) (if (isTRUE(v < 0)) "$-$" else "$+$") else "" }
  L <- c("\\begin{table}[htb]", "\\centering", "\\small",
         sprintf("\\caption{%s: causal edges reported by \\texttt{MR.RGM} (PIP $> %.2f$), their PIP under \\texttt{MR.RGM+}, and whether each other method reports them under its rule ($+$/$-$: reported, with the sign of its estimate; blank: not reported; --: the method failed or could not estimate the pair). Abbreviations: NoConf = \\texttt{MR.RGM\\_NoConf}, IVW / SMed / WMed = \\texttt{MendelianRandomization} IVW, simple and weighted median, TSRI = two-stage residual inclusion (\\texttt{OneSampleMR}).%s}", data_name, thr, note),
         sprintf("\\label{tab:%s_edge_support}", prefix), "\\begin{adjustbox}{max width=\\textwidth}",
         sprintf("\\begin{tabular}{|l|%s}", paste(rep("c|", length(hdr) - 1), collapse = "")), "\\hline", paste(paste(hdr, collapse = " & "), "\\\\"), "\\hline")
  for (r in seq_len(nrow(rep_e)))
    L <- c(L, paste(paste(c(sprintf("\\textit{%s} $\\rightarrow$ \\textit{%s}", tex_esc(rep_e$from[r]), tex_esc(rep_e$to[r])),
                            sprintf("%.2f", rep_e$`MR.RGM_score`[r]), if (!is.null(fit_plus)) sprintf("%.2f", rep_e$`MR.RGM+_score`[r]),
                            vapply(others, function(m) cellm(r, m), "")), collapse = " & "), "\\\\"))
  L <- c(L, "\\hline", "\\end{tabular}", "\\end{adjustbox}", "\\end{table}")
  writeLines(display_text(L), file.path(out_dir, paste0(prefix, "_edge_support.tex")))
  invisible(tab)
}

## ---------------------------------------------------------------------------
## Several chains of one model: run, check convergence, pool
##   run_one(seed) returns an RGM() fit (with $warnings). The chains differ in their random seeds; every chain starts
##   from the package's starting point (A = 0, least-squares instrument effects, residual error covariance).
##   Convergence (Vehtari, Gelman, Simpson, Carpenter and Buerkner 2021, package posterior): rank-normalized split-R-hat
##   and bulk effective sample size (ESS) of the log-likelihood and of every causal edge indicator; the PIPs of every pair
##   of chains. The chains are converged if every split-R-hat is <= RHAT_MAX (1.01).
##   Pooled fit: posterior means averaged over the chains (equal numbers of draws, so this is the posterior mean of the
##   pooled draws), GammaPst and LLPst concatenated, warnings of every chain. Pooled results are interpreted only for
##   converged chains (OneK1K_Results.R).
## ---------------------------------------------------------------------------
RHAT_MAX <- 1.01
fit_chains <- function(run_one, seeds) {
  if (!requireNamespace("posterior", quietly = TRUE)) stop("Package 'posterior' is required: install.packages(\"posterior\")")
  fits <- lapply(seeds, run_one)
  K <- length(fits); p <- nrow(fits[[1]]$GammaEst); m <- dim(fits[[1]]$GammaPst)[3]
  off <- which(row(fits[[1]]$GammaEst) != col(fits[[1]]$GammaEst))
  ## convergence
  LL <- sapply(fits, function(f) f$LLPst)                                   # draws x chains
  edge_draws <- function(i) sapply(fits, function(f) matrix(f$GammaPst, p * p, m)[i, ])
  rh_e <- vapply(off, function(i) { x <- edge_draws(i); if (var(as.vector(x)) == 0) NA_real_ else posterior::rhat(x) }, numeric(1))
  es_e <- vapply(off, function(i) { x <- edge_draws(i); if (var(as.vector(x)) == 0) NA_real_ else posterior::ess_bulk(x) }, numeric(1))
  pip <- lapply(fits, function(f) f$GammaEst[off])
  pr  <- utils::combn(K, 2)
  diffs <- apply(pr, 2, function(j) max(abs(pip[[j[1]]] - pip[[j[2]]])))
  cors  <- apply(pr, 2, function(j) suppressWarnings(cor(pip[[j[1]]], pip[[j[2]]])))
  rhat_ll <- posterior::rhat(LL); rhat_max <- suppressWarnings(max(rh_e, na.rm = TRUE))
  RH <- matrix(NA_real_, p, p); RH[off] <- rh_e
  conv <- list(seeds = seeds, chains = K, draws_per_chain = m,
               rhat_loglik = rhat_ll, rhat_edges_max = rhat_max, rhat_edges = RH,
               ess_bulk_loglik = posterior::ess_bulk(LL), ess_bulk_edges_min = suppressWarnings(min(es_e, na.rm = TRUE)),
               pip_pairwise_max_diff = max(diffs), pip_pairwise_min_cor = min(cors),
               pip_chains = lapply(fits, `[[`, "GammaEst"), effect_chains = lapply(fits, `[[`, "AEst"), loglik_chains = LL,
               accept_A_chains = vapply(fits, function(f) if (is.null(f$AccptA)) NA_real_ else f$AccptA, numeric(1)),
               converged = is.finite(rhat_ll) && rhat_ll <= RHAT_MAX && (!is.finite(rhat_max) || rhat_max <= RHAT_MAX),
               criterion = sprintf("rank-normalized split-R-hat <= %.2f for the log-likelihood and every causal edge indicator (Vehtari et al. 2021)", RHAT_MAX))
  ## pooled fit
  pooled <- fits[[1]]
  for (nm in names(pooled)) {
    if (nm %in% c("GammaPst", "LLPst", "settings", "diagnostics", "warnings")) next
    v <- lapply(fits, `[[`, nm)
    if (all(vapply(v, function(x) is.numeric(x) && identical(dim(x), dim(v[[1]])) && length(x) == length(v[[1]]), logical(1))))
      pooled[[nm]] <- Reduce(`+`, v) / K
  }
  pooled$GammaPst <- array(unlist(lapply(fits, `[[`, "GammaPst")), dim = c(p, p, m * K), dimnames = dimnames(fits[[1]]$GammaPst))
  pooled$LLPst    <- as.vector(LL)
  pooled$warnings <- unique(unlist(lapply(fits, `[[`, "warnings")))
  message(sprintf("%d chains: split-R-hat log-likelihood %.3f, largest over edges %.3f; bulk ESS log-likelihood %.0f, smallest over edges %.0f; between chains: largest PIP difference %.3f, smallest PIP correlation %.3f -> %s",
                  K, conv$rhat_loglik, conv$rhat_edges_max, conv$ess_bulk_loglik, conv$ess_bulk_edges_min,
                  conv$pip_pairwise_max_diff, conv$pip_pairwise_min_cor, if (conv$converged) "converged (pooled)" else "NOT converged"))
  list(fit = pooled, convergence = conv)
}
