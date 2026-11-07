#' Posterior probability of a network motif
#'
#' @description \code{NetworkMotif} quantifies the uncertainty of a given network structure among the response
#'              variables: it returns the proportion of posterior draws in which every edge of the motif is present,
#'              that is, the posterior probability that all edges of the motif are present simultaneously. The input
#'              draws are the \code{GammaPst} output of \code{\link{RGM}}.
#'
#' @param Gamma A p * p binary matrix (numeric 0/1 or logical) with a zero diagonal, specifying the motif: \code{Gamma[i, j] = 1} for the directed edge from response j
#'              to response i (the same orientation as \code{A} and \code{GammaEst} in \code{\link{RGM}}). Entries equal to
#'              0 are ignored, so the motif can be any subgraph (a single edge, a cascade, a feed-forward loop, a
#'              feedback loop, ...); only the presence of the specified edges is required, not the absence of others.
#' @param GammaPst A p * p * n_pst array of posterior draws of the causal edge indicators, the \code{GammaPst} component
#'              returned by \code{\link{RGM}}.
#'
#' @return A number between 0 and 1: the proportion of posterior draws that contain all edges of \code{Gamma}.
#'
#'
#'
#' @export
#'
#' @examples
#' # Simulate a cascade Y1 -> Y2 -> Y3 with two own instruments per response
#' set.seed(9154)
#' n = 2000; p = 3; k = 6
#' A = matrix(0, p, p); A[2, 1] = 0.3; A[3, 2] = -0.3     # A[i, j]: effect of Y_j on Y_i
#' D = matrix(0, p, k)
#' for (i in 1:p) D[i, (2 * i - 1):(2 * i)] = 1
#' X = matrix(rnorm(n * k), n, k)
#' Y = (X %*% t(D) + matrix(rnorm(n * p), n, p)) %*% t(solve(diag(p) - A))
#' X = scale(X); Y = scale(Y)                              # the model has no intercept
#'
#' # Fit MR.RGM (a short chain to keep the example fast; the defaults are the recommended length)
#' fit = RGM(X = X, Y = Y, D = D, nIter = 2000, nBurnin = 1000, Thin = 2)
#'
#' # Posterior probability of the cascade Y1 -> Y2 -> Y3 (both edges present)
#' Gamma = matrix(0, p, p); Gamma[2, 1] = 1; Gamma[3, 2] = 1
#' NetworkMotif(Gamma = Gamma, GammaPst = fit$GammaPst)
#'
#' # Posterior probability of the feedback loop Y1 -> Y2 -> Y3 -> Y1
#' Gamma[1, 3] = 1
#' NetworkMotif(Gamma = Gamma, GammaPst = fit$GammaPst)
#'
#' @references
#' Ni, Y., Ji, Y., & M\enc{ü}{ue}ller, P. (2018).
#' Reciprocal graphical models for integrative gene regulatory network analysis.
#' \emph{Bayesian Analysis},
#' \strong{13(4)}, 1095-1110.
#' \doi{10.1214/17-BA1087}.
NetworkMotif = function(Gamma, GammaPst) {

  # Check whether GammaPst is a numeric array with three dimensions
  if (!is.numeric(GammaPst) || !is.array(GammaPst) || length(dim(GammaPst)) != 3) {

    # Stop with an error
    stop("GammaPst must be a numeric array with three dimensions.")

  }

  # Calculate number of rows of GammaPst
  p = nrow(GammaPst)

  # Check whether number of rows of GammaPst is same as number of columns of GammaPst
  if(ncol(GammaPst) != p){

    # Stop with an error
    stop("Number of rows and columns of GammaPst should be equal.")

  }

  # A logical motif matrix is accepted
  if (is.logical(Gamma) && is.matrix(Gamma)) Gamma = Gamma * 1

  # Check whether Gamma is a numeric matrix
  if(!is.numeric(Gamma) || !is.matrix(Gamma)){

    # Stop with an error
    stop("Gamma should be a numeric matrix.")

  }

  # Check whether Gamma is a square matrix with dimension p * p
  if((nrow(Gamma) != p) || (ncol(Gamma) != p)) {

    # Stop with an error
    stop("Gamma should be a square matrix with number of rows and columns equal to number of rows of GammaPst.")

  }

  # Entries of Gamma other than 0 and 1 are not meaningful
  if (!all(Gamma %in% c(0, 1))) {

    # Stop with an error
    stop("Gamma should contain only 0 and 1 (1 marks the edges of the motif).")

  }

  # The motif cannot contain self-loops
  if (any(diag(Gamma) != 0)) stop("Gamma must have a zero diagonal (no self-loops).")

  # The draws must be binary, with at least one draw
  if (dim(GammaPst)[3] < 1) stop("GammaPst contains no posterior draws.")
  if (!all(GammaPst %in% c(0, 1))) stop("GammaPst should contain only 0 and 1 (the GammaPst component returned by RGM()).")

  # A motif without edges is present in every draw
  if (!any(Gamma == 1)) return(1)

  # Return Network Motif
  return(NetworkMotif_cpp(Gamma = Gamma, Gamma_Pst = GammaPst))

}
