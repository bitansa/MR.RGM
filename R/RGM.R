#' Fitting Multivariate Bidirectional Mendelian Randomization Networks
#' Using Bayesian Directed Cyclic Graphical Models
#'
#' @description \code{RGM} fits the Bayesian directed cyclic graphical model
#'   \deqn{Y = A Y + B X + C U + E,\qquad E \sim N_p(0, \Sigma^{\ast}),}
#'   where \code{Y} are \code{p} response variables (traits), \code{X} are \code{k} instrumental variables
#'   (for example genetic variants), \code{U} are \code{l} optional measured covariates, \code{A} (\code{p x p},
#'   zero diagonal) contains the direct causal effects among the responses (\code{A[i, j]} is the effect of
#'   \code{Y_j} on \code{Y_i}; directed cycles are allowed), \code{B} (\code{p x k}) the instrument effects,
#'   \code{C} (\code{p x l}) the covariate effects and \eqn{\Sigma^{\ast}} the error covariance, whose off-diagonal
#'   entries represent unmeasured confounding between responses.
#'
#'   Two versions of the method differ only in how the instruments are mapped to the responses (in the accompanying
#'   methods paper they are called IV-RGM and IV-RGM+, since the model applies to instrumental variables in general):
#'   \describe{
#'     \item{\strong{MR.RGM} (known instrument map)}{\code{D} is a binary \code{p x k} matrix with
#'       \code{D[i, j] = 1} if instrument \code{j} is allowed to affect response \code{i} (for example a
#'       cis-eQTL for its own gene). Entries of \code{B} with \code{D = 0} are fixed at zero (the exclusion
#'       restriction) and the allowed entries receive a slab prior \code{N(0, eta)} without any selection:
#'       every allowed instrument stays in the model. This is the default when \code{D} contains zeros
#'       (\code{selectB = FALSE}).}
#'     \item{\strong{MR.RGM+} (instrument selection)}{\code{D} is a matrix of ones: every instrument may affect
#'       every response, and each entry of \code{B} receives a spike-and-slab prior whose inclusion indicator is
#'       sampled, so that the instrument--response map is learned from the data (horizontal pleiotropy). This is
#'       the default when \code{D} has no zeros (\code{selectB = TRUE}).}
#'   }
#'   The argument \code{selectB} overrides these defaults (see its description).
#'
#'   Inference is by MCMC on the sufficient statistics of the Gaussian model, so the function accepts
#'   individual-level data (\code{X}, \code{Y}, optionally \code{U}), the second-moment matrices
#'   (\code{Syy}, \code{Syx}, \code{Sxx}, optionally \code{Syu}, \code{Sxu}, \code{Suu}) computed from one sample,
#'   or the matrices \code{Sxx}, \code{Beta} and \code{SigmaHat} of marginal instrument--response regressions;
#'   responses alone (\code{Y} or \code{Syy}) fit the model without instruments. The marginal-regression input does
#'   not contain the covariances between the responses: they are reconstructed under independent errors
#'   (no unmeasured confounding), so this input requires \code{SigmaStarModel = "diagonal"}, no covariates and, for
#'   every response, an instrument that affects only that response, and it is exact only for uncorrelated
#'   instruments. Use individual-level data or the second-moment matrices to model confounding.
#'   The model has no intercept and the stochastic-search prior on \eqn{\Sigma^{\ast}} assumes responses on roughly
#'   unit scale: center and scale the columns of \code{Y}, \code{X} and \code{U} (or compute the summary
#'   matrices from centered and scaled data) before fitting.
#'
#'   The output contains posterior means of all parameters, posterior inclusion probabilities of the causal
#'   edges (\code{GammaEst}), of the instrument effects (\code{PhiEst}, MR.RGM+ only) and of the confounding
#'   links (\code{ZEst}, SSSL only), the graph estimated by thresholding the inclusion probabilities at 0.5,
#'   the posterior draws of the causal edge indicators (\code{GammaPst}, used by \code{\link{NetworkMotif}}),
#'   the log-likelihood trace, acceptance rates and convergence/identification diagnostics. Warnings are issued for
#'   signs of weak identification or poor convergence (see \code{diagnostics} in the Value section and the Note).
#'
#'   \strong{MCMC algorithm.} Every iteration updates, in turn, the instrument effects \code{B}, the covariate
#'   effects \code{C}, the error covariance \eqn{\Sigma^{\ast}} (with its confounding indicators) and the causal
#'   effects \code{A}, each with its prior hyperparameters. Given the other parameters, the likelihood is quadratic in
#'   \code{B} and its (spike-and-slab or slab) prior is normal given the inclusion indicators, so each row of
#'   \code{B} is drawn exactly from its multivariate normal full conditional (Gibbs; \code{gibbsB = TRUE}). The
#'   entries of \code{A} enter the likelihood also through \eqn{\log|\det(I - A)|} and are updated one at a time by
#'   random-walk Metropolis--Hastings with adaptive proposals: each entry has its own proposal variance, started at
#'   \eqn{2.38^2} times an approximation of its posterior variance and tuned during the burn-in towards an
#'   acceptance rate of 0.44 (adaptive Metropolis-within-Gibbs of Roberts and Rosenthal, 2009, with larger adaptation steps early in the burn-in); the variances are frozen after the burn-in, so the kept
#'   draws come from a fixed Markov kernel. With more than 100 allowed instruments for a response, its row of \code{B}
#'   is drawn entry by entry (each entry exactly from its scalar full conditional) instead of jointly. The same adaptive scheme updates \code{B} when it cannot be drawn by Gibbs
#'   (Threshold prior with selection, or \code{gibbsB = FALSE}). The chain starts at \code{A = 0}, at the
#'   least-squares effects of each response on its allowed instruments (known map; \code{B = 0} with selection) and
#'   at the corresponding residual covariance.
#'
#'   \strong{Chain length and convergence.} One call runs one chain. The defaults (50,000 iterations, 10,000 burn-in,
#'   thinning 10) are the settings validated in the simulation study of the accompanying paper for networks of up to
#'   about 20 responses; larger networks or weak instruments may need longer chains. For a reported analysis, run
#'   several chains (for example four) with different seeds
#'   (\code{set.seed()} before each call), check their agreement (for example the rank-normalized split-\eqn{\hat R}
#'   of \code{LLPst} and of the edge indicators in \code{GammaPst}, computed with the package \pkg{posterior}, and the
#'   posterior inclusion probabilities of the chains), and pool the chains only if they agree. The function also
#'   reports diagnostics and warnings of its own (see \code{diagnostics} in the value).
#'
#' @param X A numeric matrix of dimension n * k of instrumental variables (rows = observations). Default NULL.
#' @param Y A numeric matrix of dimension n * p of response variables (rows = observations). Default NULL.
#' @param U A numeric matrix of dimension n * l of measured covariates (rows = observations). Default NULL.
#' @param Syy A p * p matrix, `t(Y) %*% Y / n`, computed from centered (and preferably scaled) responses. Default NULL.
#' @param Syx A p * k matrix, `t(Y) %*% X / n`. Default NULL.
#' @param Sxx A k * k matrix, `t(X) %*% X / n`. Default NULL.
#' @param Syu A p * l matrix, `t(Y) %*% U / n`. Default NULL.
#' @param Sxu A k * l matrix, `t(X) %*% U / n`. Required together with \code{Syu} and \code{Suu} when instruments and covariates are both given as summary matrices. Default NULL.
#' @param Suu An l * l matrix, `t(U) %*% U / n`. Default NULL.
#' @param Beta A p * k matrix of regression coefficients of each response on each instrument separately (marginal regressions on centered data). Used with \code{Sxx} and \code{SigmaHat}; this input treats the marginal effects as joint effects, which is exact for uncorrelated instruments, and reconstructs the covariances between the responses under independent errors. It requires \code{SigmaStarModel = "diagonal"}, no covariates, and a map \code{D} in which every response has an instrument that affects only that response (so it is not available for MR.RGM+ with \code{D} all ones). Default NULL.
#' @param SigmaHat A p * k matrix of the residual variances of the same marginal regressions: the residual sum of squares divided by n (not by the residual degrees of freedom). Default NULL.
#' @param D A binary p * k matrix (numeric 0/1 or logical): \code{D[i, j] = 1} if instrument j is allowed to affect response i. A matrix of ones gives MR.RGM+; a matrix with zeros gives MR.RGM (see Description). Required whenever instruments are given; not used without instruments.
#' @param n A positive integer, the number of observations. Required when summary-level matrices are the input; with individual-level data it is \code{nrow(Y)}.
#' @param nIter A positive integer, the number of MCMC iterations. Default 50,000 (the setting of the simulation study of the accompanying paper).
#' @param nBurnin A non-negative integer smaller than \code{nIter}, the number of burn-in iterations. Default 10,000. The proposal variances of the Metropolis updates are tuned during the burn-in (in batches of 50 iterations); a burn-in of a few thousand iterations or more is recommended. With the spike-and-slab prior and a known instrument map (\code{D} with zeros, \code{selectB = FALSE}), the inclusion indicators of \code{A} are held at 1 during the first half of the burn-in, so that every causal effect first moves under its slab prior before selection starts.
#' @param Thin A positive integer, the thinning interval of the kept draws (draws are kept at iterations \code{nBurnin + Thin, nBurnin + 2 Thin, ...}). It must not exceed \code{nIter - nBurnin}. Default 10 (4,000 kept draws with the default chain length).
#' @param prior The prior on the graph structure: \code{"Spike and Slab"} (default) or \code{"Threshold"}.
#' @param SigmaStarModel The model for the error covariance \eqn{\Sigma^{\ast}}: \code{"SSSL"} (default; stochastic search structure learning prior with correlated errors and inclusion indicators for the covariances), \code{"IW"} (inverse-Wishart prior \eqn{IW(p + 1, I_p)}, correlated errors) or \code{"diagonal"} (independent errors with inverse-gamma priors \eqn{IG(aSigma, bSigma)} on the variances).
#' @param aRho,bRho Positive shape parameters of the Beta prior on the inclusion probability of a causal edge (spike-and-slab prior on \code{A}). Defaults 3 and 1 (prior inclusion probability 0.75).
#' @param nu1 A positive scalar, the spike variance of the spike-and-slab prior on \code{A} relative to the slab variance (spike-and-slab prior only). Default 0.001.
#' @param aPsi,bPsi Positive shape parameters of the Beta prior on the inclusion probability of an instrument effect (spike-and-slab prior with selection of the instrument effects, e.g. MR.RGM+, only). Defaults 0.5 and 0.5.
#' @param nu2 A positive scalar, the spike variance of the spike-and-slab prior on \code{B} relative to the slab variance (spike-and-slab prior with selection, e.g. MR.RGM+, only). Default 0.0001.
#' @param aSigma,bSigma Positive shape and rate parameters of the inverse-gamma prior on the error variances when \code{SigmaStarModel = "diagonal"}. Defaults 0.01 and 0.01.
#' @param PropVarA,PropVarB Starting variances of the adaptive random-walk Metropolis--Hastings proposals for the entries of \code{A} and (when \code{B} is not drawn by Gibbs) \code{B}; the proposals are tuned during the burn-in. \code{NULL} (default): an automatic start, \eqn{2.38^2 / (n \Omega_{ii} S_{jj} + 1)} for entry \code{(i, j)}, where \eqn{\Omega} is the inverse of the starting error covariance and \code{S} is \code{Syy} (for \code{A}) or \code{Sxx} (for \code{B}); it scales with the sample size. A positive number is used as the starting variance of every entry.
#' @param gibbsB Logical. \code{TRUE} (default): the instrument effects \code{B} are drawn exactly from their normal full conditional (Gibbs), row by row; with the Threshold prior this applies when the instrument effects are not selected (\code{selectB = FALSE}). \code{FALSE}: adaptive random-walk Metropolis--Hastings for every entry of \code{B}.
#' @param selectB Logical or NULL. Whether the instrument effects in \code{B} are selected by a spike-and-slab (or threshold) prior. \code{NULL} (default) means \code{all(D == 1)}: selection for MR.RGM+ (\code{D} all ones) and no selection for MR.RGM (\code{D} with zeros). \code{TRUE} with a \code{D} that contains zeros selects among the allowed entries only; \code{FALSE} with \code{D} all ones gives every entry a slab prior without selection. The value used is returned in \code{settings}.
#' @param tauC A positive scalar, the column variance of the matrix-normal prior \code{C ~ MN(0, SigmaStar, tauC I_l)} on the covariate effects. Default 2.
#' @param omega1,omega2 Positive slab and spike variances of the SSSL prior on the off-diagonal entries of \eqn{\Sigma^{\ast}} (\code{SigmaStarModel = "SSSL"}). \code{omega1} must exceed \code{omega2}. Defaults 5 and 0.005; they assume responses on unit scale.
#' @param varpi The prior inclusion probability of a confounding link under the SSSL prior, in (0, 1). Default 0.5.
#' @param lambda A positive scalar; under the SSSL prior the diagonal entries of \eqn{\Sigma^{\ast}} have an exponential prior with rate \code{lambda / 2} (the parametrization of Wang, 2015). Default 1.
#'
#' @return A list with the following components (components that do not apply to the chosen prior or model are omitted):
#'
#' \item{AEst}{p * p matrix, posterior mean of \code{A} (\code{AEst[i, j]} is the estimated direct effect of \code{Y_j} on \code{Y_i}).}
#' \item{BEst}{p * k matrix, posterior mean of \code{B}.}
#' \item{CEst}{p * l matrix, posterior mean of \code{C} (with covariates).}
#' \item{zAEst}{p * p binary matrix, \code{GammaEst > 0.5} (median probability model).}
#' \item{zBEst}{p * k binary matrix, \code{PhiEst > 0.5}; equal to \code{D} without selection.}
#' \item{Graph}{An \code{igraph} object of the estimated causal network among the responses (edges of \code{zAEst}, labelled by \code{AEst}).}
#' \item{A0Est, B0Est}{Posterior means of the unthresholded \code{A} and \code{B} (Threshold prior).}
#' \item{GammaEst}{p * p matrix of posterior inclusion probabilities of the causal edges (\code{GammaEst[i, j]} for the edge \code{Y_j -> Y_i}).}
#' \item{TauEst}{p * p matrix, posterior mean of the slab variances of \code{A}.}
#' \item{RhoEst}{p * p matrix, posterior mean of the edge inclusion probabilities \code{rho} (spike-and-slab prior).}
#' \item{PhiEst}{p * k matrix of posterior inclusion probabilities of the instrument effects when they are selected (MR.RGM+); equal to \code{D} otherwise.}
#' \item{EtaEst}{p * k matrix, posterior mean of the slab variances of \code{B}.}
#' \item{PsiEst}{p * k matrix, posterior mean of the inclusion probabilities \code{psi} of the instrument effects (spike-and-slab prior with selection only).}
#' \item{tAEst, tBEst}{Posterior means of the thresholds (Threshold prior; \code{tBEst} with selection only).}
#' \item{SigmaEst}{The estimated error covariance: a p * p matrix for \code{"SSSL"} and \code{"IW"}, a vector of length p of error variances for \code{"diagonal"}.}
#' \item{ZEst}{p * p matrix of posterior inclusion probabilities of the confounding links (\code{"SSSL"} only; the diagonal is 1).}
#' \item{AccptA, AccptB}{Acceptance rates (in percent, over the iterations after burn-in) of the Metropolis--Hastings updates of \code{A} and \code{B}; \code{AccptB} is omitted when \code{B} is drawn by Gibbs.}
#' \item{AccpttA, AccpttB}{Acceptance rates (in percent, after burn-in) of the threshold updates (Threshold prior; \code{AccpttB} only with selection of the instrument effects).}
#' \item{PropVarA, PropVarB}{The tuned (frozen) proposal variances of the entries of \code{A} (p * p, zero diagonal) and, when \code{B} is updated by Metropolis--Hastings, of \code{B} (p * k, zero where \code{D = 0}).}
#' \item{LLPst}{Numeric vector of the log-likelihood at the kept draws.}
#' \item{GammaPst}{p * p * nPst array of the kept draws of the causal edge indicators (input of \code{\link{NetworkMotif}}).}
#' \item{diagnostics}{A list of convergence and identification diagnostics, each with a warning when it signals a problem: \code{instrumentF} (only without instrument selection, \code{selectB = FALSE}: first-stage F statistic of each response on its allowed instruments, computed from the second-moment matrices without adjusting for the covariates, \code{NA} where it cannot be computed; warning below 10), \code{spectralRadiusA} (spectral radius of \code{AEst}; warning at 1 or above, where the estimated feedback loops do not converge and the equilibrium interpretation fails), \code{halvesMaxDiff} (largest difference of the causal inclusion probabilities between the two halves of the kept draws; with at least 200 kept draws, warning when, for at least one edge, the difference exceeds 0.2 and four times that edge's batch-means Monte Carlo standard error; returned only with at least 200 kept draws; a heuristic, which can also fire for edges that mix slowly), \code{halvesMCSE} (the Monte Carlo standard error of that largest difference), \code{acceptanceA} (= \code{AccptA}; warning below 10 percent) and \code{errorVarianceRatio} (ratio of each estimated error variance to the second moment of its response; warning above 5).}
#' \item{settings}{A list of the settings used: \code{input} (\code{"Y"}, \code{"XY"}, \code{"S"} or \code{"Beta"}), \code{prior}, \code{SigmaStarModel}, \code{instruments}, \code{covariates}, \code{selectB}, \code{gibbsB} (whether \code{B} was drawn by Gibbs), \code{n}, \code{p}, \code{k} (number of instruments), \code{l} (number of covariates), \code{nIter}, \code{nBurnin}, \code{Thin}, the hyperparameters, the starting proposal variances (\code{"auto"} for the automatic start) and the package version. The Threshold prior also uses a fixed upper bound of 0.1 on the thresholds and a truncated-normal proposal with standard deviation 0.1 for them.}
#'
#' @export
#'
#' @examples
#' # ------------------------------------------------------------------------------------------
#' # Simulated data: Y1 -> Y2 -> Y3 and Y3 -> Y1 (a feedback loop), two own instruments per
#' # response, one covariate, and an unmeasured confounder shared by Y1 and Y3
#' set.seed(9154)
#' n = 2000; p = 3; k = 6; l = 1
#' A = matrix(0, p, p); A[2, 1] = 0.3; A[3, 2] = -0.3; A[1, 3] = 0.2   # A[i, j]: effect of Y_j on Y_i
#' D = matrix(0, p, k)                                                  # known instrument map
#' for (i in 1:p) D[i, (2 * i - 1):(2 * i)] = 1
#' X = matrix(rnorm(n * k), n, k); U = matrix(rnorm(n * l), n, l); W = rnorm(n)
#' E = cbind(W, 0, W) + matrix(rnorm(n * p), n, p)                      # confounding of Y1 and Y3
#' Y = (X %*% t(D) + U %*% matrix(0.5, l, p) + E) %*% t(solve(diag(p) - A))
#' colnames(Y) = paste0("Y", 1:p); colnames(X) = paste0("SNP", 1:k); colnames(U) = "age"
#' X = scale(X); Y = scale(Y); U = scale(U)                            # the model has no intercept
#'
#' # Example 1: MR.RGM (known map D) from individual-level data, spike-and-slab prior and
#' # SSSL error covariance (the defaults); a short chain to keep the example fast, whereas the
#' # default chain length is the recommended one
#' fit = RGM(X = X, Y = Y, U = U, D = D, nIter = 2000, nBurnin = 1000, Thin = 2)
#' round(fit$GammaEst, 2)   # posterior inclusion probabilities of the causal edges
#' round(fit$AEst, 2)       # posterior means of the causal effects
#' round(fit$ZEst, 2)       # posterior inclusion probabilities of the confounding links
#' fit$diagnostics          # weak-instrument, stability and convergence diagnostics
#' plot(fit$Graph, main = "Estimated causal network")
#' plot(fit$LLPst, type = "l", xlab = "Kept draw", ylab = "Log-likelihood")
#'
#' # Example 2: the same model from the second-moment matrices of the same sample
#' fit_s = RGM(Syy = crossprod(Y) / n, Syx = crossprod(Y, X) / n, Sxx = crossprod(X) / n,
#'             Syu = crossprod(Y, U) / n, Sxu = crossprod(X, U) / n, Suu = crossprod(U) / n,
#'             D = D, n = n, nIter = 2000, nBurnin = 1000, Thin = 2)
#' round(fit_s$GammaEst, 2)
#'
#' # Example 3: MR.RGM+ (D = 1): every instrument may affect every response and the instrument
#' # effects are selected by a spike-and-slab prior (horizontal pleiotropy)
#' fit_plus = RGM(X = X, Y = Y, U = U, D = matrix(1, p, k), nIter = 2000, nBurnin = 1000, Thin = 2)
#' round(fit_plus$PhiEst, 2)   # posterior inclusion probabilities of the instrument effects
#' fit_plus$settings$selectB   # TRUE
#'
#' # Example 4: threshold prior on A and a diagonal error covariance
#' fit_t = RGM(X = X, Y = Y, D = D, prior = "Threshold", SigmaStarModel = "diagonal",
#'             nIter = 2000, nBurnin = 1000, Thin = 2)
#' fit_t$zAEst
#'
#' # Example 5: marginal-regression input (Sxx, Beta, SigmaHat), with uncorrelated instruments,
#' # no covariates and independent errors (SigmaStarModel = "diagonal" is required)
#' Y0 = scale((X %*% t(D) + matrix(rnorm(n * p), n, p)) %*% t(solve(diag(p) - A)))
#' Beta = SigmaHat = matrix(0, p, k)
#' for (i in 1:p) for (j in 1:k) {
#'   b = sum(X[, j] * Y0[, i]) / sum(X[, j]^2)                          # centered data, no intercept
#'   Beta[i, j] = b; SigmaHat[i, j] = sum((Y0[, i] - b * X[, j])^2) / n  # residual sum of squares / n
#' }
#' fit_b = RGM(Sxx = crossprod(X) / n, Beta = Beta, SigmaHat = SigmaHat, D = D, n = n,
#'             SigmaStarModel = "diagonal", nIter = 2000, nBurnin = 1000, Thin = 2)
#' round(fit_b$GammaEst, 2)
#'
#' @note \strong{Weak identification.} With \code{D = 1} (MR.RGM+), the instrument effects are identified only
#'   through sparsity and an instrument-diversity condition, which require instruments that explain a
#'   non-negligible share of the variance of several responses. With weak, response-specific instruments
#'   (for example cis-eQTLs) the posterior can stay close to the prior; in that case \code{RGM} warns (error
#'   variances far above the second moments of the responses) and the instrument map \code{D} should encode
#'   which instrument affects which response (MR.RGM). Compare the posterior inclusion probabilities with their
#'   prior values (0.75 for causal edges under the default \code{aRho = 3, bRho = 1}; 0.5 for confounding
#'   links) and the causal inclusion probabilities computed from the two halves of \code{GammaPst}.
#'
#'   \strong{Identification of MR.RGM.} With a known map, the causal effects are identified when, for every
#'   response, the instruments excluded from its equation satisfy the classical rank condition; a sufficient
#'   condition is that every response has at least one instrument that affects only that response. \code{RGM}
#'   warns when this sufficient condition fails and no selection is performed.
#'
#' @references
#' Roberts, G. O., & Rosenthal, J. S. (2009).
#' Examples of adaptive MCMC.
#' \emph{Journal of Computational and Graphical Statistics},
#' \strong{18(2)}, 349-367.
#' \doi{10.1198/jcgs.2009.06134}.
#'
#' Wang, H. (2015).
#' Scaling it up: stochastic search structure learning in graphical models.
#' \emph{Bayesian Analysis},
#' \strong{10(2)}, 351-377.
#' \doi{10.1214/14-BA916}.
#'
#' Ni, Y., Ji, Y., & M\enc{ü}{ue}ller, P. (2018).
#' Reciprocal graphical models for integrative gene regulatory network analysis.
#' \emph{Bayesian Analysis},
#' \strong{13(4)}, 1095-1110.
#' \doi{10.1214/17-BA1087}.
#'
#' Sarkar, B. & Ni, Y. (2025).
#' MR.RGM: an R package for fitting Bayesian multivariate bidirectional Mendelian randomization networks.
#' \emph{Bioinformatics}, \strong{41(4)}, btaf130.
#' \doi{10.1093/bioinformatics/btaf130}.
RGM = function(X = NULL, Y = NULL, U = NULL, Syy = NULL, Syx = NULL, Sxx = NULL, Syu = NULL, Sxu = NULL, Suu = NULL,
               Beta = NULL, SigmaHat = NULL, D, n, nIter = 50000, nBurnin = 10000, Thin = 10,
               prior = c("Spike and Slab", "Threshold"), SigmaStarModel = c("SSSL", "IW", "diagonal"),
               aRho = 3, bRho = 1, nu1 = 0.001, aPsi = 0.5, bPsi = 0.5, nu2 = 0.0001, aSigma = 0.01, bSigma = 0.01,
               PropVarA = NULL, PropVarB = NULL, gibbsB = TRUE,
               selectB = NULL, tauC = 2, omega1 = 5, omega2 = 0.005, varpi = 0.5, lambda = 1) {

  # Resolve the prior and the error-covariance model (the first option is the default)
  prior = match.arg(prior)
  SigmaStarModel = match.arg(SigmaStarModel)

  # ---------------------------------------------------------------------------------------------------------------
  # 1. Which input was given?  "Y" (responses only), "XY" (individual-level), "S" (second moments) or "Beta"
  # ---------------------------------------------------------------------------------------------------------------
  if (!is.null(X) && is.null(Y)) stop("X was given without Y: give (X, Y), or the summary matrices (Syy, Syx, Sxx).")
  supplied = function(...) vapply(list(...), Negate(is.null), logical(1))
  if (!is.null(Y) && any(supplied(Syy, Syx, Sxx, Beta, SigmaHat))) {
    stop("Give either individual-level data (X, Y, U) or summary-level matrices, not both.")
  }
  if (!is.null(Beta) && any(supplied(Syy, Syx))) stop("With the input (Sxx, Beta, SigmaHat), do not supply Syy or Syx.")
  if (xor(is.null(Beta), is.null(SigmaHat))) stop("Beta and SigmaHat must be given together.")
  if (is.null(Y) && is.null(Beta) && !is.null(Syy) && xor(is.null(Syx), is.null(Sxx))) stop("Syx and Sxx must be given together.")
  if (!missing(D) && is.logical(D) && is.matrix(D)) D = D * 1
  mode = if ((!is.null(Y) && is.null(X)) || (!is.null(Syy) && is.null(Syx) && is.null(Sxx) && is.null(Beta))) {
    "Y"
  } else if (!is.null(X) && !is.null(Y)) {
    "XY"
  } else if (!is.null(Syy) && !is.null(Syx) && !is.null(Sxx)) {
    "S"
  } else if (!is.null(Sxx) && !is.null(Beta) && !is.null(SigmaHat)) {
    "Beta"
  } else {
    stop("Please give (Y) or (Syy) or (X, Y) or (Syy, Syx, Sxx) or (Sxx, Beta, SigmaHat) as input.")
  }
  instruments = mode != "Y"

  # ---------------------------------------------------------------------------------------------------------------
  # 2. Sufficient statistics of the responses and instruments
  # ---------------------------------------------------------------------------------------------------------------
  if (mode == "Y") {

    if (!is.null(Y)) {
      check_matrix(Y, "Y")
      n = nrow(Y); p = ncol(Y)
      check_centered(Y, "Y")
      Syy = crossprod(Y) / n
    } else {
      check_n(n)
      check_matrix(Syy, "Syy")
      p = ncol(Syy)
      if (nrow(Syy) != p) stop("Syy should be a square matrix.")
    }
    k = 0

  } else if (mode == "XY") {

    check_matrix(X, "X"); check_matrix(Y, "Y")
    n = nrow(X)
    if (nrow(Y) != n) stop("Number of rows of X and Y should be equal.")
    p = ncol(Y); k = ncol(X)
    check_centered(X, "X"); check_centered(Y, "Y")
    Syy = crossprod(Y) / n
    Syx = crossprod(Y, X) / n
    Sxx = crossprod(X) / n

  } else if (mode == "S") {

    check_n(n)
    check_matrix(Syy, "Syy"); check_matrix(Syx, "Syx"); check_matrix(Sxx, "Sxx")
    p = ncol(Syy); k = ncol(Sxx)
    if (nrow(Syy) != p) stop("Syy should be a square matrix.")
    if (nrow(Sxx) != k) stop("Sxx should be a square matrix.")
    if (nrow(Syx) != p || ncol(Syx) != k) stop("Syx should have as many rows as Syy and as many columns as Sxx.")

  } else {   # mode == "Beta"

    check_n(n)
    if (SigmaStarModel != "diagonal") {
      stop("With the input (Sxx, Beta, SigmaHat) the covariances between the responses are reconstructed under ",
           "independent errors, so unmeasured confounding cannot be modeled: use SigmaStarModel = \"diagonal\", or give ",
           "individual-level data or the second-moment matrices (Syy, Syx, Sxx).")
    }
    if (!is.null(U) || any(supplied(Syu, Sxu, Suu))) {
      stop("Covariates cannot be combined with the input (Sxx, Beta, SigmaHat): adjust the marginal regressions for the ",
           "covariates instead, or give the second-moment matrices including Syu, Sxu and Suu.")
    }
    check_matrix(Sxx, "Sxx"); check_matrix(Beta, "Beta"); check_matrix(SigmaHat, "SigmaHat")
    p = nrow(Beta); k = ncol(Sxx)
    if (nrow(Sxx) != k) stop("Sxx should be a square matrix.")
    if (ncol(Beta) != k) stop("Number of columns of Beta should be equal to number of columns of Sxx.")
    if (nrow(SigmaHat) != p || ncol(SigmaHat) != k) stop("SigmaHat should have the dimensions of Beta.")
    Sxx = symmetrize(Sxx, "Sxx")
    check_D(D, p, k)

    # The reconstruction of Syy from marginal regressions needs, for every response, an instrument that
    # affects only that response (the first such instrument is used).
    valid_IV = apply(D, 1, function(row) any(row == 1 & colSums(D) == 1))
    if (!all(valid_IV)) {
      stop("With (Sxx, Beta, SigmaHat) as input every response needs an instrument that affects only that response (a column of D with a single 1).")
    }
    first_valid_IV = apply(D, 1, function(row) which(row == 1 & colSums(D) == 1)[1])

    # Syx from the marginal coefficients
    Syx = t(t(Beta) * diag(Sxx))

    # A from the ratios of marginal coefficients (indirect least squares on the valid instruments)
    # (I - A) Beta_New = diag(b): row i of the inverse of Beta_New, divided by its diagonal entry, is row i of I - A
    Beta_New = matrix(Beta[, first_valid_IV], nrow = p, ncol = p)
    Binv = tryCatch(solve(Beta_New), error = function(e) stop("The marginal coefficients of the instruments that affect ",
                    "only one response form a singular matrix: the causal effects cannot be reconstructed from Beta.", call. = FALSE))
    AEst0 = diag(p) - Binv / diag(Binv)
    diag(AEst0) = 0
    Mult_Mat = solve(diag(p) - AEst0)

    # Second moments of the responses and the error variances implied by the marginal regressions
    Syy_Diag = SigmaHat[, 1] + 2 * Beta[, 1] * Syx[, 1] - Beta[, 1]^2 * Sxx[1, 1]
    Cor_X = stats::cov2cor(Sxx)
    Cor_YX = t(t(Beta * (1 / sqrt(Syy_Diag))) * sqrt(diag(Sxx)))
    Cor_X_Det = det(Cor_X)
    Error_Sumsq = rep(0, p)
    for (i in 1:p) {
      R = cbind(c(1, Cor_YX[i, ]), rbind(Cor_YX[i, ], Cor_X))
      Error_Sumsq[i] = Syy_Diag[i] * det(R) / Cor_X_Det
    }
    Sigma0 = solve(Mult_Mat^2, Error_Sumsq)
    if (any(!is.finite(AEst0)) || any(!is.finite(Sigma0)) || any(Sigma0 <= 0)) {
      stop("Beta and SigmaHat are not consistent with the model: the implied error variances are not positive. ",
           "Check that SigmaHat holds the residual sums of squares divided by n of the same marginal regressions.")
    }
    Syy = Beta %*% Sxx %*% t(Beta) + Mult_Mat %*% (t(Mult_Mat) * Sigma0)

  }

  if (p < 2) stop("RGM() needs at least two responses.")
  if (!instruments && (!missing(D) || !is.null(selectB))) warning("No instruments were given: D and selectB are not used.", call. = FALSE)

  # ---------------------------------------------------------------------------------------------------------------
  # 3. Covariates: either U (individual-level input only) or the complete set of summary matrices
  # ---------------------------------------------------------------------------------------------------------------
  any_cov_summary = !is.null(Syu) || !is.null(Suu) || !is.null(Sxu)
  covariates = !is.null(U) || any_cov_summary
  l = 0
  if (covariates) {

    if (!is.null(U) && any_cov_summary) stop("Give the covariates either as U or as summary matrices (Syu, Sxu, Suu), not both.")

    if (!is.null(U)) {

      if (is.null(Y)) {
        stop("The covariate summary matrices cannot be computed from summary-level data and U. ",
             "Please provide (Y, U) or (X, Y, U), or the summary matrices Syu", if (instruments) ", Sxu" else "", " and Suu.")
      }
      check_matrix(U, "U")
      if (nrow(U) != n) stop("Number of rows of U and Y should be equal.")
      l = ncol(U)
      check_centered(U, "U")
      Syu = crossprod(Y, U) / n
      Suu = crossprod(U) / n
      if (instruments) Sxu = crossprod(X, U) / n

    } else {

      needed = c("Syu", "Suu", if (instruments) "Sxu")
      given = c(Syu = !is.null(Syu), Suu = !is.null(Suu), Sxu = !is.null(Sxu))[needed]
      if (!all(given)) {
        stop("Covariate summary matrices are incomplete: ", paste(needed, collapse = ", "), " are all required (missing: ",
             paste(needed[!given], collapse = ", "), ").")
      }
      check_matrix(Syu, "Syu"); check_matrix(Suu, "Suu")
      l = ncol(Suu)
      if (nrow(Suu) != l) stop("Suu should be a square matrix.")
      if (nrow(Syu) != p || ncol(Syu) != l) stop("Syu should have p rows and as many columns as Suu.")
      if (instruments) {
        check_matrix(Sxu, "Sxu")
        if (nrow(Sxu) != k || ncol(Sxu) != l) stop("Sxu should have k rows and as many columns as Suu.")
      } else {
        Sxu = NULL
      }

    }

  }

  # ---------------------------------------------------------------------------------------------------------------
  # 4. Symmetry and positive semi-definiteness of the second-moment matrices (the samplers assume both)
  # ---------------------------------------------------------------------------------------------------------------
  Syy = symmetrize(Syy, "Syy")
  if (instruments && mode != "Beta") Sxx = symmetrize(Sxx, "Sxx")   # Beta input: already done in step 2
  if (covariates) Suu = symmetrize(Suu, "Suu")

  # ---------------------------------------------------------------------------------------------------------------
  # 5. Instrument map D and whether the instrument effects are selected
  # ---------------------------------------------------------------------------------------------------------------
  if (!is.null(selectB) && (!is.logical(selectB) || length(selectB) != 1 || is.na(selectB))) {
    stop("selectB should be TRUE, FALSE or NULL.")
  }
  if (instruments) {

    if (mode != "Beta") check_D(D, p, k)   # Beta input: already checked in step 2

    if (is.null(selectB)) selectB = all(D == 1)

    # Without selection, identification of A rests on the excluded instruments: warn when the sufficient
    # condition (an instrument that affects only that response, for every response) fails.
    if (!selectB) {
      valid_IV = apply(D, 1, function(row) any(row == 1 & colSums(D) == 1))
      if (!all(valid_IV)) {
        warning("No instrument affects only response(s) ", paste(which(!valid_IV), collapse = ", "), ", so the causal effects of ",
                "these responses on the other responses may not be identified. Check the rank condition on the excluded instruments.",
                call. = FALSE)
      }
    }

  } else {
    selectB = NA
  }

  # ---------------------------------------------------------------------------------------------------------------
  # 6. Hyperparameters and MCMC settings
  # ---------------------------------------------------------------------------------------------------------------
  check_positive(list(aRho = aRho, bRho = bRho, nu1 = nu1, aSigma = aSigma, bSigma = bSigma,
                      tauC = tauC, omega1 = omega1, omega2 = omega2, lambda = lambda))
  if (instruments) check_positive(list(aPsi = aPsi, bPsi = bPsi, nu2 = nu2))
  if (nu1 >= 1) warning("nu1 is the spike variance relative to the slab variance and should be smaller than 1.", call. = FALSE)
  if (instruments && nu2 >= 1) warning("nu2 is the spike variance relative to the slab variance and should be smaller than 1.", call. = FALSE)
  if (!is.null(PropVarA) && length(PropVarA) != 1) stop("PropVarA should be NULL or a single positive number.")
  if (!is.null(PropVarB) && length(PropVarB) != 1) stop("PropVarB should be NULL or a single positive number.")
  if (!is.null(PropVarA)) check_positive(list(PropVarA = PropVarA))
  if (!is.null(PropVarB)) check_positive(list(PropVarB = PropVarB))
  if (!is.logical(gibbsB) || length(gibbsB) != 1 || is.na(gibbsB)) stop("gibbsB should be TRUE or FALSE.")
  if (!is.numeric(varpi) || length(varpi) != 1 || !is.finite(varpi) || varpi <= 0 || varpi >= 1) stop("varpi should be a number in (0, 1).")
  if (SigmaStarModel == "SSSL" && omega1 <= omega2) stop("omega1 (slab variance) should be larger than omega2 (spike variance).")
  if (instruments && !is.null(PropVarB) && gibbsB && (prior == "Spike and Slab" || isFALSE(selectB))) {
    warning("PropVarB is not used: the instrument effects B are drawn by Gibbs (gibbsB = TRUE).", call. = FALSE)
  }

  is_count = function(x, min) is.numeric(x) && length(x) == 1 && is.finite(x) && x == round(x) && x >= min && x <= .Machine$integer.max
  if (!is_count(nIter, 1)) stop("nIter should be a positive integer.")
  if (!is_count(nBurnin, 0) || nBurnin >= nIter) stop("nBurnin should be a non-negative integer strictly smaller than nIter.")
  if (!is_count(Thin, 1) || Thin > (nIter - nBurnin)) stop("Thin should be a positive integer not larger than nIter - nBurnin.")

  # Warn if the SSSL hyperparameters do not match the scale of the responses
  check_sssl_scale(Syy, SigmaStarModel)

  # ---------------------------------------------------------------------------------------------------------------
  # 7. Run the sampler
  # ---------------------------------------------------------------------------------------------------------------
  spike = prior == "Spike and Slab"
  sampler = paste0("RGM_", if (spike) "SpikeSlab" else "Threshold", "_",
                   c(SSSL = "SSSL", IW = "IW", diagonal = "Diag")[[SigmaStarModel]],
                   if (!instruments) "_Star" else "", if (covariates) "_Covariates" else "")

  # PropVarA / PropVarB = NULL: automatic starting proposal variances (signalled to the samplers by -1)
  args = list(S_YY = Syy, n = n, nIter = nIter, nBurnin = nBurnin, Thin = Thin, nu_1 = nu1,
              Prop_VarA = if (is.null(PropVarA)) -1 else PropVarA)
  if (spike) args = c(args, list(a_rho = aRho, b_rho = bRho))
  if (instruments) {
    args = c(args, list(S_YX = Syx, S_XX = Sxx, D = D, nu_2 = nu2, Prop_VarB = if (is.null(PropVarB)) -1 else PropVarB,
                        selectB = selectB, gibbsB = gibbsB))
    if (spike) args = c(args, list(a_psi = aPsi, b_psi = bPsi))
  }
  if (covariates) {
    args = c(args, list(S_YU = Syu, S_UU = Suu, TAU = tauC))
    if (instruments) args = c(args, list(S_XU = Sxu))
  }
  if (SigmaStarModel == "diagonal") args = c(args, list(a_sigma = aSigma, b_sigma = bSigma))
  if (SigmaStarModel == "SSSL") args = c(args, list(v0_square = omega2, v1_square = omega1, Pi = varpi, lambda = lambda))

  Output = do.call(sampler, args)

  # ---------------------------------------------------------------------------------------------------------------
  # 8. Assemble the result
  # ---------------------------------------------------------------------------------------------------------------
  name_map = c(A_Est = "AEst", B_Est = "BEst", C_Est = "CEst", zA_Est = "zAEst", zB_Est = "zBEst",
               A0_Est = "A0Est", B0_Est = "B0Est", Gamma_Est = "GammaEst", Tau_Est = "TauEst", Rho_Est = "RhoEst",
               Phi_Est = "PhiEst", Eta_Est = "EtaEst", Psi_Est = "PsiEst", tA_Est = "tAEst", tB_Est = "tBEst",
               Sigma_Est = "SigmaEst", Z_Est = "ZEst", AccptA = "AccptA", AccptB = "AccptB",
               Accpt_tA = "AccpttA", Accpt_tB = "AccpttB", PropVarA = "PropVarA", PropVarB = "PropVarB",
               LL_Pst = "LLPst", Gamma_Pst = "GammaPst")
  out = list()
  for (nm in names(name_map)) if (!is.null(Output[[nm]])) out[[name_map[[nm]]]] = Output[[nm]]

  # The diagonal model returns the error variances as a 1 x p matrix: report a vector; the log-likelihood trace too
  if (SigmaStarModel == "diagonal") out$SigmaEst = as.numeric(out$SigmaEst)
  if (!is.null(out$LLPst)) out$LLPst = as.numeric(out$LLPst)

  # Names of the responses, instruments and covariates, taken from the input (column names of Y, X, U or the
  # dimnames of the summary matrices); without names in the input the outputs have no dimnames
  first_names = function(len, ...) {
    for (x in list(...)) if (!is.null(x) && length(x) == len) return(as.character(x))
    NULL
  }
  yn = first_names(p, colnames(Y), colnames(Syy), rownames(Syy), rownames(Syx), rownames(Beta), if (instruments) rownames(D))
  xn = if (instruments) first_names(k, colnames(X), colnames(Syx), colnames(Sxx), colnames(Beta), colnames(D)) else NULL
  un = if (covariates) first_names(l, colnames(U), colnames(Syu), colnames(Suu), colnames(Sxu)) else NULL
  if (!is.null(yn)) {
    for (nm in c("AEst", "zAEst", "A0Est", "GammaEst", "TauEst", "RhoEst", "ZEst", "PropVarA")) {
      if (is.matrix(out[[nm]]) && all(dim(out[[nm]]) == p)) dimnames(out[[nm]]) = list(yn, yn)
    }
    if (is.matrix(out$SigmaEst)) dimnames(out$SigmaEst) = list(yn, yn) else names(out$SigmaEst) = yn
    if (!is.null(out$GammaPst)) dimnames(out$GammaPst) = list(yn, yn, NULL)
  }
  if (!is.null(yn) || !is.null(xn)) {
    for (nm in c("BEst", "zBEst", "B0Est", "PhiEst", "EtaEst", "PsiEst", "PropVarB")) {
      if (is.matrix(out[[nm]]) && nrow(out[[nm]]) == p && ncol(out[[nm]]) == k) dimnames(out[[nm]]) = list(yn, xn)
    }
  }
  if ((!is.null(yn) || !is.null(un)) && is.matrix(out$CEst) && nrow(out$CEst) == p && ncol(out$CEst) == l) {
    dimnames(out$CEst) = list(yn, un)
  }

  # Graph among the responses, placed after the adjacency matrices
  out = append(out, list(Graph = create_causal_graph_Y(out$AEst, out$zAEst)), after = max(which(names(out) %in% c("zAEst", "zBEst"))))

  # B drawn exactly by Gibbs: no acceptance rate or proposal variance for B
  B_gibbs = instruments && gibbsB && (spike || !selectB)
  if (B_gibbs) {
    out$AccptB = NULL
    out$PropVarB = NULL
  }

  # Without selection of the instrument effects there is nothing to report for psi and the threshold tB
  if (instruments && !selectB) {
    out$PsiEst = NULL
    out$tBEst = NULL
    out$AccpttB = NULL
  }

  out$settings = list(input = mode, prior = prior, SigmaStarModel = SigmaStarModel, instruments = instruments, covariates = covariates,
                      selectB = selectB, n = n, p = p, k = k, l = l, nIter = nIter, nBurnin = nBurnin, Thin = Thin,
                      aRho = aRho, bRho = bRho, nu1 = nu1, aPsi = aPsi, bPsi = bPsi, nu2 = nu2, aSigma = aSigma, bSigma = bSigma,
                      PropVarA = if (is.null(PropVarA)) "auto" else PropVarA,
                      PropVarB = if (!instruments || B_gibbs) NA else if (is.null(PropVarB)) "auto" else PropVarB,
                      gibbsB = if (instruments) B_gibbs else NA, tauC = tauC, omega1 = omega1, omega2 = omega2, varpi = varpi,
                      lambda = lambda, package_version = as.character(utils::packageVersion("MR.RGM")))

  finish_fit(out, Syy, Syx = if (instruments) Syx else NULL, Sxx = if (instruments) Sxx else NULL, D = if (instruments) D else NULL)

}



# ------------------------------------------------------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------------------------------------------------------

# A numeric matrix without missing values
check_matrix = function(M, name) {
  if (!is.numeric(M) || !is.matrix(M)) stop(name, " should be a numeric matrix.")
  if (anyNA(M) || any(!is.finite(M))) stop(name, " contains missing or non-finite values.")
  invisible(NULL)
}

# n: a positive integer (summary-level input)
check_n = function(n) {
  if (missing(n) || !is.numeric(n) || length(n) != 1 || n != round(n) || n <= 0) stop("n (number of observations) should be a positive integer.")
  invisible(NULL)
}

# D: binary p x k
check_D = function(D, p, k) {
  if (missing(D) || !is.numeric(D) || !is.matrix(D) || nrow(D) != p || ncol(D) != k || !all(D %in% c(0, 1))) {
    stop("D must be a numeric p x k matrix with entries 0 or 1 (rows = responses, columns = instruments).")
  }
  if (sum(D) == 0) stop("D has no entry equal to 1: no instrument is allowed to affect any response.")
  invisible(NULL)
}

# Strictly positive scalars (x: a named list, one element per argument)
check_positive = function(x) {
  bad = vapply(x, function(v) !is.numeric(v) || length(v) != 1 || !is.finite(v) || v <= 0, logical(1))
  if (any(bad)) stop(paste(names(x)[bad], collapse = ", "), " should be a single positive number.")
  invisible(NULL)
}

# The model has no intercept: warn when a column mean is not negligible relative to its spread
check_centered = function(M, name, tol = 0.1) {
  s = apply(M, 2, stats::sd)
  m = abs(colMeans(M))
  if (any(s > 0 & m > tol * s)) {
    warning("The model has no intercept, but column(s) ", paste(which(s > 0 & m > tol * s), collapse = ", "), " of ", name,
            " are not centered. Center (and scale) the columns before fitting, for example with scale().", call. = FALSE)
  }
  invisible(NULL)
}

# Second-moment matrices must be symmetric (the incremental likelihood updates of the samplers assume it);
# warn when the input is clearly asymmetric or has a negative eigenvalue, and use the symmetric part.
symmetrize = function(S, name) {
  asym = max(abs(S - t(S)))
  scale_S = max(abs(S), 1e-12)
  if (asym > 1e-6 * scale_S) warning(name, " is not symmetric (largest asymmetry ", signif(asym, 3), "); its symmetric part is used.", call. = FALSE)
  S = (S + t(S)) / 2
  if (ncol(S) <= 500) {   # the eigenvalue check is skipped for very large instrument sets
    ev = min(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
    if (ev < -1e-8 * scale_S) warning(name, " is not positive semi-definite (smallest eigenvalue ", signif(ev, 3), ").", call. = FALSE)
  }
  S
}

# Internal: the SSSL prior on the error covariance uses fixed spike and slab variances
# (defaults 0.005 and 5) and an exponential prior on the variances, which assume that the
# responses are on roughly unit scale. On very different scales the data cannot move the
# inclusion probabilities of the error covariances away from their prior value.
check_sssl_scale = function(Syy, SigmaStarModel) {

  if (SigmaStarModel == "SSSL") {

    # Second moments of the responses
    v = diag(Syy)

    if (any(v > 100) || any(v < 0.01)) {

      warning("SigmaStarModel = \"SSSL\" assumes responses on roughly unit scale (second moments of the responses range from ",
              signif(min(v), 3), " to ", signif(max(v), 3), "). Center and scale the responses (and instruments/covariates) ",
              "before fitting; otherwise the inclusion probabilities of the error covariances (ZEst) stay close to their prior value.",
              call. = FALSE)

    }

  }

  invisible(NULL)

}



# Internal: post-fit check for weak identification. If the estimated error variances are far larger
# than the second moments of the responses, the likelihood is nearly flat along a direction of the
# parameter space (typically MR.RGM+ with instruments that explain little variance, or responses
# without any instrument), and the posterior mostly reflects the prior.
finish_fit = function(out, Syy, Syx = NULL, Sxx = NULL, D = NULL) {

  warn = function(...) warning(..., call. = FALSE)
  st = out$settings

  # 1. Error variances much larger than the second moments of the responses: weak identification
  v = if (is.matrix(out$SigmaEst) && nrow(out$SigmaEst) == ncol(out$SigmaEst)) diag(out$SigmaEst) else as.numeric(out$SigmaEst)
  ratio = v / diag(Syy)
  out$diagnostics$errorVarianceRatio = ratio
  if (any(is.finite(ratio) & ratio > 5)) {
    warn("Estimated error variances exceed the second moments of the responses by a factor of up to ",
         signif(max(ratio, na.rm = TRUE), 2), " (diagnostics$errorVarianceRatio): the model may be weakly identified",
         if (isTRUE(st$instruments)) paste0(", for example with D all ones (selectB = TRUE) and instruments that explain little variance, ",
                                            "or with a response without any instrument. Check the posterior inclusion probabilities against ",
                                            "their prior, compare the two halves of GammaPst, and consider a known instrument map D (selectB = FALSE).")
         else ". Without instruments the causal effects are identified only through the prior; check the posterior inclusion probabilities against their prior.")
  }

  # 2. Weak instruments (known map, MR.RGM): first-stage F statistic of each response on its allowed instruments
  if (!is.null(D) && isFALSE(st$selectB)) {
    Fst = vapply(seq_len(nrow(D)), function(j) {
      L = which(D[j, ] == 1); kL = length(L); df2 = st$n - kL - 1
      if (kL == 0 || df2 <= 0) return(NA_real_)
      s = Syx[j, L]; G = Sxx[L, L, drop = FALSE]
      R2 = tryCatch(as.numeric(crossprod(s, solve(G, s))) / Syy[j, j], error = function(e) NA_real_)
      if (!is.finite(R2) || R2 >= 1) return(NA_real_)
      (R2 / kL) / ((1 - R2) / df2)
    }, numeric(1))
    out$diagnostics$instrumentF = Fst
    if (any(is.finite(Fst) & Fst < 10)) {
      warn("Weak instruments: the first-stage F statistic on the allowed instruments is below 10 for response(s) ",
           paste(which(is.finite(Fst) & Fst < 10), collapse = ", "), " (minimum ", signif(min(Fst, na.rm = TRUE), 3),
           "; diagnostics$instrumentF). These instruments identify the causal effects of these responses on the other responses, ",
           "which are therefore weakly identified: interpret them with caution, or use stronger or additional instruments.")
    }
  }

  # 3. Stability of the estimated network: with spectral radius of AEst < 1 the feedback loops converge (the
  #    equilibrium interpretation of the model); a value >= 1 is possible in an identified model, but with instrument
  #    selection (MR.RGM+) it usually signals weak identification or an observationally equivalent solution
  rho_A = max(Mod(eigen(out$AEst, only.values = TRUE)$values))
  out$diagnostics$spectralRadiusA = rho_A
  if (is.finite(rho_A) && rho_A >= 1) {
    warn("The spectral radius of the estimated causal effect matrix AEst is ", signif(rho_A, 3), " (>= 1): the estimated ",
         "feedback loops do not converge, so the equilibrium interpretation of the model does not hold.",
         if (isTRUE(st$selectB)) paste0(" With instrument selection (selectB = TRUE) this usually signals weak identification or ",
                                        "an observationally equivalent solution: consider a known instrument map D (selectB = FALSE) or more data.")
         else " Check diagnostics$instrumentF and diagnostics$errorVarianceRatio, and run a longer chain or several chains.")
  }

  # 4. Convergence: causal inclusion probabilities in the two halves of the kept draws. A difference is flagged only
  #    when it exceeds 0.2 and 4 times its Monte Carlo standard error (batch means: 10 batches per half), so that
  #    ordinary autocorrelation of uncertain edges does not trigger the warning.
  GP = out$GammaPst
  if (!is.null(GP) && dim(GP)[3] >= 200) {
    m = dim(GP)[3]; h = m %/% 2; off = which(row(out$GammaEst) != col(out$GammaEst))
    half_stats = function(idx) {
      x = matrix(GP, nrow = length(out$GammaEst))[off, idx, drop = FALSE]          # edges x draws
      bm = sapply(split(seq_along(idx), cut(seq_along(idx), 10, labels = FALSE)), function(b) rowMeans(x[, b, drop = FALSE]))
      list(mean = rowMeans(x), se = apply(matrix(bm, nrow = length(off)), 1, stats::sd) / sqrt(10))
    }
    a = half_stats(1:h); b = half_stats((h + 1):m)
    d = abs(a$mean - b$mean); se = sqrt(a$se^2 + b$se^2)
    out$diagnostics$halvesMaxDiff = max(d)
    out$diagnostics$halvesMCSE = se[which.max(d)]
    flag = d > 0.2 & d > 4 * se
    if (any(flag)) {
      warn("The posterior inclusion probabilities of ", sum(flag), " causal edge(s) differ by up to ", signif(max(d[flag]), 2),
           " between the two halves of the kept draws, beyond their Monte Carlo error: the chain may not have converged. ",
           "Increase nIter and nBurnin, or run several chains with different seeds and compare them.")
    }
  }

  # 5. Random-walk updates of A: acceptance rate after the (adaptive) burn-in
  if (!is.null(out$AccptA)) {
    out$diagnostics$acceptanceA = out$AccptA
    if (is.finite(out$AccptA) && out$AccptA < 10) {
      warn("The acceptance rate of the Metropolis updates of A after burn-in is ", signif(out$AccptA, 2), "%: the ",
           "proposal scales have not adapted well. Increase nBurnin (the proposals are tuned during the burn-in).")
    }
  }

  out

}
