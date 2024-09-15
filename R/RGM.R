#' Fitting Bayesian Multivariate Bidirectional Mendelian Randomization Networks
#'
#' @description The RGM function transforms causal inference by merging Mendelian randomization and network-based methods, enabling the creation of comprehensive causal graphs within complex biological systems. RGM accommodates varied data contexts with three input options: individual-level data (X, Y matrices), summary-level data including Syy, Syx, and Sxx matrices, and intricate data with challenging cross-correlations, utilizing Sxx, Beta, and SigmaHat matrices.
#'              For the latter input, data centralization is necessary. Users can select any of these data formats to suit their needs and don’t have to specify all of them, allowing flexibility based on data availability. Crucial inputs encompass "d" (instrument count per response) and "n" (total observations, only required for summary level data), amplified by customizable parameters that refine analysis. Additionally, users can tailor the analysis by setting parameters such as "nIter" (number of MCMC iterations), "nBurnin" (number of discarded samples during burn-in for convergence), and "Thin" (thinning of posterior samples). These customizable parameters enhance the precision and relevance of the analysis.
#'              RGM provides essential causal effect/strength estimates between response variables and between response and instrument variables. Moreover, it furnishes adjacency matrices, visually mapping causal graph structures. These outputs empower researchers to untangle intricate relationships within biological networks, fostering a holistic understanding of complex systems.
#'
#' @param X A matrix of dimension n * k. Each row represents a distinct observation, and each column corresponds to a specific instrumental variable. The default value is set to NULL.
#' @param Y A matrix of dimension n * p. Each row represents a specific observation, and each column corresponds to a particular response variable. The default value is set to NULL.
#' @param Syy A matrix of dimension p * p, where "p" is the number of response variables. It is calculated as t(Y) %*% Y / n, where "Y" represents the response data matrix and "n" is the number of observations.
#' @param Syx A matrix of dimension p * k, where "p" is the number of response variables, and "k" is the number of instrumental variables. It is calculated as t(Y) %*% X / n, where "Y" represents the response data matrix, "X" represents the instrumental data matrix, and "n" is the number of observations.
#' @param Sxx A matrix of dimension k * k, where "k" is the number of instrumental variables. It is derived as t(X) %*% X / n, where "X" represents the instrumental data matrix and "n" is the number of observations.
#' @param Beta A matrix of dimension p * k, where each row corresponds to a specific response variable and each column pertains to an instrumental variable. Each entry represents the regression coefficient of the response variable on the instrumental variable. When using Beta as input, ensure that both Y (response data) and X (instrument data) are centered before calculating Beta, Sxx, and SigmaHat.
#' @param SigmaHat A matrix of dimension p * k. Each row corresponds to a specific response variable, and each column pertains to an instrumental variable. Each entry represents the mean square error of the regression between the response and the instrumental variable. As with Beta, ensure that both Y and X are centered before calculating SigmaHat.
#' @param D A binary indicator matrix of dimension p * k, where each row corresponds to a response variable, and each column corresponds to an instrumental variable. The entry `D[i, j]` is 1 if instrumental variable j affects response variable i, and 0 otherwise. For each response variable, there must be at least one instrumental variable that affects only that response (i.e., for each row in D, there must be at least one column with 1, and that column must have zeros in all other rows). If you use `Syy`, `Beta`, and `SigmaHat` as inputs, this condition must be satisfied to run this algorithm. If this condition is not met, an error will be thrown. However, if using `X`, `Y` or `Syy`, `Syx`, `Sxx` as inputs, a warning will be issued if the condition is violated, but the method will still proceed.
#' @param n A positive integer input representing the count of data points or observations in the dataset. This input is only required when summary level data is used as input.
#' @param nIter A positive integer input representing the number of MCMC (Markov Chain Monte Carlo) sampling iterations. The default value is set to 10,000.
#' @param nBurnin A non-negative integer input representing the number of samples to be discarded during the burn-in phase of MCMC sampling. It's important that nBurnin is less than nIter. The default value is set to 2000.
#' @param Thin A positive integer input denoting the thinning factor applied to posterior samples. Thinning reduces the number of samples retained from the MCMC process for efficiency. Thin should not exceed (nIter - nBurnin). The default value is set to 1.
#' @param prior A parameter representing the prior assumption on the graph structure. It offers two options: "Threshold" or "Spike and Slab". The default value is "Spike and Slab".
#' @param aRho A positive scalar input representing the first parameter of a Beta distribution. The default value is set to 3.
#' @param bRho A positive scalar input representing the second parameter of a Beta distribution. The default value is set to 1.
#' @param nu1 A positive scalar input representing the multiplication factor in the variance of the spike part in the spike and slab distribution of matrix A. The default value is set to 0.001.
#' @param aPsi A positive scalar input corresponding to the first parameter of a Beta distribution. The default value is set to 0.5.
#' @param bPsi  A positive scalar input corresponding to the second parameter of a Beta distribution. The default value is set to 0.5.
#' @param nu2 A positive scalar input corresponding to the multiplication factor in the variance of the spike part in the spike and slab distribution of matrix B. The default value is set to 0.0001.
#' @param aSigma A positive scalar input corresponding to the first parameter of an Inverse Gamma distribution, which is associated with the variance of the model. The default value is set to 0.01.
#' @param bSigma A positive scalar input corresponding to the second parameter of an Inverse Gamma distribution, which is associated with the variance of the model. The default value is set to 0.01.
#' @param PropVarA A positive scalar input representing the variance of the normal distribution used for proposing terms within the A matrix. The default value is set to 0.01.
#' @param PropVarB A positive scalar input representing the variance of the normal distribution used for proposing terms within the B matrix. The default value is set to 0.01.
#'
#' @return
#'
#' \item{AEst}{A matrix of dimensions p * p, representing the estimated causal effects or strengths between the response variables.}
#' \item{BEst}{A matrix of dimensions p * k, representing the estimated causal effects or strengths of the instrument variables on the response variables. Each row corresponds to a specific response variable, and each column corresponds to a particular instrument variable.}
#' \item{zAEst}{A binary adjacency matrix of dimensions p * p, indicating the graph structure between the response variables. Each entry in the matrix represents the presence (1) or absence (0) of a causal link between the corresponding response variables.}
#' \item{zBEst}{A binary adjacency matrix of dimensions p * k, illustrating the graph structure between the response variables and the instrument variables. Each row corresponds to a specific response variable, and each column corresponds to a particular instrument variable. The presence of a causal link is denoted by 1, while the absence is denoted by 0.}
#' \item{A0Est}{A matrix of dimensions p * p, representing the estimated causal effects or strengths between response variables before thresholding. This output is particularly relevant for cases where the "Threshold" prior assumption is utilized.}
#' \item{B0Est}{A matrix of dimensions p * k, representing the estimated causal effects or strengths between the response variables and the instrument variables before thresholding. This output is particularly relevant for cases where the "Threshold" prior assumption is utilized. Each row corresponds to a specific response variable, and each column corresponds to a particular instrument variable.}
#' \item{GammaEst}{A matrix of dimensions p * p, representing the estimated probabilities of edges between response variables in the graph structure. Each entry in the matrix indicates the probability of a causal link between the corresponding response variables.}
#' \item{TauEst}{A matrix of dimensions p * p, representing the estimated variances of causal interactions between response variables. Each entry in the matrix corresponds to the variance of the causal effect between the corresponding response variables.}
#' \item{PhiEst}{A matrix of dimensions p * k, representing the estimated probabilities of edges between response and instrument variables in the graph structure. Each row corresponds to a specific response variable, and each column corresponds to a particular instrument variable.}
#' \item{EtaEst}{A matrix of dimensions p * k, representing the estimated variances of causal interactions between response and instrument variables. Each row corresponds to a specific response variable, and each column corresponds to a particular instrument variable.}
#' \item{tAEst}{A scalar value representing the estimated thresholding value of causal interactions between response variables. This output is relevant when using the "Threshold" prior assumption.}
#' \item{tBEst}{A scalar value representing the estimated thresholding value of causal interactions between response and instrument variables. This output is applicable when using the "Threshold" prior assumption.}
#' \item{SigmaEst}{A vector of length p, representing the estimated variances of each response variable. Each element in the vector corresponds to the variance of a specific response variable.}
#' \item{AccptA}{The percentage of accepted entries in the A matrix, which represents the causal interactions between response variables. This metric indicates the proportion of proposed changes that were accepted during the sampling process.}
#' \item{AccptB}{The percentage of accepted entries in the B matrix, which represents the causal interactions between response and instrument variables. This metric indicates the proportion of proposed changes that were accepted during the sampling process.}
#' \item{AccpttA}{The percentage of accepted thresholding values for causal interactions between response variables when using the "Threshold" prior assumption. This metric indicates the proportion of proposed thresholding values that were accepted during the sampling process.}
#' \item{AccpttB}{The percentage of accepted thresholding values for causal interactions between response and instrument variables when using the "Threshold" prior assumption. This metric indicates the proportion of proposed thresholding values that were accepted during the sampling process.}
#' \item{LLPst}{A vector containing the posterior log-likelihoods of the model. Each element in the vector represents the log-likelihood of the model given the observed data and the estimated parameters.}
#' \item{RhoEst}{A matrix of dimensions p * p, representing the estimated Bernoulli success probabilities of causal interactions between response variables when using the "Spike and Slab" prior assumption. Each entry in the matrix corresponds to the success probability of a causal interaction between the corresponding response variables.}
#' \item{PsiEst}{A matrix of dimensions p * k, representing the estimated Bernoulli success probabilities of causal interactions between response and instrument variables when using the "Spike and Slab" prior assumption. Each row in the matrix corresponds to a specific response variable, and each column corresponds to a particular instrument variable.}
#' \item{GammaPst}{An array containing the posterior samples of the network structure among the response variables.}
#'
#'
#'
#'
#'
#'
#'
#' @export
#'
#' @examples
#'
#'
#' # ---------------------------------------------------------
#'
#' # Example 1:
#' # Run RGM based on individual level data with Threshold prior based on the model Y = AY + BX + E
#'
#' # Data Generation
#' set.seed(9154)
#'
#' # Number of data points
#' n = 10000
#'
#' # Number of response variables and number of instrument variables
#' p = 3
#' k = 4
#'
#' # Initialize causal interaction matrix between response variables
#' A = matrix(sample(c(-0.1, 0.1), p^2, replace = TRUE), p, p)
#'
#' # Diagonal entries of A matrix will always be 0
#' diag(A) = 0
#'
#' # Make the network sparse
#' A[sample(which(A != 0), length(which(A != 0)) / 2)] = 0
#'
#' # Create D matrix (Indicator matrix where each row corresponds to a response variable
#' # and each column corresponds to an instrument variable)
#' D = matrix(0, nrow = p, ncol = k)
#'
#' # Manually assign values to D matrix
#' D[1, 1:2] = 1  # First response variable is influenced by the first 2 instruments
#' D[2, 3] = 1    # Second response variable is influenced by the 3rd instrument
#' D[3, 4] = 1    # Third response variable is influenced by the 4th instrument
#'
#'
#' # Initialize B matrix
#' B = matrix(0, p, k)  # Initialize B matrix with zeros
#'
#' # Calculate B matrix based on D matrix
#' for (i in 1:p) {
#'   for (j in 1:k) {
#'     if (D[i, j] == 1) {
#'       B[i, j] = 1  # Set B[i, j] to 1 if D[i, j] is 1
#'     }
#'   }
#' }
#'
#' # Define Sigma matrix
#' Sigma = diag(p)
#'
#' # Compute Mult_Mat
#' Mult_Mat = solve(diag(p) - A)
#'
#' # Calculate Variance
#' Variance = Mult_Mat %*% Sigma %*% t(Mult_Mat)
#'
#' # Generate instrument data matrix
#' X = matrix(rnorm(n * k, 0, 1), nrow = n, ncol = k)
#'
#' # Initialize response data matrix
#' Y = matrix(0, nrow = n, ncol = p)
#'
#' # Generate response data matrix based on instrument data matrix
#' for (i in 1:n) {
#'   Y[i, ] = MASS::mvrnorm(n = 1, Mult_Mat %*% B %*% X[i, ], Variance)
#' }
#'
#' # Define a function to create smaller arrowheads
#' smaller_arrowheads = function(graph) {
#'   igraph::E(graph)$arrow.size = 1  # Adjust the arrow size value as needed
#'   return(graph)
#' }
#'
#' # Print true causal interaction matrices between response variables
#' # and between response and instrument variables
#' A
#' B
#'
#' # Plot the true graph structure between response variables
#' plot(smaller_arrowheads(igraph::graph_from_adjacency_matrix((A != 0) * 1,
#'   mode = "directed")), layout = igraph::layout_in_circle, main = "True Graph")
#'
#' # Apply RGM on individual level data for Threshold Prior
#' Output = RGM(X = X, Y = Y, D = D, prior = "Threshold")
#'
#' # Get the graph structure between response variables
#' Output$zAEst
#'
#' # Plot the estimated graph structure between response variables
#' plot(smaller_arrowheads(igraph::graph_from_adjacency_matrix(Output$zAEst,
#'   mode = "directed")), layout = igraph::layout_in_circle, main = "Estimated Graph")
#'
#' # Get the estimated causal strength matrix between response variables
#' Output$AEst
#'
#' # Get the graph structure between response and instrument variables
#' Output$zBEst
#'
#' # Get the estimated causal strength matrix between response and instrument variables
#' Output$BEst
#'
#' # Plot posterior log-likelihood
#' plot(Output$LLPst, type = 'l', xlab = "Number of Iterations", ylab = "Log-likelihood")
#'
#' # -----------------------------------------------------------------
#' # Example 2:
#' # Run RGM based on Syy, Syx and Sxx with Spike and Slab prior based on the model Y = AY + BX + E
#'
#' # Data Generation
#' set.seed(9154)
#'
#' # Number of data points
#' n = 10000
#'
#' # Number of response variables and number of instrument variables
#' p = 3
#' k = 4
#'
#' # Initialize causal interaction matrix between response variables
#' A = matrix(sample(c(-0.1, 0.1), p^2, replace = TRUE), p, p)
#'
#' # Diagonal entries of A matrix will always be 0
#' diag(A) = 0
#'
#' # Make the network sparse
#' A[sample(which(A!=0), length(which(A!=0))/2)] = 0
#'
#'
#' # Create D matrix (Indicator matrix where each row corresponds to a response variable
#' # and each column corresponds to an instrument variable)
#' D = matrix(0, nrow = p, ncol = k)
#'
#' # Manually assign values to D matrix
#' D[1, 1:2] = 1  # First response variable is influenced by the first 2 instruments
#' D[2, 3] = 1    # Second response variable is influenced by the 3rd instrument
#' D[3, 4] = 1    # Third response variable is influenced by the 4th instrument
#'
#'
#' # Initialize B matrix
#' B = matrix(0, p, k)  # Initialize B matrix with zeros
#'
#' # Calculate B matrix based on D matrix
#' for (i in 1:p) {
#'   for (j in 1:k) {
#'     if (D[i, j] == 1) {
#'       B[i, j] = 1  # Set B[i, j] to 1 if D[i, j] is 1
#'     }
#'   }
#' }
#'
#'
#' Sigma = diag(p)
#'
#' Mult_Mat = solve(diag(p) - A)
#'
#' Variance = Mult_Mat %*% Sigma %*% t(Mult_Mat)
#'
#' # Generate instrument data matrix
#' X = matrix(rnorm(n * k, 0, 1), nrow = n, ncol = k)
#'
#' # Initialize response data matrix
#' Y = matrix(0, nrow = n, ncol = p)
#'
#' # Generate response data matrix based on instrument data matrix
#' for (i in 1:n) {
#'
#'     Y[i, ] = MASS::mvrnorm(n = 1, Mult_Mat %*% B %*% X[i, ], Variance)
#'
#' }
#'
#'
#' # Calculate summary level data
#' Syy = t(Y) %*% Y / n
#' Syx = t(Y) %*% X / n
#' Sxx = t(X) %*% X / n
#'
#'
#' # Print true causal interaction matrices between response variables
#' # and between response and instrument variables
#' A
#' B
#'
#' # Plot the true graph structure between response variables
#' plot(smaller_arrowheads(igraph::graph_from_adjacency_matrix(((A != 0) * 1),
#'  mode = "directed")), layout = igraph::layout_in_circle, main = "True Graph")
#'
#'
#' # Apply RGM on summary level data for Spike and Slab Prior
#' Output = RGM(Syy = Syy, Syx = Syx, Sxx = Sxx,
#'           D = D, n = 10000, prior = "Spike and Slab")
#'
#' # Get the graph structure between response variables
#' Output$zAEst
#'
#' # Plot the estimated graph structure between response variables
#' plot(smaller_arrowheads(igraph::graph_from_adjacency_matrix(Output$zAEst,
#'  mode = "directed")), layout = igraph::layout_in_circle, main = "Estimated Graph")
#'
#' # Get the estimated causal strength matrix between response variables
#' Output$AEst
#'
#' # Get the graph structure between response and instrument variables
#' Output$zBEst
#'
#' # Get the estimated causal strength matrix between response and instrument variables
#' Output$BEst
#'
#' # Plot posterior log-likelihood
#' plot(Output$LLPst, type = 'l', xlab = "Number of Iterations", ylab = "Log-likelihood")
#'
#'
#'
#'
#' # -----------------------------------------------------------------
#' # Example 3:
#' # Run RGM based on Sxx, Beta and SigmaHat with Spike and Slab prior
#' # based on the model Y = AY + BX + E
#'
#' # Data Generation
#' set.seed(9154)
#'
#' # Number of datapoints
#' n = 10000
#'
#' # Number of response variables and number of instrument variables
#' p = 3
#' k = 4
#'
#' # Initialize causal interaction matrix between response variables
#' A = matrix(sample(c(-0.1, 0.1), p^2, replace = TRUE), p, p)
#'
#' # Diagonal entries of A matrix will always be 0
#' diag(A) = 0
#'
#' # Make the network sparse
#' A[sample(which(A!=0), length(which(A!=0))/2)] = 0
#'
#'
#' # Create D matrix (Indicator matrix where each row corresponds to a response variable
#' # and each column corresponds to an instrument variable)
#' D = matrix(0, nrow = p, ncol = k)
#'
#' # Manually assign values to D matrix
#' D[1, 1:2] = 1  # First response variable is influenced by the first 2 instruments
#' D[2, 3] = 1    # Second response variable is influenced by the 3rd instrument
#' D[3, 4] = 1    # Third response variable is influenced by the 4th instrument
#'
#'
#' # Initialize B matrix
#' B = matrix(0, p, k)  # Initialize B matrix with zeros
#'
#' # Calculate B matrix based on D matrix
#' for (i in 1:p) {
#'   for (j in 1:k) {
#'     if (D[i, j] == 1) {
#'       B[i, j] = 1  # Set B[i, j] to 1 if D[i, j] is 1
#'     }
#'   }
#' }
#'
#'
#' Sigma = diag(p)
#'
#' Mult_Mat = solve(diag(p) - A)
#'
#' Variance = Mult_Mat %*% Sigma %*% t(Mult_Mat)
#'
#' # Generate DNA expressions
#' X = matrix(rnorm(n * k, 0, 1), nrow = n, ncol = k)
#'
#' # Initialize response data matrix
#' Y = matrix(0, nrow = n, ncol = p)
#'
#' # Generate response data matrix based on instrument data matrix
#' for (i in 1:n) {
#'
#'     Y[i, ] = MASS::mvrnorm(n = 1, Mult_Mat %*% B %*% X[i, ], Variance)
#'
#' }
#'
#'
#' # Centralize Data
#' Y = t(t(Y) - colMeans(Y))
#' X = t(t(X) - colMeans(X))
#'
#' # Calculate Sxx
#' Sxx = t(X) %*% X / n
#'
#' # Generate Beta matrix and SigmaHat
#' Beta = matrix(0, nrow = p, ncol = k)
#' SigmaHat = matrix(0, nrow = p, ncol = k)
#'
#' for (i in 1:p) {
#'
#'    for (j in 1:k) {
#'
#'        fit = lm(Y[, i] ~ X[, j])
#'
#'        Beta[i, j] =  fit$coefficients[2]
#'
#'        SigmaHat[i, j] = sum(fit$residuals^2) / n
#'
#'        }
#'
#'    }
#'
#'
#' # Print true causal interaction matrices between response variables
#' # and between response and instrument variables
#' A
#' B
#'
#'
#' # Plot the true graph structure between response variables
#' plot(smaller_arrowheads(igraph::graph_from_adjacency_matrix(((A != 0) * 1),
#'  mode = "directed")), layout = igraph::layout_in_circle, main = "True Graph")
#'
#'
#' # Apply RGM based on Sxx, Beta and SigmaHat for Spike and Slab Prior
#' Output = RGM(Sxx = Sxx, Beta = Beta, SigmaHat = SigmaHat,
#'           D = D, n = 10000, prior = "Spike and Slab")
#'
#' # Get the graph structure between response variables
#' Output$zAEst
#'
#' # Plot the estimated graph structure between response variables
#' plot(smaller_arrowheads(igraph::graph_from_adjacency_matrix(Output$zAEst,
#'  mode = "directed")), layout = igraph::layout_in_circle, main = "Estimated Graph")
#'
#' # Get the estimated causal strength matrix between response variables
#' Output$AEst
#'
#' # Get the graph structure between response and instrument variables
#' Output$zBEst
#'
#' # Get the estimated causal strength matrix between response and instrument variables
#' Output$BEst
#'
#' # Plot posterior log-likelihood
#' plot(Output$LLPst, type = 'l', xlab = "Number of Iterations", ylab = "Log-likelihood")
#'
#'
#'
#'
#'
#'
#' @references
#' Ni, Y., Ji, Y., & Müller, P. (2018).
#' Reciprocal graphical models for integrative gene regulatory network analysis.
#' \emph{Bayesian Analysis},
#' \strong{13(4)}, 1095-1110.
#' \doi{10.1214/17-BA1087}.
RGM = function(X = NULL, Y = NULL, Syy = NULL, Syx = NULL, Sxx = NULL, Beta = NULL, SigmaHat = NULL, D, n, nIter = 10000, nBurnin = 2000, Thin = 1, prior = c("Threshold", "Spike and Slab"), aRho = 3, bRho = 1, nu1 = 0.001, aPsi = 0.5, bPsi = 0.5, nu2 = 0.0001, aSigma = 0.01, bSigma = 0.01, PropVarA = 0.01, PropVarB = 0.01){

  # Check whether Y or Syy is given as data input
  if((!is.null(Y) && is.null(X)) || (!is.null(Syy) && is.null(Syx) && is.null(Sxx))){

    # Check whether input Y is given
    if(!is.null(Y)){

      # Check whether Y is a numeric matrix
      if(!is.numeric(Y) || !is.matrix(Y)){

        # Print an error message
        stop("Y should be numeric matrix.")

      }

      # Calculate number of datapoints
      n = nrow(Y)

      # Calculate number of response variables from Y matrix
      p = ncol(Y)

      # Calculate Syy
      Syy = t(Y) %*% Y / n

    } else {

      # Check whether n is a positive integer
      if(!is.numeric(n) || n != round(n) || n <= 0){

        # Print an error message
        stop("Number of datapoints should be a positive integer.")

      }

      # Check whether Syy is numeric matrix
      if(!is.numeric(Syy) || !is.matrix(Syy)){

        # Print an error message
        stop("Syy should be a numeric matrix.")

      }

      # Calculate number of response variables from Syy matrix
      p = ncol(Syy)

      # Check whether number of rows of Syy is equal to p
      if(nrow(Syy) != p){

        # Print an error message
        stop("Syy should be a square matrix.")

      }


    }


    # Check whether all the beta and inverse gamma parameters are positive or not
    if(!is.numeric(aRho) || aRho < 0 || !is.numeric(bRho) || bRho < 0  || !is.numeric(aSigma) || aSigma < 0 || !is.numeric(bSigma) || bSigma < 0){

      # Print an error message
      stop("All the beta and inverse gamma parameters should be positive.")

    }


    # Check whether the variance terms are positive or not
    if(!is.numeric(nu1) || nu1 <= 0 || !is.numeric(PropVarA) || PropVarA <= 0){

      # Print an error message
      stop("All the variance terms should be positive.")

    }


    # Check whether nIter is a postive integer
    if(!is.numeric(nIter) || nIter != round(nIter) || nIter <= 0){

      # Print an error message
      stop("Number of iterations should be a positive integer.")

    }

    # Check whether nBurnin is a non-negative integer and strictly less than nIter
    if(!is.numeric(nBurnin) || nBurnin != round(nBurnin) || nBurnin < 0 || nBurnin >= nIter){

      # Print an error message
      stop("Number of Burnin points should be a nonnegative integer strictly less than number of iterations.")

    }

    # Check whether Thin is a positive integer less than equal to (nIter - nBurnin)
    if(!is.numeric(Thin) || Thin != round(Thin) || Thin <= 0 || Thin > (nIter - nBurnin)){

      # Print an error message
      stop("Thin should be a positive integer less than or equal to the difference between number of iterations and number of burnin points.")

    }

    # Apply RGM for Spike and Slab prior and Threshold prior
    if ("Spike and Slab" %in% prior){

      # Run the algorithm for Spike and Slab prior
      Output = RGM_SpikeSlab1(Syy, n, nIter = nIter, nBurnin = nBurnin, Thin = Thin,
                              a_rho = aRho, b_rho = bRho, nu_1 = nu1,
                              a_sigma = aSigma, b_sigma = bSigma, Prop_VarA = PropVarA)



      # Return outputs
      return(list(AEst = Output$A_Est, zAEst = Output$zA_Est,
                  GammaEst = Output$Gamma_Est, TauEst = Output$Tau_Est, RhoEst = Output$Rho_Est,
                  SigmaEst = Output$Sigma_Est,
                  AccptA = Output$AccptA,
                  LLPst = Output$LL_Pst, GammaPst = Output$Gamma_Pst))



    } else if ("Threshold" %in% prior){

      # Run the algorithm for Threshold prior
      Output = RGM_Threshold1(Syy, n, nIter = nIter, nBurnin = nBurnin, Thin = Thin, nu_1 = nu1,
                              a_sigma = aSigma, b_sigma = bSigma, Prop_VarA = PropVarA)




      # Return outputs
      return(list(AEst = Output$A_Est, zAEst = Output$zA_Est,
                  A0Est = Output$A0_Est, GammaEst = Output$Gamma_Est, TauEst = Output$Tau_Est,
                  tAEst = Output$tA_Est,
                  SigmaEst = Output$Sigma_Est,
                  AccptA = Output$AccptA, AccpttA = Output$Accpt_tA,
                  LLPst = Output$LL_Pst, GammaPst = Output$Gamma_Pst))





    } else {

      # Print an error message
      stop("Please specify a prior among Threshold prior and Spike and Slab prior.")


    }


  } else if((!is.null(X) && !is.null(Y)) || (!is.null(Syy) && !is.null(Syx) && !is.null(Sxx)) || (!is.null(Sxx) && !is.null(Beta) && !is.null(SigmaHat))){


    # Check whether input X and Y are given
    if(!is.null(X) && !is.null(Y)){

      # Check whether X and Y are both numeric matrices
      if(!is.numeric(X) || !is.numeric(Y) || !is.matrix(X) || !is.matrix(Y)){

        # Print an error message
        stop("X and Y should be numeric matrices.")

      }

      # Calculate number of datapoints
      n = nrow(X)

      # Check whether number of rows of X and Y are equal or not
      if(nrow(Y) != n){

        # Print an error message
        stop("Number of rows of X and Y should be equal.")

      }

      # Calculate number of response variables from Y matrix
      p = ncol(Y)

      # Calculate number of instrument variables from X matrix
      k = ncol(X)

      # Calculate Syy, Syx, Sxx
      Syy = t(Y) %*% Y / n
      Syx = t(Y) %*% X / n
      Sxx = t(X) %*% X / n

    } else if(!is.null(Syy) && !is.null(Syx) && !is.null(Sxx)){

      # Check whether n is a positive integer
      if(!is.numeric(n) || n != round(n) || n <= 0){

        # Print an error message
        stop("Number of datapoints should be a positive integer.")

      }

      # Check whether Syy ,Syx and Sxx are numeric matrices
      if(!is.numeric(Syy) || !is.numeric(Syx) || !is.numeric(Sxx) || !is.matrix(Syy) || !is.matrix(Syx) || !is.matrix(Sxx)){

        # Print an error message
        stop("Syy, Syx and Sxx should be numeric matrices.")

      }

      # Calculate number of response variables from Syy matrix
      p = ncol(Syy)

      # Calculate number of instrument variables from Sxx matrix
      k = ncol(Sxx)

      # Check whether number of rows of Syy is equal to p
      if(nrow(Syy) != p){

        # Print an error message
        stop("Syy should be a square matrix.")

      }

      # Check whether number of rows of Sxx is equal to k
      if(nrow(Sxx) != k){

        # Print an error message
        stop("Sxx should be a square matrix.")

      }

      # Check whether number of rows of Syx is equal to p and number of columns of Syx is equal to k
      if(nrow(Syx) != p || ncol(Syx) != k){

        # Print an error message
        stop("Number of rows of Syx should be equal to number of rows of Syy and number of columns of Syx should be equal to number of columns of Sxx.")

      }



    } else {

      # Check whether n is a positive integer
      if(!is.numeric(n) || n != round(n) || n <= 0){

        # Print an error message
        stop("Number of datapoints should be a positive integer.")

      }

      # Check whether Sxx, Beta and SigmaHat are numeric matrices
      if(!is.numeric(Sxx) || !is.numeric(Beta) || !is.numeric(SigmaHat) || !is.matrix(Sxx) || !is.matrix(Beta) || !is.matrix(SigmaHat)){

        # Print an error message
        stop("Sxx, Beta and SigmaHat should be numeric matrices.")

      }

      # Calculate number of response variables from Beta matrix
      p = nrow(Beta)

      # Calculate number of instrument variables from Sxx matrix
      k = ncol(Sxx)

      # Check whether number of rows of Sxx is equal to k
      if(nrow(Sxx) != k){

        # Print an error message
        stop("Sxx should be a square matrix.")

      }

      # Check whether number of columns of Beta is equal to k
      if(ncol(Beta) != k){

        # Print an error message
        stop("Number of cloumns of Beta should be equal to number of columns of Sxx.")

      }

      # Check whether number of rows of SigmaHat is equal to p and number of columns of SigmaHat is equal to k
      if(nrow(SigmaHat) != p || ncol(SigmaHat) != k){

        # Print an error message
        stop("Number of rows of SigmaHat should be equal to number of rows of Beta and number of columns of SigmaHat should be equal to number of columns of Sxx.")

      }



      # Calculate Syx matrix
      Syx = t(t(Beta) *  diag(Sxx))

      # Criterion 1: Check if D is a numeric matrix with dimensions p * k, and entries are only 0 or 1.
      if (!is.numeric(D) || !is.matrix(D) || dim(D)[1] != p || dim(D)[2] != k || !all(D %in% c(0, 1))) {

        # Print an error message
        stop("D must be a numeric matrix with entries 0 or 1, where the number of rows should be equal to the number of Response variables (p) and the number of columns should be equal to the number of Instrumental Variables (k).")

      }

      # Criterion 2: Check if each trait has at least one valid IV
      # A valid IV for a trait i is an IV that affects only that trait.
      valid_IV_per_trait <- apply(D, 1, function(row) {
        any(row == 1 & colSums(D) == 1)
      })

      # If any trait does not have a valid IV, issue an error message
      if (!all(valid_IV_per_trait)) {

        stop("Each response variable does not have a valid IV and hence Sxx, Beta and SigmaHat can't be used as input.")

      }

      # Extract the first valid IV column index for each trait
      first_valid_IVs <- apply(D, 1, function(row) {
        valid_IV_columns <- which(row == 1 & colSums(D) == 1)
        valid_IV_columns[1]  # Return the first valid IV index
      })

      # Calculate Beta_New matrix by taking particular columns from Beta using Col_Ind
      Beta_New = matrix(Beta[, first_valid_IVs], nrow = p, ncol = p)

      # Initiate Estimate of A
      AEst = matrix(0, nrow = p, ncol = p)

      # Update AEst
      for (i in 1:p) {

        # Remove ith row and ith column of Beta_New
        Beta_ii = Beta_New[-i, -i]
        Beta_ii_det = det(Beta_ii)

        for (j in 1:p) {

          if(j != i){

            # Remove jth row and ith column of Beta_New
            beta_ji = Beta_New[-j, -i]

            AEst[i, j] = det(beta_ji) / Beta_ii_det * (-1)^(i + j + 1)

          }

        }

      }

      # Calculate (I-AEst)^(-1)
      Mult_Mat = solve(diag(p) - AEst)

      # Calculate diagonal entries of Syy
      Syy_Diag = SigmaHat[, 1] + 2 * Beta[, 1] * Syx[, 1] - Beta[, 1]^2 * Sxx[1, 1]

      # Calculate correlation matrix between X and between Y and X
      Cor_X = stats::cov2cor(Sxx)
      Cor_YX = t(t(Beta * (1 / sqrt(Syy_Diag))) * sqrt(diag(Sxx)))

      # Calculate determinant of Cor_X
      Cor_X_Det = det(Cor_X)

      # Calculate sum of squares of error while fitting Y[, i] on X
      Error_Sumsq = rep(0, p)

      for (i in 1:p) {

        # Create R
        R = cbind(c(1, Cor_YX[i,]), rbind(Cor_YX[i,], Cor_X))

        # Update Error_Sumsq
        Error_Sumsq[i] = Syy_Diag[i] * det(R) / Cor_X_Det


      }

      # calculate Sigma for the model
      Sigma = solve(Mult_Mat^2, Error_Sumsq)

      # Calculate Syy
      Syy = Beta %*% Sxx %*% t(Beta)  + Mult_Mat %*% (t(Mult_Mat) * Sigma)


    }


    # Criterion 1: Check if D is a numeric matrix with dimensions p * k, and entries are only 0 or 1.
    if (!is.numeric(D) || !is.matrix(D) || dim(D)[1] != p || dim(D)[2] != k || !all(D %in% c(0, 1))) {

      # Print an error message
      stop("D must be a numeric matrix with entries 0 or 1, where the number of rows should be equal to the number of Response variables (p) and the number of columns should be equal to the number of Instrumental Variables (k).")

    }

    # Criterion 2: Check if each trait has at least one valid IV
    # A valid IV for a trait i is an IV that affects only that trait.
    valid_IV_per_trait <- apply(D, 1, function(row) {
      any(row == 1 & colSums(D) == 1)
    })

    # If any trait does not have a valid IV, issue a warning
    if (!all(valid_IV_per_trait)) {

      warning("Each response variable does not have a valid IV and hence the model might be non-identifiable.")

    }


    # Check whether all the beta and inverse gamma parameters are positive or not
    if(!is.numeric(aRho) || aRho < 0 || !is.numeric(bRho) || bRho < 0 || !is.numeric(aPsi) ||aPsi < 0 || !is.numeric(bPsi) ||bPsi < 0 || !is.numeric(aSigma) || aSigma < 0 || !is.numeric(bSigma) || bSigma < 0){

      # Print an error message
      stop("All the beta and inverse gamma parameters should be positive.")

    }


    # Check whether the variance terms are positive or not
    if(!is.numeric(nu1) || nu1 <= 0 || !is.numeric(nu2) ||nu2 <= 0 || !is.numeric(PropVarA) || PropVarA <= 0 || !is.numeric(PropVarB) || PropVarB <= 0){

      # Print an error message
      stop("All the variance terms should be positive.")

    }


    # Check whether nIter is a postive integer
    if(!is.numeric(nIter) || nIter != round(nIter) || nIter <= 0){

      # Print an error message
      stop("Number of iterations should be a positive integer.")

    }

    # Check whether nBurnin is a non-negative integer and strictly less than nIter
    if(!is.numeric(nBurnin) || nBurnin != round(nBurnin) || nBurnin < 0 || nBurnin >= nIter){

      # Print an error message
      stop("Number of Burnin points should be a nonnegative integer strictly less than number of iterations.")

    }

    # Check whether Thin is a positive integer less than equal to (nIter - nBurnin)
    if(!is.numeric(Thin) || Thin != round(Thin) || Thin <= 0 || Thin > (nIter - nBurnin)){

      # Print an error message
      stop("Thin should be a positive integer less than or equal to the difference between number of iterations and number of burnin points.")

    }


    # Apply RGM for Spike and Slab prior and Threshold prior
    if ("Spike and Slab" %in% prior){

      # Run the algorithm for Spike and Slab prior
      Output = RGM_SpikeSlab2(Syy, Syx, Sxx, D, n, nIter = nIter, nBurnin = nBurnin, Thin = Thin,
                              a_rho = aRho, b_rho = bRho, nu_1 = nu1, a_psi =aPsi, b_psi =bPsi,
                              nu_2 =nu2, a_sigma = aSigma, b_sigma = bSigma, Prop_VarA = PropVarA, Prop_VarB = PropVarB)



      # Return outputs
      return(list(AEst = Output$A_Est, BEst = Output$B_Est, zAEst = Output$zA_Est, zBEst = Output$zB_Est,
                  GammaEst = Output$Gamma_Est, TauEst = Output$Tau_Est, RhoEst = Output$Rho_Est,
                  PhiEst = Output$Phi_Est, EtaEst = Output$Eta_Est, PsiEst = Output$Psi_Est,
                  SigmaEst = Output$Sigma_Est,
                  AccptA = Output$AccptA, AccptB = Output$AccptB,
                  LLPst = Output$LL_Pst, GammaPst = Output$Gamma_Pst))



    } else if ("Threshold" %in% prior){

      # Run the algorithm for Threshold prior
      Output = RGM_Threshold2(Syy, Syx, Sxx, D, n, nIter = nIter, nBurnin = nBurnin, Thin = Thin, nu_1 = nu1, nu_2 =nu2,
                              a_sigma = aSigma, b_sigma = bSigma, Prop_VarA = PropVarA, Prop_VarB = PropVarB)



      # Return outputs
      return(list(AEst = Output$A_Est, BEst = Output$B_Est, zAEst = Output$zA_Est, zBEst = Output$zB_Est,
                  A0Est = Output$A0_Est, B0Est = Output$B0_Est, GammaEst = Output$Gamma_Est, TauEst = Output$Tau_Est,
                  PhiEst = Output$Phi_Est, EtaEst = Output$Eta_Est, tAEst = Output$tA_Est, tBEst = Output$tB_Est,
                  SigmaEst = Output$Sigma_Est,
                  AccptA = Output$AccptA, AccptB = Output$AccptB, AccpttA = Output$Accpt_tA, AccpttB = Output$Accpt_tB,
                  LLPst = Output$LL_Pst, GammaPst = Output$Gamma_Pst))





    } else {

      # Print an error message
      stop("Please specify a prior among Threshold prior and Spike and Slab prior.")


    }


  } else {

    # Print an error message
    stop("Please give Y or Syy or X, Y or Syy, Syx, Sxx or Sxx, Beta and SigmaHat as input.")


  }

}
