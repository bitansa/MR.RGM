#' Fitting Multivariate Bidirectional Mendelian Randomization Models for Causal Gene Regulatory Network
#'
#' @description Using DNA or SNPs as the instrument and protein or gene expression data as the response, RGM may be used to design Multivariate Bidirectional Mendelian Randomization Models (Reciprocal Graphical Models) for Causal Gene Regulatory Network. Three input formats are supported by RGM: (i) The protein or gene expression data matrix Y and the DNA or SNP data matrix X; (ii) S_YY(t(Y)%*% Y / n), S_YX(t(Y)%*% X / n), and S_XX(t(X)%*% X / n) are summary level data; (iii) S_XX, Beta, and Sigma_Hat. If the user wishes to utilise (iii) as input, then compute these after centralising the data matrices.
#'              The graph structures between the response variables and between the response variables and the instruments will be returned by RGM as output. Additionally, the user must supply d, a vector input whose length is equal to that of the proteins or genes, and whose elements are individually positive integers representing the number of DNA or SNPs that impact a certain protein or gene, the total of which must match the number of DNA or SNPs. Additionally, the total number of observations, n, must also be given as input. The prior assumption on the model may either be "Threshold" or "Spike & Slab",
#'              with the former being the default. The user can also specify this assumption.
#'
#' @param X n * k matrix data input, each row denotes a particular observation and each column denotes a particular DNA or SNP. Default value is NULL.
#' @param Y n * p matrix data input, each row denotes a particular observation and each column denotes a particular protein or gene. Default value is NULL.
#' @param S_YY p * p matrix data input, where p denotes number of proteins or genes. It is obtained by t(Y) %*% Y / n.
#' @param S_YX p * k matrix data input, where p denotes number of proteins or genes and k denotes number of DNAs or SNPs. It is obtained by t(Y) %*% X / n.
#' @param S_XX k * k matrix data input, where k denotes number of DNAs or SNPs. It is obtained by t(X) %*% X / n.
#' @param Beta p * k matrix data input, where each row corresponds to a particular protein or gene and each column corresponds to a particular DNA or SNP. Each entry is the regression coefficient of the particular protein or gene on the particular DNA or SNP. If you want to use Beta as an input first centralize each column of Y and X and then calculate Beta, S_XX and Sigma_Hat.
#' @param Sigma_Hat p * k matrix data input, where each row corresponds to a particular protein or gene and each column corresponds to a particular DNA or SNP. Each entry is the mean square error for regressing the particular protein or gene on the particular DNA or SNP. If you want to use Sigma_Hat as an input first centralize each column of Y and X and then calculate Beta, S_XX and Sigma_Hat.
#' @param d Vector input of length p. Each element is a positive integer corresponds to number of DNA or SNP affecting a particular protein or gene, sum of which should be k.
#' @param n Positive integer input corresponding to number of datapoints.
#' @param nIter Positive integer input corresponding to number of MCMC sampling. Default value is 10,000.
#' @param nBurnin Non-negative integer input corresponding to number of samples to be discarded. nBurnin should be less than nIter. Default value is 2000.
#' @param Thin Positive integer input corresponding to thinning of posterior samples. Thin should be less than or equal to (nIter - nBurnin). Default value is 1.
#' @param prior Prior assumption on the graph structure. It can either be "Threshold" or "Spike and Slab". Default is "Threshold".
#' @param a_tau Positive scalar input. It corresponds to the first parameter of an Inverse Gamma distribution. Default value is 0.01.
#' @param b_tau Positive scalar input. It corresponds to the second parameter of an Inverse Gamma distribution. Default value is 0.01.
#' @param a_rho Positive scalar input. It corresponds to the first parameter of a Beta distribution. Default value is 0.5.
#' @param b_rho Positive scalar input. It corresponds to the second parameter of a Beta distribution. Default value is 0.5.
#' @param nu_1 Positive scalar input. It corresponds to the variance of slab part of spike and slab distribution of A. Default value is 0.0001.
#' @param a_eta Positive scalar input. It corresponds to the first parameter of an Inverse Gamma distribution. Default value is 0.01.
#' @param b_eta Positive scalar input. It corresponds to the second parameter of an Inverse Gamma distribution. Default value is 0.01.
#' @param a_psi Positive scalar input. It corresponds to the first parameter of a Beta distribution. Default value is 0.5.
#' @param b_psi Positive scalar input. It corresponds to the second parameter of a Beta distribution. Default value is 0.5.
#' @param nu_2 Positive scalar input. It corresponds to the variance of slab part of spike and slab distribution of B. Default value is 0.0001.
#' @param a_sigma Positive scalar input. It corresponds to the first parameter of an Inverse Gamma distribution corresponding to the variance of the model. Default value is 0.01.
#' @param b_sigma positive scalar input. It corresponds to the second parameter of an Inverse Gamma distribution corresponding to the variance of the model. Default value is 0.01.
#' @param Prop_VarA positive scalar input. It corresponds to the variance of the normal distribution for proposing A matrix terms. Default value is 0.01.
#' @param Prop_VarB positive scalar input. It corresponds to the variance of the normal distribution for proposing B matrix terms. Default value is 0.01.
#'
#' @return
#'
#' \item{A_Est}{p * p matrix output of causal effect or strength between the proteins or the genes.}
#' \item{B_Est}{p * k matrix output of causal effect or strength between the proteins or the genes and the DNAs or SNPs.}
#' \item{zA_Est}{p * p binary adjacency matrix representing the graph structure between the proteins or the genes.}
#' \item{zB_Est}{p * k binary adjacency matrix representing the graph structure between the proteins or the genes and the DNAs or SNPs.}
#' \item{A0_Est}{p * p matrix output of protein-protein interactions before thresholding for Threshold prior.}
#' \item{B0_Est}{p * k matrix output of protein-DNA interactions before thresholding for Threshold prior.}
#' \item{Gamma_Est}{p * p matrix output of probabilities of protein-protein interactions.}
#' \item{Tau_Est}{p * p matrix output of variances of protein-protein interactions.}
#' \item{Phi_Est}{p * k matrix output of probabilities of protein-DNA interactions.}
#' \item{Eta_Est}{p * k matrix output of variances of protein-DNA interactions.}
#' \item{tA_Est}{Scalar output of thresholding value of protein-protein interactions for Threshold prior.}
#' \item{tB_Est}{Scalar output of thresholding value of protein-DNA interactions for Threshold prior.}
#' \item{Sigma_Est}{Vector output of length p corresponding to variances of each protein.}
#' \item{AccptA}{Percentage acceptance of entries of A matrix i.e. protein-protein interactions.}
#' \item{AccptB}{Percentage acceptance of entries of B matrix i.e. protein-DNA interactions.}
#' \item{Accpt_tA}{Percentage acceptance of the thresholding value for protein-protein interactions for Threshold prior.}
#' \item{Accpt_tB}{Percentage acceptance of the thresholding value for protein-DNA interactions for Threshold prior.}
#' \item{LL_Pst}{Vector ouput of posterior log-likelihoods of the model.}
#' \item{Rho_Est}{p * p matrix output of bernoulli success probabilities of protein-protein interactions for Spike and Slab prior.}
#' \item{Psi_Est}{p * k matrix output of bernoulli success probabilities of protein-DNA interactions for Spike and Slab prior.}
#'
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
#' # -----------------------------------------------------------------
#' # Example 1:
#' # Run RGM based on individual level data with Threshold prior
#'
#' # Data Generation
#' set.seed(9154)
#'
#' # Number of datapoints
#' n = 10000
#'
#' # Number of Genes and number of DNAs
#' p = 5
#' k = 6
#'
#' # Initialize gene-gene interaction matrix
#' A = matrix(sample(c(-0.1, 0.1), p^2, replace = TRUE), p, p)
#'
#' # Diagonal entries of A matrix will always be 0
#' diag(A) = 0
#'
#' # Make the network sparse
#' A[sample(which(A!=0), length(which(A!=0))/2)] = 0
#'
#' # Initialize gene-DNA interaction matrix
#' B = matrix(0, p, k)
#'
#' # Create d vector
#' d = c(2, 1, 1, 1, 1)
#'
#'
#' # Initialize m
#' m = 1
#'
#' # Calculate B matrix based on d vector
#' for (i in 1:p) {
#'
#'  # Update ith row of B
#'  B[i, m:(m + d[i] - 1)] = 1
#'
#'  # Update m
#'  m = m + d[i]
#'
#' }
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
#' Y = matrix(0, nrow = n, ncol = p)
#'
#' for (i in 1:n) {
#'
#'     Y[i, ] = MASS::mvrnorm(n = 1, Mult_Mat %*% B %*% X[i, ], Variance)
#'
#' }
#'
#' # Print true Protein-Protein interaction matrix and Protein-DNA interaction matrix
#' A
#' B
#'
#' # Apply RGM on individual level data for Threshold Prior
#' Output = RGM(X = X, Y = Y, d = c(2, 1, 1, 1, 1), n = 10000, prior = "Threshold")
#'
#' # Get the graph structure between the proteins
#' Output$zA_Est
#'
#' # Get Protein-Protein interactions
#' Output$A_Est
#'
#' # Get the graph structure between the proteins and the DNAs
#' Output$zB_Est
#'
#' # Get Protein-DNA interactions
#' Output$B_Est
#'
#' # Plot posterior log-likelihood
#' plot(Output$LL_Pst, type = 'l', xlab = "Number of Iterations", ylab = "Log-likelihood")
#'
#'
#'
#'
#' # -----------------------------------------------------------------
#' # Example 2:
#' # Run RGM based on summary level data with Spike and Slab prior
#'
#' # Data Generation
#' set.seed(9154)
#'
#' # Number of datapoints
#' n = 10000
#'
#' # Number of Genes and number of DNAs
#' p = 5
#' k = 6
#'
#' # Initialize gene-gene interaction matrix
#' A = matrix(sample(c(-0.1, 0.1), p^2, replace = TRUE), p, p)
#'
#' # Diagonal entries of A matrix will always be 0
#' diag(A) = 0
#'
#' # Make the network sparse
#' A[sample(which(A!=0), length(which(A!=0))/2)] = 0
#'
#' # Initialize gene-DNA interaction matrix
#' B = matrix(0, p, k)
#'
#' # Create d vector
#' d = c(2, 1, 1, 1, 1)
#'
#'
#' # Initialize m
#' m = 1
#'
#' # Calculate B matrix based on d vector
#' for (i in 1:p) {
#'
#'  # Update ith row of B
#'  B[i, m:(m + d[i] - 1)] = 1
#'
#'  # Update m
#'  m = m + d[i]
#'
#' }
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
#' Y = matrix(0, nrow = n, ncol = p)
#'
#' for (i in 1:n) {
#'
#'     Y[i, ] = MASS::mvrnorm(n = 1, Mult_Mat %*% B %*% X[i, ], Variance)
#'
#' }
#'
#'
#' # Calculate summary level data
#' S_YY = t(Y) %*% Y / n
#' S_YX = t(Y) %*% X / n
#' S_XX = t(X) %*% X / n
#'
#'
#' # Print true Protein-Protein interaction matrix and Protein-DNA interaction matrix
#' A
#' B
#'
#'
#' # Apply RGM on summary level data for Spike and Slab Prior
#' Output = RGM(S_YY = S_YY, S_YX = S_YX, S_XX = S_XX,
#'           d = c(2, 1, 1, 1, 1), n = 10000, prior = "Spike and Slab")
#'
#' # Get the graph structure between the proteins
#' Output$zA_Est
#'
#' # Get Protein-Protein interactions
#' Output$A_Est
#'
#' # Get the graph structure between the proteins and the DNAs
#' Output$zB_Est
#'
#' # Get Protein-DNA interactions
#' Output$B_Est
#'
#' # Plot posterior log-likelihood
#' plot(Output$LL_Pst, type = 'l', xlab = "Number of Iterations", ylab = "Log-likelihood")
#'
#'
#'
#'
#' # -----------------------------------------------------------------
#' # Example 3:
#' # Run RGM based on Beta and Sigma_Hat with Threshold prior
#'
#' # Data Generation
#' set.seed(9154)
#'
#' # Number of datapoints
#' n = 10000
#'
#' # Number of Genes and number of DNAs
#' p = 5
#' k = 6
#'
#' # Initialize gene-gene interaction matrix
#' A = matrix(sample(c(-0.1, 0.1), p^2, replace = TRUE), p, p)
#'
#' # Diagonal entries of A matrix will always be 0
#' diag(A) = 0
#'
#' # Make the network sparse
#' A[sample(which(A!=0), length(which(A!=0))/2)] = 0
#'
#' # Initialize gene-DNA interaction matrix
#' B = matrix(0, p, k)
#'
#' # Create d vector
#' d = c(2, 1, 1, 1, 1)
#'
#'
#' # Initialize m
#' m = 1
#'
#' # Calculate B matrix based on d vector
#' for (i in 1:p) {
#'
#'  # Update ith row of B
#'  B[i, m:(m + d[i] - 1)] = 1
#'
#'  # Update m
#'  m = m + d[i]
#'
#' }
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
#' Y = matrix(0, nrow = n, ncol = p)
#'
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
#' # Calculate S_XX
#' S_XX = t(X) %*% X / n
#'
#' # Generate Beta matrix and Sigma_Hat
#' Beta = matrix(0, nrow = p, ncol = k)
#' Sigma_Hat = matrix(0, nrow = p, ncol = k)
#'
#' for (i in 1:p) {
#'
#'    for (j in 1:k) {
#'
#'        fit = lm(Y[, i] ~ X[, j])
#'
#'        Beta[i, j] =  fit$coefficients[2]
#'
#'        Sigma_Hat[i, j] = sum(fit$residuals^2) / n
#'
#'        }
#'
#'    }
#'
#'
#' # Print true Protein-Protein interaction matrix and Protein-DNA interaction matrix
#' A
#' B
#'
#'
#' # Apply RGM based on S_XX, Beta adn Sigma_Hat for Threshold Prior
#' Output = RGM(S_XX = S_XX, Beta = Beta, Sigma_Hat = Sigma_Hat,
#'           d = c(2, 1, 1, 1, 1), n = 10000, prior = "Threshold")
#'
#' # Get the graph structure between the proteins
#' Output$zA_Est
#'
#' # Get Protein-Protein interactions
#' Output$A_Est
#'
#' # Get the graph structure between the proteins and the DNAs
#' Output$zB_Est
#'
#' # Get Protein-DNA interactions
#' Output$B_Est
#'
#' # Plot posterior log-likelihood
#' plot(Output$LL_Pst, type = 'l', xlab = "Number of Iterations", ylab = "Log-likelihood")
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
RGM = function(X = NULL, Y = NULL, S_YY = NULL, S_YX = NULL, S_XX = NULL, Beta = NULL, Sigma_Hat = NULL, d, n, nIter = 10000, nBurnin = 2000, Thin = 1, prior = c("Threshold", "Spike and Slab"), a_tau = 0.01, b_tau = 0.01, a_rho = 0.5, b_rho = 0.5, nu_1 = 0.0001, a_eta = 0.01, b_eta = 0.01, a_psi = 0.5, b_psi = 0.5, nu_2 = 0.0001, a_sigma = 0.01, b_sigma = 0.01, Prop_VarA = 0.01, Prop_VarB = 0.01){

  # Check whether n is a positive integer
  if(!is.numeric(n) || n != round(n) || n <= 0){

    # Print an error message
    stop("Number of datapoints should be a positive integer.")

  }

  # Check whether input X and Y are given
  if(!is.null(X) && !is.null(Y)){

    # Check whether X and Y are both numeric matrices
    if(!is.numeric(X) || !is.numeric(Y) || !is.matrix(X) || !is.matrix(Y)){

      # Print an error message
      stop("X and Y should be numeric matrices.")

    }

    # Check whether number of rows of X and Y are n
    if(nrow(X) != n || nrow(Y) != n){

      # Print an error message
      stop("Number of rows of X and Y should be n.")

    }

    # Calculate number of nodes from Y matrix
    p = ncol(Y)

    # Calculate number of covariates from X matrix
    k = ncol(X)

    # Calculate S_YY, S_YX, S_XX
    S_YY = t(Y) %*% Y / n
    S_YX = t(Y) %*% X / n
    S_XX = t(X) %*% X / n

  } else if(!is.null(S_YY) && !is.null(S_YX) && !is.null(S_XX)){

    # Check whether S_YY ,S_YX and S_XX are numeric matrices
    if(!is.numeric(S_YY) || !is.numeric(S_YX) || !is.numeric(S_XX) || !is.matrix(S_YY) || !is.matrix(S_YX) || !is.matrix(S_XX)){

      # Print an error message
      stop("S_YY, S_YX and S_XX should be numeric matrices.")

    }

    # Calculate number of nodes from S_YY matrix
    p = ncol(S_YY)

    # Calculate number of covariates from S_XX matrix
    k = ncol(S_XX)

    # Check whether number of rows of S_YY is equal to p
    if(nrow(S_YY) != p){

      # Print an error message
      stop("S_YY should be a diagonal matrix.")

    }

    # Check whether number of rows of S_XX is equal to k
    if(nrow(S_XX) != k){

      # Print an error message
      stop("S_XX should be a diagonal matrix.")

    }

    # Check whether number of rows of S_YX is equal to p and number of columns of S_YX is equal to k
    if(nrow(S_YX) != p || ncol(S_YX) != k){

      # Print an error message
      stop("Number of rows of S_YX should be equal to number of rows of S_YY and number of columns of S_YX should be equal to number of columns of S_XX.")

    }



  } else if(!is.null(S_XX) && !is.null(Beta) && !is.null(Sigma_Hat)){

    # Check whether S_XX, Beta and Sigma_Hat are numeric matrices
    if(!is.numeric(S_XX) || !is.numeric(Beta) || !is.numeric(Sigma_Hat) || !is.matrix(S_XX) || !is.matrix(Beta) || !is.matrix(Sigma_Hat)){

      # Print an error message
      stop("S_XX, Beta and SIgma_Hat should be numeric matrices.")

    }

    # Calculate number of nodes from Beta matrix
    p = nrow(Beta)

    # Calculate number of covariates from S_XX matrix
    k = ncol(S_XX)

    # Check whether number of rows of S_XX is equal to k
    if(nrow(S_XX) != k){

      # Print an error message
      stop("S_XX should be a diagonal matrix.")

    }

    # Check whether number of columns of Beta is equal to k
    if(ncol(Beta) != k){

      # Print an error message
      stop("Number of cloumns of Beta should be equal to number of columns of S_XX.")

    }

    # Check whether number of rows of Sigma_Hat is equal to p and number of columns of Sigma_Hat is equal to k
    if(nrow(Sigma_Hat) != p || ncol(Sigma_Hat) != k){

      # Print an error message
      stop("Number of rows of Sigma_Hat should be equal to number of rows of Beta and number of columns of Sigma_Hat should be equal to number of columns of S_XX.")

    }



    # Calculate S_YX matrix
    S_YX = t(t(Beta) *  diag(S_XX))

    # Check whether d is a vector of non-negative integers of length p and sum od entries of d is equal to k
    if(!is.numeric(d) || sum(d != round(d)) != 0 || sum(d <= 0) != 0 || length(d) != p || sum(d) != k){

      # Print an error message
      stop("d should be a vector of positive integers of length equal to number of nodes and sum of entries should be equal to number of covariates.")

    }

    # Store indices to extract particular columns from beta matrix
    Col_Ind = c(0, cumsum(d))
    Col_Ind = Col_Ind[-length(Col_Ind)] + 1

    # Calculate Beta_New matrix by taking particular columns from Beta using Col_Ind
    Beta_New = matrix(Beta[, Col_Ind], nrow = p, ncol = p)

    # Initiate Estimate of A
    A_Est = matrix(0, nrow = p, ncol = p)

    # Update A_Est
    for (i in 1:p) {

      # Remove ith row and ith column of Beta_New
      Beta_ii = Beta_New[-i, -i]
      Beta_ii_det = det(Beta_ii)

      for (j in 1:p) {

        if(j != i){

          # Remove jth row and ith column of Beta_New
          beta_ji = Beta_New[-j, -i]

          A_Est[i, j] = det(beta_ji) / Beta_ii_det * (-1)^(i + j + 1)

        }

      }

    }

    # Calculate (I-A_Est)^(-1)
    Mult_Mat = solve(diag(p) - A_Est)

    # Calculate diagonal entries of S_YY
    S_YY_Diag = Sigma_Hat[, 1] + 2 * Beta[, 1] * S_YX[, 1] - Beta[, 1]^2 * S_XX[1, 1]

    # Calculate correlation matrix between X and between Y and X
    Cor_X = stats::cov2cor(S_XX)
    Cor_YX = t(t(Beta * (1 / sqrt(S_YY_Diag))) * sqrt(diag(S_XX)))

    # Calculate determinant of Cor_X
    Cor_X_Det = det(Cor_X)

    # Calculate sum of squares of error while fitting Y[, i] on X
    Error_Sumsq = rep(0, p)

    for (i in 1:p) {

      # Create R
      R = cbind(c(1, Cor_YX[i,]), rbind(Cor_YX[i,], Cor_X))

      # Update Error_Sumsq
      Error_Sumsq[i] = S_YY_Diag[i] * det(R) / Cor_X_Det


    }

    # calculate Sigma for the model
    Sigma = solve(Mult_Mat^2, Error_Sumsq)

    # Calculate S_YY
    S_YY = Beta %*% S_XX %*% t(Beta)  + Mult_Mat %*% (t(Mult_Mat) * Sigma)


  } else {

    # Print an error message
    stop("Please give X, Y or S_YY, S_YX, S_XX or S_XX, Beta and Sigma_Hat as input.")


  }



  # Check whether d is a vector of positive integers of length p and sum of entries of d is equal to k
  if(!is.numeric(d) || sum(d != round(d)) != 0 || sum(d <= 0) != 0 || length(d) != p || sum(d) != k){

    # Print an error message
    stop("d should be a vector of positive integers of length equal to number of nodes and sum of entries should be equal to number of covariates.")

  }

  # Initialize D matrix
  D = matrix(0, p, k)

  # Initialize m
  m = 1

  # Calculate D matrix based on d vector
  for (i in 1:p) {

    # Update ith row of D
    D[i, m:(m + d[i] - 1)] = 1

    # Update m
    m = m + d[i]

  }


  # Check whether the inverse gamma parameters are positive or not
  if(!is.numeric(a_tau) || a_tau < 0 || !is.numeric(b_tau) || b_tau < 0 || !is.numeric(a_rho) || a_rho < 0 || !is.numeric(b_rho) || b_rho < 0 || !is.numeric(a_eta) || a_eta < 0 || !is.numeric(b_eta) || b_eta < 0 || !is.numeric(a_psi) || a_psi < 0 || !is.numeric(b_psi) || b_psi < 0 || !is.numeric(a_sigma) || a_sigma < 0 || !is.numeric(b_sigma) || b_sigma < 0){

    # Print an error message
    stop("All the inverse gamma parameters should be positive")

  }


  # Check whether the variance terms are positive or not
  if(!is.numeric(nu_1) || nu_1 <= 0 || !is.numeric(nu_2) || nu_2 <= 0 || !is.numeric(Prop_VarA) || Prop_VarA <= 0 || !is.numeric(Prop_VarB) || Prop_VarB <= 0){

    # Print an error message
    stop("All the variance terms should be positive")

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


  # Apply RGM for Threshold prior and Spike and Slab prior
  if ("Threshold" %in% prior){

    # Run the algorithm for Threshold prior
    Output = RGM_Threshold(S_YY, S_YX, S_XX, D, n, nIter = nIter, nBurnin = nBurnin, Thin = Thin, nu_1 = nu_1, nu_2 = nu_2,
                           a_sigma = a_sigma, b_sigma = b_sigma, Prop_VarA = Prop_VarA, Prop_VarB = Prop_VarB)


    # Return outputs
    return(list(A_Est = Output$A_Est, B_Est = Output$B_Est, zA_Est = Output$zA_Est, zB_Est = Output$zB_Est,
                A0_Est = Output$A0_Est, B0_Est = Output$B0_Est, Gamma_Est = Output$Gamma_Est, Tau_Est = Output$Tau_Est,
                Phi_Est = Output$Phi_Est, Eta_Est = Output$Eta_Est, tA_Est = Output$tA_Est, tB_Est = Output$tB_Est,
                Sigma_Est = Output$Sigma_Est,
                AccptA = Output$AccptA, AccptB = Output$AccptB, Accpt_tA = Output$Accpt_tA, Accpt_tB = Output$Accpt_tB,
                LL_Pst = Output$LL_Pst))





  } else if ("Spike and Slab" %in% prior){

    # Run the algorithm for Threshold prior
    Output = RGM_SpikeSlab(S_YY, S_YX, S_XX, D, n, nIter = nIter, nBurnin = nBurnin, Thin = Thin, a_tau = a_tau, b_tau = b_tau,
                           a_rho = a_rho, b_rho = b_rho, nu_1 = nu_1, a_eta = a_eta, b_eta = b_eta, a_psi = a_psi, b_psi = b_psi,
                           nu_2 = nu_2, a_sigma = a_sigma, b_sigma = b_sigma, Prop_VarA = Prop_VarA, Prop_VarB = Prop_VarB)



    # Return outputs
    return(list(A_Est = Output$A_Est, B_Est = Output$B_Est, zA_Est = Output$zA_Est, zB_Est = Output$zB_Est,
                Gamma_Est = Output$Gamma_Est, Tau_Est = Output$Tau_Est, Rho_Est = Output$Rho_Est,
                Phi_Est = Output$Phi_Est, Eta_Est = Output$Eta_Est, Psi_Est = Output$Psi_Est,
                Sigma_Est = Output$Sigma_Est,
                AccptA = Output$AccptA, AccptB = Output$AccptB,
                LL_Pst = Output$LL_Pst))



  } else{

    # Print an error message
    stop("Please specify a prior among Threshold prior and Spike and Slab prior.")


  }




}



