#' Rho generating function
#'
#' @param Gamma matrix input
#' @param p positive scalar input
#' @param a_rho positive scalar input
#' @param b_rho positive scalar input
#'
#' @return scalar value generated from beta distribution, rbeta(1, sum(Gamma) + a_rho, p * (p - 1) - sum(Gamma) + b_rho)
#'
#' @examples
Generate_Rho = function(Gamma, p, a_rho, b_rho){

  # Calculate sum of all entries of Gamma matrix
  Gamma_sum = sum(Gamma)

  # Generate Rho from Beta distribution
  Rho = stats::rbeta(1, Gamma_sum + a_rho, p * (p - 1) - Gamma_sum + b_rho)

  # Return Rho
  return(Rho)

}





#' Psi generating function
#'
#' @param Phi matrix input
#' @param d non-negative scalar input
#' @param a_psi positive scalar input
#' @param b_psi positive scalar input
#'
#' @return scalar value generated from beta distribution, rbeta(1, sum(Phi) + a_psi, d - sum(Phi) + b_psi)
#'
#' @examples
Generate_Psi = function(Phi, d, a_psi, b_psi){

  # Calculate sum of all entries of Phi matrix
  Phi_sum = sum(Phi)

  # Generate Psi from Beta distribution
  Psi = stats::rbeta(1, Phi_sum + a_psi, d - Phi_sum + b_psi)

  # Return Rho
  return(Psi)

}







#' Eta generating function
#'
#' @param b scalar input
#' @param phi scalar input either 0 or 1
#' @param a_eta positive scalar input
#' @param b_eta positive scalar input
#' @param nu_2 positve scalar input
#'
#' @return scalar value generated from inverse gamma distribution, phi / rgamma(1, a_eta + 1/2, b^2/2 + b_eta) + (1 - phi) / rgamma(1, a_eta + 1/2, b^2/(2 * nu_2) + b_eta)
#'
#' @examples
Generate_Eta = function(b, phi, a_eta, b_eta, nu_2){

  # Generate eta based on Inverse Gamma distribution
  eta = phi / stats::rgamma(1, a_eta + 1/2, b^2/2 + b_eta) + (1 - phi) / stats::rgamma(1, a_eta + 1/2, b^2/(2 * nu_2) + b_eta)

  # Return eta
  return(eta)

}






#' Tau generating function
#'
#' @param a scalar input
#' @param gamma scalar input either 0 or 1
#' @param a_tau positive scalar input
#' @param b_tau positive scalar input
#' @param nu_1 positive scalar input
#'
#' @return scalar value generated from inverse gamma distribution, gamma / rgamma(1, a_tau + 1/2, a^2/2 + b_tau) + (1 - gamma) / rgamma(1, a_tau + 1/2, a^2/(2 * nu_1) + b_tau)
#'
#' @examples
Generate_Tau = function(a, gamma, a_tau, b_tau, nu_1){

  # Generate tau based on Inverse Gamma distribution
  tau = gamma / stats::rgamma(1, a_tau + 1/2, a^2/2 + b_tau) + (1 - gamma) / stats::rgamma(1, a_tau + 1/2, a^2/(2 * nu_1) + b_tau)

  # Return tau
  return(tau)

}






#' Phi generating function
#'
#' @param b scalar input
#' @param psi scalar input between 0 and 1
#' @param eta positive scalar input
#' @param nu_2 positive scalar input
#'
#' @return scalar value generated from a bernoulli distribution i.e. 0 or 1
#'
#' @examples
Generate_Phi = function(b, psi, eta, nu_2){

  # Calculate C1 and C2
  C1 = psi / sqrt(eta) * exp(- b^2 / (2 * eta))

  C2 = (1 - psi) / sqrt(nu_2 * eta) * exp(- b^2 / (2 * nu_2 * eta))

  # Generate phi based on bernoulli distribution
  phi = stats::rbinom(1, 1, C1/(C1+C2))

  # Return phi
  return(phi)

}





#' Gamma generating function
#'
#' @param a scalar input
#' @param rho scalar input between 0 and 1
#' @param tau  positive scalar input
#' @param nu_1 positive scalar input
#'
#' @return scalar value generated from a bernoulli distribution i.e. 0 or 1
#'
#' @examples
Generate_Gamma = function(a, rho, tau, nu_1){

  # Calculate D1 and D2
  D1 = rho / sqrt(tau) * exp(- a^2 / (2 * tau))

  D2 = (1 - rho) / sqrt(nu_1 * tau) * exp(- a^2 / (2 * nu_1 * tau))

  # Generate gamma based on bernoulli distribution
  gamma = stats::rbinom(1, 1, D1/(D1+D2))

  # Return gamma
  return(gamma)

}









#' Sigma generating function
#'
#' @param n positive integer input
#' @param z_sum non-negative scalar input
#' @param a_sigma positive scalar input
#' @param b_sigma positive scalar input
#'
#' @return scalar value generated from inverse gamma distribution, 1 / rgamma(1, n/2 + a_sigma, z_sum/2 + b_sigma)
#'
#' @examples
Generate_Sigma = function(n, z_sum, a_sigma, b_sigma){

  # Generate sigma based on Inverse gamma distribution
  sigma = 1 / stats::rgamma(1, n/2 + a_sigma, z_sum/2 + b_sigma)

  # Return sigma
  return(sigma)

}









#' B matrix entries generating function
#'
#' @param xz scalar input
#' @param b_vec vector input
#' @param b scalar input
#' @param x_mat matrix input
#' @param x_vec vector input
#' @param sigma positive scalar input
#' @param eta positive scalar input
#' @param phi scalar input 0 or 1
#' @param nu2 positive scalar input
#'
#' @return scalar value generated from normal distribution
#'
#' @examples
Generate_B = function(xz, b_vec, b, x_mat, x_vec, sigma, eta, phi, nu2){

  # Calculate numerator of mean
  mean = xz - sum(tcrossprod(b_vec, x_mat) * x_vec) / sigma + sum(x_vec^2) * b / sigma

  # Calculate variance
  variance = phi / (sum(x_vec) / sigma + 1 / eta) + (1 - phi) / (sum(x_vec) / sigma + 1 / (eta * nu2))

  # Generate b
  b = stats::rnorm(1, mean = mean * variance, sqrt(variance))

  # Return beta
  return(b)

}









#' Target function for a particular A entry and corresponding gamma
#'
#' @param X n * K matrix input
#' @param Y n * p matrix input
#' @param A p * p matrix input
#' @param a scalar input
#' @param N positive scalar input
#' @param Sigma_Inv p * p matrix input
#' @param p positive integer input
#' @param B p * K matrix input
#' @param gamma scalar input, 0 or 1
#' @param tau positive scalar input
#' @param rho scalar input in between 0 and 1
#' @param nu_1 positive scalar input
#'
#' @return scalar target value
#'
#' @examples
Target_Agamma = function(X, Y, A, a, N, Sigma_Inv, p, B, gamma, tau, rho, nu_1){

  # Calculate (I_p - A) and (I_p - A)^{-1} * B
  Mult_Mat = diag(p) - A
  Inv_MatB = solve(Mult_Mat) %*% B


  # Calculate mean matrix and variance matrix inverse
  Mean_mat = tcrossprod(Inv_MatB, X)
  Var_mat_inv = crossprod(Mult_Mat, Sigma_Inv %*% Mult_Mat)

  # Initiate Sum
  Sum = 0

  for (i in 1:N) {

    # Calculate difference
    Diff = Y[i, ] - t(Mean_mat[, i])

    # Calculate Sum
    Sum = Sum + tcrossprod(Diff %*% Var_mat_inv, Diff)

  }

  # Calculate target value
  Target = N/2 * determinant(Var_mat_inv, logarithm = TRUE)$modulus -1/2 * Sum + gamma * (-a^2/(2 * tau) - log(sqrt(tau))) + (1 - gamma) * (-a^2/(2 * nu_1 * tau) - log(sqrt(nu_1 * tau))) + gamma * log(rho) + (1 - gamma) * log(1 - rho)


  # Return Target
  return(Target)

}







#' A matrix entries and gamma generating function
#'
#'
#' @inheritParams Target_Agamma
#' @param i integer input in between 1 and p
#' @param j integer input in between 1 and p
#' @param prop_var1 postive scalar input
#'
#' @return scalar a value and scalar gamma value, 0 or 1
#'
#' @examples
Generate_Agamma = function(X, Y, A, i, j, Sigma_Inv, N, p, B, gamma, tau, rho, nu_1, prop_var1){

  # Value to update
  a = A[i, j]

  # Proposed value
  a_new = stats::rnorm(1, a, prop_var1)

  # New A matrix with proposed a value
  A_new = A
  A_new[i, j] = a_new

  # Calculate r
  r = Target_Agamma(X, Y, A_new, a_new, N, Sigma_Inv, p, B, 1-gamma, tau, rho, nu_1) - Target_Agamma(X, Y, A, a, N, Sigma_Inv, p, B, gamma, tau, rho, nu_1)



  # Generate uniform u
  u = stats::runif(1, 0, 1)

  if(!is.na(r)){

    # Check whether r is big or not
    # min(1, r) >= u
    if(min(0, r) >= log(u)){

      a = a_new

      gamma = 1 - gamma

    }

  }

  return(list(a = a, gamma = gamma))

}







#' Target function for a particular A entry
#'
#' @inheritParams Target_Agamma
#'
#' @return scalar target value
#'
#' @examples
Target_A = function(X, Y, A, a, N, Sigma_Inv, p, B, gamma, tau, nu_1){

  # Calculate (I_p - A)^{-1} and (I_p - A)^{-1} * B
  Mult_Mat = diag(p) - A

  Inv_Mat = solve(Mult_Mat)
  Inv_MatB = Inv_Mat %*% B

  # Calculate mean matrix and variance matrix inverse
  Mean_mat = tcrossprod(Inv_MatB, X)
  Var_mat_inv = crossprod(Mult_Mat, Sigma_Inv %*% Mult_Mat)

  # Initiate Sum
  Sum = 0

  for (i in 1:N) {

    # Calculate difference
    Diff = Y[i, ] - t(Mean_mat[, i])

    # Calculate Sum
    Sum = Sum + tcrossprod(Diff %*% Var_mat_inv, Diff)

  }

  # Calculate target value
  Target = (det(Var_mat_inv)^(N/2)) * exp(-1/2 * Sum) * (gamma * exp(-a^2/(2 * tau)) / sqrt(tau) + (1 - gamma) * exp(- a^2/(2 * nu_1 * tau)) / sqrt(nu_1 * tau))


  # Return Target
  return(Target)

}








#' A matrix entries generating function
#'
#' @inheritParams Generate_Agamma
#'
#' @return scalar a value
#'
#' @examples
Generate_A = function(X, Y, A, i, j, Sigma_Inv, N, p, B, gamma, tau, nu_1, prop_var1){

  # Value to update
  a = A[i, j]

  # Proposed value
  a_new = stats::rnorm(1, a, prop_var1)

  # New A matrix with proposed a value
  A_new = A
  A_new[i, j] = a_new

  # Calculate r
  r = Target_A(X, Y, A_new, a_new, Sigma_Inv, N, p, B, gamma, tau, nu_1)  / Target_A(X, Y, A, a, Sigma_Inv, N, p, B, gamma, tau, nu_1)

  # Generate uniform u
  u = stats::runif(1, 0, 1)

  if(!is.na(r)){

    # Check whether r is big or not
    if(min(1, r) >= u){

      a = a_new

    }

  }

  return(a)

}








#' Target function for a particular A entry and corresponding gamma
#'
#' @inheritParams Generate_Agamma
#' @inheritParams Generate_B
#' @param Mult_Inv_Y p * n matrix input
#' @param psi scalar input in between 0 and 1
#' @param nu_2 positive scalar input
#'
#' @return scalar target value
#'
#' @examples
Target_Bphi = function(X, Y, B, Sigma_Inv, Mult_Inv_Y, b, phi, eta, psi, nu_2){

  # Calculate Z vector
  Z = Mult_Inv_Y - tcrossprod(B, X)


  # Calculate Sum
  Sum = sum(rowSums(Z^2) * Sigma_Inv)

  # Calculate target value
  Target =  -1/2 * Sum + phi * (-b^2/(2 * eta)) + (1 - phi) * (-b^2/(2 * nu_2 * eta)) + phi * log(psi) + (1 - phi) * log(1 - psi)


  # Return Target
  return(Target)

}






#' B matrix entries and phi generating function
#'
#' @inheritParams Target_Bphi
#'
#' @param i positive integer input in between 1 and p
#' @param j positive integer input in between 1 and K
#' @param prop_var2 positive scalar input
#'
#' @return scalar b value and scalar phi value, 0 or 1
#'
#' @examples
Generate_Bphi = function(X, Y, B, i, j, Sigma_Inv, Mult_Inv_Y, phi, eta, psi, nu_2, prop_var2){

  # Value to update
  b = B[i, j]

  # Proposed value
  b_new = stats::rnorm(1, b, prop_var2)

  # New B matrix with proposed b value
  B_new = B
  B_new[i, j] = b_new

  # Calculate r
  r = Target_Bphi(X, Y, B_new, Sigma_Inv, Mult_Inv_Y, b_new, 1 - phi, eta, psi, nu_2) - Target_Bphi(X, Y, B, Sigma_Inv, Mult_Inv_Y, b, phi, eta, psi, nu_2)



  # Generate uniform u
  u = stats::runif(1, 0, 1)

  if(!is.na(r)){

    # Check whether r is big or not
    # min(1, r) >= u
    if(min(0, r) >= log(u)){

      b = b_new

      phi = 1 - phi

    }

  }

  return(list(b = b, phi = phi))

}
























