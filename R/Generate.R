#' Rho generating function
#'
#' @param Gamma matrix input
#' @param p positive scalar input
#' @param a_rho positive scalar input, first parameter of beta distribution
#' @param b_rho positive scalar input, second parameter of beta distribution
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
#' @param a_psi positive scalar input, first parameter of beta distribution
#' @param b_psi positive scalar input, second parameter of beta distribution
#'
#' @return scalar value generated from beta distribution, rbeta(1, sum(Phi) + a_psi, d - sum(Phi) + b_psi)
#'
#' @examples
Generate_Psi = function(Phi, d, a_psi, b_psi){

  # Calculate sum of all entries of Phi matrix
  Phi_sum = sum(Phi)

  # Generate Psi from Beta distribution
  Psi = stats::rbeta(1, Phi_sum + a_psi, d - Phi_sum + b_psi)

  # Return Psi
  return(Psi)

}







#' Eta generating function
#'
#' @param b scalar input
#' @param phi scalar input either 0 or 1
#' @param a_eta positive scalar input, first parameter of inverse gamma distribution
#' @param b_eta positive scalar input, second parameter of inverse gamma distribution
#' @param nu_2 positive scalar input generated from inverse gamma distribution
#'
#' @return scalar value generated from inverse gamma distribution, phi / rgamma(1, a_eta + 1/2, b^2/2 + b_eta) + (1 - phi) / rgamma(1, a_eta + 1/2, b^2/(2 * nu_2) + b_eta)
#'
#' @examples
Generate_Eta = function(b, phi, a_eta, b_eta, nu_2){

  # Generate Eta based on Inverse Gamma distribution
  if(phi == 1){

    Eta = 1 / stats::rgamma(1, a_eta + 1/2, b^2/2 + b_eta)

  } else {

    Eta = 1 / stats::rgamma(1, a_eta + 1/2, b^2/(2 * nu_2) + b_eta)

  }

  # Return Eta
  return(Eta)

}






#' Tau generating function
#'
#' @param a scalar input
#' @param gamma scalar input either 0 or 1
#' @param a_tau positive scalar input, first parameter of inverse gamma distribution
#' @param b_tau positive scalar input, second parameter of inverse gamma distribution
#' @param nu_1 positive scalar input
#'
#' @return scalar value generated from inverse gamma distribution, gamma / rgamma(1, a_tau + 1/2, a^2/2 + b_tau) + (1 - gamma) / rgamma(1, a_tau + 1/2, a^2/(2 * nu_1) + b_tau)
#'
#' @examples
Generate_Tau = function(a, gamma, a_tau, b_tau, nu_1){


  # Generate Tau based on Inverse Gamma distribution
  if(gamma == 1){

    Tau = 1 / stats::rgamma(1, a_tau + 1/2, a^2/2 + b_tau)

  } else {

    Tau = 1 / stats::rgamma(1, a_tau + 1/2, a^2/(2 * nu_1) + b_tau)

  }

  # Return Tau
  return(Tau)

}





#' Sigma generating function
#'
#' @param n positive integer input
#' @param z_sum non-negative scalar input
#' @param a_sigma positive scalar input, first parameter of inverse gamma distribution
#' @param b_sigma positive scalar input, second parameter of inverse gamma distribution
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









#' Target function for a particular A entry and corresponding gamma
#'
#' @param X n * K matrix input
#' @param Y n * p matrix input
#' @param A p * p matrix input
#' @param a scalar input
#' @param N positive scalar input
#' @param Sigma_Inv p * 1 vector input
#' @param p positive integer input
#' @param B p * K matrix input
#' @param gamma scalar input, 0 or 1
#' @param tau positive scalar input
#' @param rho scalar input in between 0 and 1
#' @param nu_1 positive scalar input
#'
#' @return scalar target value corresponding to a particular a and gamma value
#'
#' @examples
Target_Agamma = function(X, Y, A, a, N, Sigma_Inv, p, B, gamma, tau, rho, nu_1){

  # Calculate (I_p - A)
  Mult_Mat = diag(p) - A

  # Calculate difference
  Diff = tcrossprod(Mult_Mat, Y) - tcrossprod(B, X)

  # Calculate Sum term inside exponential in likelihood
  Sum = sum(Sigma_Inv * rowSums(Diff^2))

  # Calculate target value
  Target = N * determinant(Mult_Mat, logarithm = TRUE)$modulus - Sum / 2 - gamma * (a^2 / (2 * tau)) - (1 - gamma) *  (a^2 / (2 * nu_1 * tau)) + gamma * log(rho) + (1 - gamma) * log(1 - rho)

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
#' @param prop_var2 postive scalar input
#'
#' @return scalar a value and scalar gamma value, 0 or 1
#'
#' @examples
Generate_Agamma = function(X, Y, A, i, j, Sigma_Inv, N, p, B, gamma, tau, rho, nu_1, prop_var1, prop_var2){

  # Value to update
  a = A[i, j]

  # Proposed value
  a_new = stats::rnorm(1, a, prop_var1)

  # New A matrix with proposed a value
  A_new = A
  A_new[i, j] = a_new

  # Calculate target values with a and a_new
  Target1 = Target_Agamma(X, Y, A_new, a_new, N, Sigma_Inv, p, B, 1-gamma, tau, rho, nu_1)
  Target2 = Target_Agamma(X, Y, A, a, N, Sigma_Inv, p, B, gamma, tau, rho, nu_1)

  # Calculate r
  r = Target1 - Target2


  # Generate uniform u
  u = stats::runif(1, 0, 1)

  # Compare u and r
  if(r >= log(u)){

    # Update a and gamma value and Target1
      a = a_new

      gamma = 1 - gamma

      Target2 = Target1

  }


  # Update a keeping gamma fixed
  A[i, j] = a

  # Proposed value
  a_new = stats::rnorm(1, a, prop_var2)

  # New A matrix with proposed a value
  A_new = A
  A_new[i, j] = a_new

  # Calculate target value with a_new
  Target1 = Target_Agamma(X, Y, A_new, a_new, N, Sigma_Inv, p, B, gamma, tau, rho, nu_1)


  # Calculate r
  r = Target1 - Target2


  # Generate uniform u
  u = stats::runif(1, 0, 1)

  # Compare u and r
  if(r >= log(u)){

    # Update a value
    a = a_new

  }


  # Return a and gamma value
  return(list(a = a, gamma = gamma))

}





#' Target function for a particular A entry and corresponding gamma
#'
#' @inheritParams Generate_Agamma
#' @param MultMat_Y p * n matrix input
#' @param b scalar input
#' @param phi scalar input, 0 or 1
#' @param eta positive scalar input
#' @param psi scalar input in between 0 and 1
#' @param nu_2 positive scalar input
#'
#' @return scalar target value
#'
#' @examples
Target_Bphi = function(X, Y, B, Sigma_Inv, MultMat_Y, b, phi, eta, psi, nu_2){

  # Calculate Difference vector
  Diff = MultMat_Y - tcrossprod(B, X)


  # Calculate Sum
  Sum = sum(rowSums(Diff^2) * Sigma_Inv)

  # Calculate Target value
  Target =  -Sum/2 - phi * (b^2 / (2 * eta)) - (1 - phi) * (b^2/(2 * nu_2 * eta)) + phi * log(psi) + (1 - phi) * log(1 - psi)


  # Return Target value
  return(Target)

}




#' B matrix entries and phi generating function
#'
#' @inheritParams Target_Bphi
#'
#' @param i positive integer input in between 1 and p
#' @param j positive integer input in between 1 and K
#' @param prop_var1 positive scalar input
#' @param prop_var2 positive scalar input
#'
#' @return scalar b value and scalar phi value, 0 or 1
#'
#' @examples
Generate_Bphi = function(X, Y, B, i, j, Sigma_Inv, MultMat_Y, phi, eta, psi, nu_2, prop_var1, prop_var2){

  # Value to update
  b = B[i, j]

  # Proposed value
  b_new = stats::rnorm(1, b, prop_var1)

  # New B matrix with proposed b value
  B_new = B
  B_new[i, j] = b_new

  # Calculate target values with b and b_new
  Target1 = Target_Bphi(X, Y, B_new, Sigma_Inv, MultMat_Y, b_new, 1 - phi, eta, psi, nu_2)
  Target2 = Target_Bphi(X, Y, B, Sigma_Inv, MultMat_Y, b, phi, eta, psi, nu_2)


  # Calculate r
  r = Target1 - Target2



  # Generate uniform u
  u = stats::runif(1, 0, 1)

  # Compare r with u
  if(r > log(u)){

    # Update b, phi and Target1
    b = b_new

    phi = 1 - phi

    Target2 = Target1

  }

  # Update b keeping phi fixed
  B[i, j] = b

  # Proposed value
  b_new = stats::rnorm(1, b, prop_var2)

  # New B matrix with proposed b value
  B_new = B
  B_new[i, j] = b_new

  # Calculate target value with b_new
  Target1 = Target_Bphi(X, Y, B_new, Sigma_Inv, MultMat_Y, b_new, phi, eta, psi, nu_2)


  # Calculate r
  r = Target1 - Target2


  # Generate uniform u
  u = stats::runif(1, 0, 1)

  # Compare u and r
  if(r >= log(u)){

    # Update b value
    b = b_new

  }

  # Return b and phi
  return(list(b = b, phi = phi))

}










#' Log-likelihood generating function
#'
#' @inheritParams Generate_Agamma
#' @param p positive integer input
#'
#' @return Log likelihood value corresponding to a particular A, B, X, Y and Sigma_Inv
#'
#' @examples
LL = function(A, B, X, Y, Sigma_Inv, p, N){

  # Calculate I_P - A
  Mult_Mat = diag(p) - A

  # Calculate difference as (I_p - A) %*% t(Y) - B %*% t(X)
  Diff = tcrossprod(Mult_Mat, Y) - tcrossprod(B, X)

  # Initiate Sum
  Sum = sum(Sigma_Inv * rowSums(Diff^2))

  # Calculate log-likelihood
  LL = N * determinant(Mult_Mat, logarithm = TRUE)$modulus - N / 2 * sum(log(1/Sigma_Inv)) - Sum / 2 - N / 2 * log(2 * pi)

  # Return LL
  return(LL)

}

















