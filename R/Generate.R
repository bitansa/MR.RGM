# Generate rho
Generate_Rho = function(Gamma, p, a_rho, b_rho){

  # Calculate sum of all entries of Gamma matrix
  Gamma_sum = sum(Gamma)

  # Generate Rho from Beta distribution
  Rho = rbeta(1, Gamma_sum + a_rho, p * (p - 1) - Gamma_sum + b_rho)

  # Return Rho
  return(Rho)

}


# Generate psi
Generate_Psi = function(Phi, d, a_psi, b_psi){

  # Calculate sum of all entries of Phi matrix
  Phi_sum = sum(Phi)

  # Generate psi from Beta distribution
  Psi = rbeta(1, Phi_sum + a_psi, d - Phi_sum + b_psi)

  # Return Rho
  return(Psi)

}


# Generate eta
Generate_Eta = function(b, phi, a_eta, b_eta, nu_2){

  # Generate eta based on Inverse Gamma distribution
  eta = phi / rgamma(1, a_eta + 1/2, b^2/2 + b_eta) + (1 - phi) / rgamma(1, a_eta + 1/2, b^2/(2 * nu_2) + b_eta)

  # Return eta
  return(eta)

}


# Generate tau
Generate_Tau = function(a, gamma, a_tau, b_tau, nu_1){

  # Generate tau based on Inverse Gamma distribution
  tau = gamma / rgamma(1, a_tau + 1/2, a^2/2 + b_tau) + (1 - gamma) / rgamma(1, a_tau + 1/2, a^2/(2 * nu_1) + b_tau)

  # Return tau
  return(tau)

}
