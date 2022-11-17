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
