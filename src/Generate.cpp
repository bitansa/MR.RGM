#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// Generate Rho
// [[Rcpp::export]]
double Generate_Rho_c(const arma::mat& Gamma, double p, double a_rho, double b_rho){

  // Calculate sum of all entries of Gamma matrix
  double  Gamma_sum = accu(Gamma);

  // Generate Rho from Beta distribution
  NumericVector RHO = Rcpp::rbeta(1, Gamma_sum + a_rho, p * (p - 1) - Gamma_sum + b_rho);

  // Store Rho
  double Rho = RHO(0);

  // Return Rho
  return Rho;

}


// Generate psi
// [[Rcpp::export]]
double Generate_Psi_c(const arma::mat& Phi, double d, double a_psi, double b_psi){

  // Calculate sum of all entries of Phi matrix
  double Phi_sum = accu(Phi);

  // Generate Psi from Beta distribution
  NumericVector PSI = Rcpp::rbeta(1, Phi_sum + a_psi, d - Phi_sum + b_psi);

  // Store Psi
  double Psi = PSI(0);

  // Return Psi
  return(Psi);

}



// Generate eta
// [[Rcpp::export]]
double Generate_Eta_c(double b, double phi, double a_eta, double b_eta, double nu_2){

  NumericVector eta;

  if(phi == 1){

    eta = 1 / Rcpp::rgamma(1, a_eta + 1 / double(2), 1 / (b * b / 2 + b_eta));

  } else {

    eta = 1 / Rcpp::rgamma(1, a_eta + 1 / double(2), 1 / (b * b / (2 * nu_2) + b_eta));

  }

  // Store eta
  double Eta = eta(0);


  // Return eta
  return(Eta);

}



// Generate tau
// [[Rcpp::export]]
double Generate_Tau_c(double a, double gamma, double a_tau, double b_tau, double nu_1){

  NumericVector tau;

  if(gamma == 1){

    tau = 1 / Rcpp::rgamma(1, a_tau + 1 / double(2), 1 / (a * a / 2 + b_tau));

  } else {

    tau = 1 / Rcpp::rgamma(1, a_tau + 1 / double(2), 1 / (a * a / (2 * nu_1) + b_tau));

  }

  // Store Tau
  double Tau = tau(0);


  // Return tau
  return(Tau);

}




// Generate sigma
// [[Rcpp::export]]
double Generate_Sigma_c(double n, double z_sum, double a_sigma, double b_sigma){

  // Generate sigma based on Inverse gamma distribution
  NumericVector sigma = 1 / Rcpp::rgamma(1, n/2 + a_sigma, 1/(z_sum/2 + b_sigma));

  // Store Sigma
  double Sigma = sigma(0);


  // Return sigma
  return(Sigma);

}



// Target function for A and gamma
// [[Rcpp::export]]
double Target_Agamma_c(const arma::mat& X, const arma::mat& Y, const arma::mat& A, const arma::mat& diag_p, double a, double N, const arma::colvec& Sigma_Inv, const arma::mat& B, double gamma, double tau, double rho, double nu_1){


  // Calculate (I_p - A)
  arma::mat Mult_Mat =  diag_p - A;

  // Calculate difference
  arma::mat Diff = Mult_Mat * Y.t() - B * X.t();

  // Calculate square difference
  Diff = Diff % Diff;

  // Initiate rowsum of difference^2
  arma::colvec Diff_sum = sum(Diff, 1);


  // Initiate Sum
  double Sum = std::inner_product(Sigma_Inv.begin(), Sigma_Inv.end(), Diff_sum.begin(), 0.0);

  // Calculate Target
  double Target = N * real(arma::log_det(Mult_Mat)) - Sum / 2 - gamma * (a * a/(2 * tau)) - (1 - gamma) *  (a * a/(2 * nu_1 * tau)) + gamma * log(rho) + (1 - gamma) * log(1 - rho);


  // Return Target
  return(Target);


}


// Generate an entry of A matrix and gamma
// [[Rcpp::export]]
NumericVector Generate_Agamma_C(const arma::mat& X, const arma::mat& Y, const arma::mat& A, const arma::mat& diag_p, double i, double j, const arma::colvec& Sigma_Inv, double N, const arma::mat& B, double gamma, double tau, double rho, double nu_1, double prop_var1){

  // Initiate Output
  NumericVector Output(2);

  // Value to update
  double a = A(i-1, j-1);


  // Proposed value
  NumericVector a_proposed = Rcpp::rnorm(1, a, std::sqrt(prop_var1));

  // Propose a_new
  double a_new = a_proposed(0);


  // New A matrix with proposed a value
  arma::mat A_new = A;

  // Update A_new
  A_new(i - 1, j - 1) = a_new;

  // Calculate r
  double r = Target_Agamma_c(X, Y, A_new, diag_p, a_new, N, Sigma_Inv, B, 1 - gamma, tau, rho, nu_1) - Target_Agamma_c(X, Y, A, diag_p, a, N, Sigma_Inv, B, gamma, tau, rho, nu_1);

  // Generate uniform u
  NumericVector rand_unif = Rcpp::runif(1, 0, 1);

  // Generate u
  double u = rand_unif(0);


  // Check whether r is big or not
  if(r > log(u)){

    a = a_new;

    gamma = 1 - gamma;

  }

  // Update Output
  Output(0) = a;
  Output(1) = gamma;


  // Return Output
  return(Output);

}


