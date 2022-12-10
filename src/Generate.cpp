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



// Target function for B and phi
// [[Rcpp::export]]
double Target_Bphi_c(const arma::mat& X, const arma::mat& Y, const arma::mat& B, const arma::colvec& Sigma_Inv, const arma::mat& MultMat_Y, double b, double phi, double eta, double psi, double nu_2){


  // Calculate Diff vector
  arma::mat Diff = MultMat_Y - B * X.t();


  // Calculate square difference
  Diff = Diff % Diff;

  // Initiate rowsum of difference^2
  arma::colvec Diff_sum = sum(Diff, 1);


  // Initiate Sum
  double Sum = std::inner_product(Sigma_Inv.begin(), Sigma_Inv.end(), Diff_sum.begin(), 0.0);



  // Calculate Target
  double Target =  - Sum / 2 - phi * (b * b / (2 * eta)) - (1 - phi) * (b * b /(2 * nu_2 * eta)) + phi * log(psi) + (1 - phi) * log(1 - psi);


  // Return Target
  return(Target);

}







// Generate an entry of B matrix and phi
// [[Rcpp::export]]
NumericVector Generate_Bphi_c(const arma::mat& X, const arma::mat& Y, const arma::mat& B, double i, double j, const arma::colvec& Sigma_Inv, const arma::mat& MultMat_Y, double phi, double eta, double psi, double nu_2, double prop_var2){


  // Initiate Output
  NumericVector Output(2);

  // Value to update
  double b = B(i-1, j-1);


  // Proposed value
  NumericVector b_proposed = Rcpp::rnorm(1, b, std::sqrt(prop_var2));

  // Propose b_new
  double b_new = b_proposed(0);

  // New B matrix with proposed b value
  arma::mat B_new = B;

  // Update B_new
  B_new(i - 1, j - 1) = b_new;

  // Calculate r
  double r = Target_Bphi_c(X, Y, B_new, Sigma_Inv, MultMat_Y, b_new, 1 - phi, eta, psi, nu_2) - Target_Bphi_c(X, Y, B, Sigma_Inv, MultMat_Y, b, phi, eta, psi, nu_2);


  // Generate uniform u
  NumericVector rand_unif = Rcpp::runif(1, 0, 1);

  // Generate u
  double u = rand_unif(0);


  // Check whether r is big or not
  if(r > log(u)){

    b = b_new;

    phi = 1 - phi;

  }

  // Update Output
  Output(0) = b;
  Output(1) = phi;


  return(Output);


}


// Calculate log-likelihood
// [[Rcpp::export]]
double LL_c(const arma::mat& A, const arma::mat& B, const arma::mat& X, const arma::mat& Y, const arma::colvec& Sigma_Inv, const arma::mat& diag_p, double N){

  // Calculate I_P - A
  arma::mat Mult_Mat = diag_p - A;

  // Calculate difference
  arma::mat Diff = Mult_Mat * Y.t() - B * X.t();

  // Calculate square difference
  Diff = Diff % Diff;

  // Initiate rowsum of difference^2
  arma::colvec Diff_sum = sum(Diff, 1);

  // Initiate Sum
  double Sum = std::inner_product(Sigma_Inv.begin(), Sigma_Inv.end(), Diff_sum.begin(), 0.0);


  // Calculate log-likelihood
  double LL = N * real(arma::log_det(Mult_Mat)) - N / 2 * accu(log(1/Sigma_Inv)) - Sum / 2;

  // Return LL
  return(LL);

}



// Generate A and B matrix
// [[Rcpp::export]]
Rcpp::List Get_AB_c(const arma::mat& A0, const arma::mat& B0, const arma::mat& X, const arma::mat& Y, const arma::mat& D, double d, const arma::mat& diag_p, double a_tau = 0.1, double b_tau = 0.1, double a_rho = 0.5, double b_rho = 0.5, double nu_1 = 0.0001, double a_eta = 0.1, double b_eta = 0.1, double a_psi = 0.5, double b_psi = 0.5, double nu_2 = 0.0001, double  a_sigma = 0.1, double b_sigma = 0.1, double Prop_varA = 5, double Prop_VarB = 5, double niter = 25000){


  // Calculate number of datapoints and number of nodes from Y matrix
  int n = Y.n_rows;
  int p = Y.n_cols;

  // Calculate number of columns of X
  int K = X.n_cols;

  // Initiate matrix A and B
  arma::mat A = A0;
  arma::mat B = B0;

  // Initiate Sigma_Inv
  arma::colvec Sigma_Inv = Rcpp::rgamma(p, a_sigma, 1 / b_sigma);

  // Initialize A_Update and B_Update matrix
  arma::mat A_Update = arma::zeros(p * p, niter);
  arma::mat B_Update = arma::zeros(p * K, niter);

  // Initialize scalar Rho, Psi
  double Rho = Rcpp::rbeta(1, a_rho, b_rho)(0);
  double Psi = Rcpp::rbeta(1, a_psi, b_psi)(0);

  // Initialize Gamma and Phi matrix
  arma::mat Gamma = arma::zeros(p, p);
  arma::mat Phi = arma::zeros(p, K);
  arma::mat Tau = arma::zeros(p, p);
  arma::mat Eta = arma::zeros(p, K);

  // Update Gamma and tau
  for(int i = 0; i < p; i++){

    for(int j = 0; j < p; j++){

      if(i != j){

        Gamma(i, j) = Rcpp::rbinom(1, 1, Rho)(0);

        Tau(i, j) = 1 / Rcpp::rgamma(1, a_tau, 1 / b_tau)(0);

      }
    }

  }


  // Update Phi and eta
  for(int i = 0; i < p; i++){

    for(int j = 0; j < K; j++){

      if(D(i, j) != 0){

        Phi(i, j) = Rcpp::rbinom(1, 1, Psi)(0);

        Eta(i, j) = 1 / Rcpp::rgamma(1, a_eta, 1/b_eta)(0);

      }
    }

  }


  // Initialize Gamma_Update and Phi_Update matrix
  arma::mat Gamma_Update = arma::zeros(p * p, niter);
  arma::mat Phi_Update = arma::zeros(p * K, niter);

  // Initialize LL_1
  arma::colvec LL_1 = arma::zeros(niter);


  // Run a loop to get the updates
  for(int i = 1; i <= niter; i++){

    // Update Rho
    Rho = Generate_Rho_c(Gamma, p, a_rho, b_rho);

    // Update Psi
    Psi = Generate_Psi_c(Phi, d, a_psi, b_psi);


    // Update Tau, A and gamma
    for(int j = 0; j < p; j++){

      for(int l = 0; l < p; l++){

        if(j != l){

          // Calculate Tau
          Tau(j, l) = Generate_Tau_c(A(j, l), Gamma(j, l), a_tau, b_tau, nu_1);

          // Generate a and gamma jointly
          NumericVector Out = Generate_Agamma_C(X, Y, A, diag_p, j + 1, l + 1, Sigma_Inv, n, B, Gamma(j, l), Tau(j, l), Rho, nu_1, Prop_varA);

          // Update A matrix
          A(j, l) = Out(0);

          // Update A_Update matrix
          A_Update(p * l + j, i - 1) = Out(0);

          // Updata Gamma matrix
          Gamma(j, l) = Out(1);

          // Update Gamma_Update mmatrix
          Gamma_Update(p * l + j, i - 1) = Out(1);


        }
      }

    }


    // Calculate MultMat_Y as (I_p - A) %*% t(Y)
    arma::mat MultMat_Y = (diag_p - A) * Y.t();


    // Update Eta, B and Phi
    for(int j = 0; j < p; j++){

      for(int l = 0; l < K; l++){

        if(D(j, l) != 0){

          // Calculate Eta
          Eta(j, l) = Generate_Eta_c(B(j, l), Phi(j, l), a_eta, b_eta, nu_2);

          // Generate b and phi jointly
          NumericVector Out = Generate_Bphi_c(X, Y, B, j + 1, l + 1, Sigma_Inv, MultMat_Y, Phi(j, l), Eta(j, l), Psi, nu_2, Prop_VarB);

          // Update B matrix
          B(j, l) = Out(0);

          // Update B_Update matrix
          B_Update(p * l + j, i - 1) = Out(0);

          // Updata Phi matrix
          Phi(j, l) = Out(1);

          // Updata Phi_Update matrix
          Phi_Update(p * l + j, i - 1) = Out(1);


        }
      }

    }


    // Calculate I_P - A
    arma::mat Mult_Mat = diag_p - A;

    // Calculate difference
    arma::mat Diff = Mult_Mat * Y.t() - B * X.t();

    // Calculate square difference
    Diff = Diff % Diff;

    // Initiate rowsum of difference^2
    arma::colvec Diff_sum = sum(Diff, 1);

    // Update Sigma_inv
    for(int j = 0; j < p; j++){

      Sigma_Inv(j) = 1 / Generate_Sigma_c(n, Diff_sum(j), a_sigma, b_sigma);

    }


    // Update LL_1
    LL_1(i - 1) = LL_c(A, B, X, Y, Sigma_Inv, diag_p, n);


  }

  return Rcpp::List::create(Rcpp::Named("A_Update") = A_Update, Rcpp::Named("B_Update") = B_Update, Rcpp::Named("Gamma_Update") = Gamma_Update, Rcpp::Named("Phi_Update") = Phi_Update, Rcpp::Named("LL_1") = LL_1);


}

