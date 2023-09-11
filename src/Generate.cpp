#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// Sample Rho
// [[Rcpp::export]]
double Sample_Rho(double Gamma, double a_rho, double b_rho) {

  // Sample Rho from beta distribution
  double Rho = Rcpp::rbeta(1, Gamma + a_rho, 1 - Gamma + b_rho)(0);

  // Return Rho
  return Rho;

}


// Sample Psi
// [[Rcpp::export]]
double Sample_Psi(double Phi, double a_psi, double b_psi) {

  // Sample Psi from beta distribution
  double Psi = Rcpp::rbeta(1, Phi + a_psi, 1 - Phi + b_psi)(0);

  // Return Psi
  return Psi;

}


// Sample Tau
// [[Rcpp::export]]
double Sample_Tau(double a, double gamma, double tau, double nu_1) {

  // Sample Epsilon based on old tau
  double Epsilon = 1 / Rcpp::rgamma(1, 1, 1 / (1 + 1 / tau))(0);

  // Sample Tau from inverse gamma distribution
  double Tau;

  // Check whether gamma is 0 or 1
  if (gamma == 1) {

    // Sample Tau based on a and Epsilon
    Tau = 1 / Rcpp::rgamma(1, 1, 1 / (a * a / 2 + 1 / Epsilon))(0);

  } else {

    // Sample Tau based on a, nu_1 and Epsilon
    Tau = 1 / Rcpp::rgamma(1, 1, 1 / (a * a / (2 * nu_1) + 1 / Epsilon))(0);

  }

  // Return Tau
  return Tau;

}



// Sample Eta
// [[Rcpp::export]]
double Sample_Eta(double b, double phi, double eta, double nu_2) {

  // Sample Epsilon based on old eta
  double Epsilon = 1 / Rcpp::rgamma(1, 1, 1 / (1 + 1 / eta))(0);

  // Sample Eta from inverse gamma distribution
  double Eta;

  // Check whether phi is 0 or 1
  if (phi == 1) {

    // Sample Eta based on b and Epsilon
    Eta = 1 / Rcpp::rgamma(1, 1, 1 / (b * b / 2 + 1 / Epsilon))(0);

  } else {

    // Sample Eta based on b, nu_2 and Epsilon
    Eta = 1 / Rcpp::rgamma(1, 1, 1 / (b * b / (2 * nu_2) + 1 / Epsilon))(0);

  }

  // Return Eta
  return Eta;

}



// Sample Gamma
// [[Rcpp::export]]
double Sample_Gamma(double a, double tau, double rho, double nu_1) {

  // Calculate acceptance probability
  double p = exp(-0.5 * (a * a / tau)) * rho / (exp(-0.5 * (a * a / tau)) * rho + 1 / sqrt(nu_1) * exp(-0.5 * (a * a / (nu_1 * tau))) * (1 - rho));

  // Sample Gamma from binomial distribution
  double Gamma = Rcpp::rbinom(1, 1, p)(0);

  // Return Gamma
  return Gamma;

}



// Sample Phi
// [[Rcpp::export]]
double Sample_Phi(double b, double eta, double psi, double nu_2) {

  // Calculate acceptance probability
  double p = exp(-0.5 * (b * b / eta)) * psi / (exp(-0.5 * (b * b / eta)) * psi + 1 / sqrt(nu_2) * exp(-0.5 * (b * b / (nu_2 * eta))) * (1 - psi));

  // Sample Phi from binomial distribution
  double Phi = Rcpp::rbinom(1, 1, p)(0);

  // Return Phi
  return Phi;

}


// Sample Sigma
// [[Rcpp::export]]
double Sample_Sigma(double n, double z_sum, double a_sigma, double b_sigma) {

  // Sample Sigma from inverse gamma distribution
  double Sigma = 1.0 / Rcpp::rgamma(1, n / 2.0 + a_sigma, 1.0 / (z_sum / 2.0 + b_sigma))(0);

  // Return Sigma
  return Sigma;

}


// Calculate target value for a particular A
// [[Rcpp::export]]
double Target_A(double a, double N, double gamma, double tau, double nu_1, double Trace3, double Trace4, double Trace5, double Trace6, double logdet) {

  // Calculate Sum term inside exponential in likelihood
  double Sum = Trace3 + Trace4 + Trace5 + Trace6;

  // Calculate Target value
  double Target = N * logdet - Sum / 2 - gamma * (a * a / (2 * tau)) - (1 - gamma) * (0.5 * log(nu_1) + a * a / (2 * nu_1 * tau));

  // Return Target
  return Target;

}


// Sample a particular entry of matrix A
// [[Rcpp::export]]
Rcpp::List Sample_A(const arma::mat& S_YY, const arma::mat& S_YX, const arma::mat& A, const arma::mat& A_Pseudo, double i, double j, const arma::colvec& Sigma_Inv, double N, double p, const arma::mat& B, double gamma, double tau, double nu_1, double prop_var1, double tA, double Trace3, double Trace4, double Trace5, double Trace6, arma::mat InvMat, double logdet) {

  // Value to update
  double a = A_Pseudo(i, j);

  // Proposed value
  double a_new = Rcpp::rnorm(1, a, sqrt(prop_var1))(0);

  // Create a copy of matrix A
  arma::mat A_new = A;

  // Modify the copy with the proposed a value
  A_new(i, j) = (fabs(a_new) > tA) * a_new;

  // Modify logdet
  double logdet_new = logdet + log(fabs(1 + ((fabs(a) > tA) * a - (fabs(a_new) > tA) * a_new) * InvMat(j, i)));


  // Calculate new trace values
  double Trace3_New = Trace3 - N * ((fabs(a_new) > tA) * a_new - (fabs(a) > tA) * a) * Sigma_Inv(i) * S_YY(i, j);
  double Trace4_New = Trace4 - N * ((fabs(a_new) > tA) * a_new - (fabs(a) > tA) * a) * Sigma_Inv(i) * S_YY(i, j);
  double Trace5_New = Trace5 + N * ((fabs(a_new) > tA) * a_new - (fabs(a) > tA) * a) * Sigma_Inv(i) * arma::trace(A.row(i) * S_YY.col(j) + S_YY.row(j) * A_new.row(i).t());
  double Trace6_New = Trace6 + 2 * N * ((fabs(a_new) > tA) * a_new - (fabs(a) > tA) * a) * Sigma_Inv(i) * arma::trace(B.row(i) * S_YX.row(j).t());

  // Calculate target values with a and a_new
  double Target1 = Target_A(a_new, N, gamma, tau, nu_1, Trace3_New, Trace4_New, Trace5_New, Trace6_New, logdet_new);
  double Target2 = Target_A(a, N, gamma, tau, nu_1, Trace3, Trace4, Trace5, Trace6, logdet);

  // Calculate r i.e. the differnce between two target values
  double r = Target1 - Target2;

  // Sample u from Uniform(0, 1)
  double u = Rcpp::runif(1, 0, 1)(0);

  // Compare u and r
  if (r >= log(u)) {

    // Update a, trace values, logdet and InvMat
    a = a_new;

    Trace3 = Trace3_New;
    Trace4 = Trace4_New;
    Trace5 = Trace5_New;
    Trace6 = Trace6_New;

    logdet = logdet_new;

    InvMat = InvMat - (((fabs(a) > tA) * a - (fabs(a_new) > tA) * a_new) / (1 + ((fabs(a) > tA) * a - (fabs(a_new) > tA) * a_new) * InvMat(j, i))) * (InvMat.col(i) * InvMat.row(j));

  }

  // Return a, trace values, logdet and InvMat
  return Rcpp::List::create(Rcpp::Named("a") = a, Rcpp::Named("Trace3") = Trace3,
                            Rcpp::Named("Trace4") = Trace4, Rcpp::Named("Trace5") = Trace5,
                            Rcpp::Named("Trace6") = Trace6, Rcpp::Named("logdet") = logdet,
                            Rcpp::Named("InvMat") = InvMat);

}



// Calculate target value for a particular B
// [[Rcpp::export]]
double Target_B(double b, double phi, double eta, double nu_2, double Trace1, double Trace2) {

  // Calculate Sum
  double Sum = Trace1 + Trace2;

  // Calculate Target value
  double Target = - Sum / 2 - phi * (b * b / (2 * eta)) - (1 - phi) * (0.5 * log(nu_2) + b * b / (2 * nu_2 * eta));

  // Return Target value
  return Target;

}


// Sample a particular entry of matrix B
// [[Rcpp::export]]
Rcpp::List Sample_B(const arma::mat& S_YX, const arma::mat& S_XX, const arma::mat& B, const arma::mat& B_Pseudo, double i, double j, const arma::colvec& Sigma_Inv, const arma::mat& MultMat, double N, double phi, double eta, double nu_2, double prop_var2, double tB, double Trace1, double Trace2) {

  // Value to update
  double b = B_Pseudo(i, j);

  // Proposed value
  double b_new = Rcpp::rnorm(1, b, sqrt(prop_var2))(0);

  // Create a copy of matrix B
  arma::mat B_new = B;

  // Modify the copy with the proposed b value
  B_new(i, j) = (fabs(b_new) > tB) * b_new;

  // Calculate new trace values
  double Trace1_New = Trace1 - 2 * N * arma::trace(Sigma_Inv(i) * ((fabs(b_new) > tB) * b_new - (fabs(b) > tB) * b) * (MultMat.row(i) * S_YX.col(j)));
  double Trace2_New = Trace2 + N * arma::trace(Sigma_Inv(i) * ((fabs(b_new) > tB) * b_new - (fabs(b) > tB) * b) * (B.row(i) * S_XX.col(j) + S_XX.row(j) * B_new.row(i).t()));

  // Calculate target values with b and b_new
  double Target1 = Target_B(b_new, phi, eta, nu_2, Trace1_New, Trace2_New);
  double Target2 = Target_B(b, phi, eta, nu_2, Trace1, Trace2);

  // Calculate r i.e. the difference between two target values
  double r = Target1 - Target2;

  // Sample u from Uniform(0, 1)
  double u = Rcpp::runif(1, 0, 1)(0);

  // Compare u and r
  if (r >= log(u)) {

    // Update b and trace values
    b = b_new;

    Trace1 = Trace1_New;
    Trace2 = Trace2_New;

  }

  // Return b and trace values
  return Rcpp::List::create(Rcpp::Named("b") = b, Rcpp::Named("Trace1") = Trace1,
                            Rcpp::Named("Trace2") = Trace2);

}



// Sample random number from truncated normal
// [[Rcpp::export]]
double Sample_tn(double mu, double sigma, double a, double b) {

  // Calculate alpha and beta
  double alpha = (a - mu) / sigma;
  double beta = (b - mu) / sigma;

  // Calculate CDF
  double cdf_alpha = arma::normcdf(alpha, 0.0, 1.0);
  double cdf_beta = arma::normcdf(beta, 0.0, 1.0);

  // Sample from truncated normal with mean mu and sd sigma
  double u = Rcpp::runif(1, cdf_alpha, cdf_beta)(0);
  double x = R::qnorm(u, 0.0, 1.0, true, false) * sigma + mu;

  // Return x
  return x;

}


// Pdf of truncated normal distribution
// [[Rcpp::export]]
double tn_pdf(double x, double mu, double sigma, double a, double b) {

  // Calculate alpha and beta
  double alpha = (a - mu) / sigma;
  double beta = (b - mu) / sigma;

  // Calculate density
  double d = arma::normpdf(x, mu, sigma) / (arma::normcdf(beta, 0.0, 1.0) - arma::normcdf(alpha, 0.0, 1.0));

  // Return d
  return d;

}


// Calculate log-likelihood
// [[Rcpp::export]]
double LL(const arma::mat& A, const arma::mat& B, const arma::mat& S_YY, const arma::mat& S_YX, const arma::mat& S_XX, const arma::colvec& Sigma_Inv, double p, double N) {

  // Calculate (I_p - A)
  const arma::mat& Mult_Mat = arma::eye(p, p) - A;

  // Calculate Sum
  double Sum = N * arma::trace(S_YY * Mult_Mat.t() * arma::diagmat(Sigma_Inv) * Mult_Mat) - 2 * N * arma::trace(S_YX * B.t() * arma::diagmat(Sigma_Inv) * Mult_Mat)
                   + N * arma::trace(S_XX * B.t() * arma::diagmat(Sigma_Inv) * B);

  // Calculate log-likelihood
  double LL = N * real(arma::log_det(Mult_Mat)) - N / 2 * accu(log(1/Sigma_Inv)) - Sum / 2 - N / 2 * log(2 * arma::datum::pi);

  // Return log-likelihood
  return LL;

}







// Do MCMC sampling with threshold prior
// [[Rcpp::export]]
Rcpp::List RGM_Threshold(const arma::mat& S_YY, const arma::mat& S_YX, const arma::mat& S_XX, const arma::mat& D, double n, int nIter, int nBurnin, int Thin, double nu_1 = 0.0001, double nu_2 = 0.0001, double a_sigma = 0.01, double b_sigma = 0.01, double Prop_VarA = 0.01, double Prop_VarB = 0.01){


  // Calculate number of nodes from S_YY matrix
  int p = S_YY.n_cols;

  // Calculate number of columns of S_XX
  int k = S_XX.n_cols;

  // Initialize A, B, A_Pseudo and B_Pseudo matrix
  arma::mat A = arma::zeros(p, p);
  arma::mat B = arma::zeros(p, k);
  arma::mat A_Pseudo = arma::zeros(p, p);
  arma::mat B_Pseudo = arma::zeros(p, k);

  // Initialize Sigma_Inv
  arma::colvec Sigma_Inv = Rcpp::rgamma(p, a_sigma, 1 / b_sigma);


  // Initialize Gamma, Phi, Tau and Eta matrix
  arma::mat Gamma = arma::ones(p, p);
  arma::mat Phi = arma::ones(p, k);
  arma::mat Tau = arma::ones(p, p);
  arma::mat Eta = arma::ones(p, k);

  // Make the diagonals of Gamma and Tau matrix to be 0
  Gamma.diag().zeros();
  Tau.diag().zeros();

  // Make Phi[i, j] = 0 and Eta[i, j] = 0 if D[i, j] = 0
  Phi = Phi % D;
  Eta = Eta % D;

  // Initialize tA, tB, t0 and t_sd
  double tA = 0;
  double tB = 0;
  double t0 = 1;
  double t_sd = 0.1;

  // Initialize acceptance counter
  double AccptA = 0;
  double AccptB = 0;
  double Accpt_tA = 0;
  double Accpt_tB = 0;

  // Calculate number of posterior samples
  int nPst = std::floor((nIter - nBurnin) / Thin);

  // Initiate Itr to index the posterior samples
  int Itr = 0;

  // Initialize posterior arrays and matrices
  arma::cube A_Pst = arma::zeros(p, p, nPst);
  arma::cube A0_Pst = arma::zeros(p, p, nPst);
  arma::cube B_Pst = arma::zeros(p, k, nPst);
  arma::cube B0_Pst = arma::zeros(p, k, nPst);
  arma::cube Gamma_Pst = arma::zeros(p, p, nPst);
  arma::cube Tau_Pst = arma::zeros(p, p, nPst);
  arma::cube Phi_Pst = arma::zeros(p, k, nPst);
  arma::cube Eta_Pst = arma::zeros(p, k, nPst);
  arma::colvec tA_Pst = arma::zeros(nPst);
  arma::colvec tB_Pst = arma::zeros(nPst);
  arma::cube Sigma_Pst = arma::zeros(1, p, nPst);

  // Initialize LogLikelihood vector
  arma::colvec LL_Pst = arma::zeros(nPst);


  // Run a loop to do MCMC sampling
  for(int i = 1; i <= nIter; i++){

    // Update B
    // Calculate I_p - A
    arma::mat MultMat = arma::eye(p, p) - A;

    // Calculate Trace
    double Trace1 = - 2 * n * arma::trace(S_YX * B.t() * arma::diagmat(Sigma_Inv) * MultMat);
    double Trace2 = n * arma::trace(S_XX * B.t() * arma::diagmat(Sigma_Inv) * B);

    // Update Eta based on corresponding b  and then Update b based on the corresponding eta
    for (int j = 0; j < p; j++) {

      for (int l = 0; l < k; l++) {

        // Don't update if the corresponding D entry is 0
        if (D(j, l) != 0) {

          // Sample Eta
          Eta(j, l) = Sample_Eta(B_Pseudo(j, l), 1, Eta(j, l), nu_2);

          // Sample b
          Rcpp::List Output1 = Sample_B(S_YX, S_XX, B, B_Pseudo, j, l, Sigma_Inv, MultMat, n, 1, Eta(j, l), nu_2, Prop_VarB, tB, Trace1, Trace2);
          double b = Output1[0];

          // Update acceptance counter
          if (B_Pseudo(j, l) != b) {

            // Increase AccptB
            AccptB = AccptB + 1;

            // Update Trace Values
            Trace1 = Output1[1];
            Trace2 = Output1[2];

          }

          // Update B_Pseudo, B and Phi
          B_Pseudo(j, l) = b;

          B(j, l) = (std::abs(B_Pseudo(j, l)) > tB) * B_Pseudo(j, l);

          Phi(j, l) = (std::abs(B_Pseudo(j, l)) > tB) * 1;

        }

      }

    }

    // Propose tB_new
    double tB_new = Sample_tn(tB, t_sd, 0, t0);

    // Create B_new
    arma::mat B_new = B_Pseudo % (arma::abs(B_Pseudo) > tB_new);

    // Calculate difference
    double Diff = LL(A, B_new, S_YY, S_YX, S_XX, Sigma_Inv, p, n) - LL(A, B, S_YY, S_YX, S_XX, Sigma_Inv, p, n) + log(tn_pdf(tB, tB_new, t_sd, 0, t0)) - log(tn_pdf(tB_new, tB, t_sd, 0, t0));

    // Compare Diff with log of a random number from Uniform(0, 1)
    if (Diff > log(Rcpp::runif(1, 0, 1)(0))) {

      // Update B, tB and Accpt_tB
      B = B_new;

      tB = tB_new;

      Accpt_tB = Accpt_tB + 1;

    }


    ////////////////////
    // Update Sigma
    for (int j = 0; j < p; j++) {

      // Calculate Sum
      double z_sum = n * arma::accu(MultMat.row(j) * S_YY * MultMat.row(j).t()) - 2 * n * arma::accu(MultMat.row(j) * S_YX * B.row(j).t()) + n * arma::accu(B.row(j) * S_XX * B.row(j).t());

      // Sample Sigma_Inv
      Sigma_Inv(j) = 1 / Sample_Sigma(n, z_sum, a_sigma, b_sigma);

    }


    ////////////////////
    // Update A
    // Calculate trace values
    double Trace3 = - n * arma::trace(S_YY * A.t() * arma::diagmat(Sigma_Inv));
    double Trace4 = - n * arma::trace(S_YY * arma::diagmat(Sigma_Inv) * A);
    double Trace5 = n * arma::trace(S_YY * A.t() * arma::diagmat(Sigma_Inv) * A);
    double Trace6 = 2 * n * arma::trace(S_YX * B.t() * arma::diagmat(Sigma_Inv) * A);

    // Calculate det(I - A) and (I - A)^(-1)
    double logdet = real(arma::log_det(MultMat));
    arma::mat InvMat = arma::inv(MultMat);

    // Update Tau based on a and then update a based on Tau
    for (int j = 0; j < p; j++) {

      for (int l = 0; l < p; l++) {

        // Don't update the diagonal entries
        if (l != j) {

          // Update Tau
          Tau(j, l) = Sample_Tau(A_Pseudo(j, l), 1, Tau(j, l), nu_1);

          // Sample a
          Rcpp::List Output2 = Sample_A(S_YY, S_YX, A, A_Pseudo, j, l, Sigma_Inv, n, p, B, 1, Tau(j, l), nu_1, Prop_VarA, tA, Trace3, Trace4, Trace5, Trace6, InvMat, logdet);
          double a = Output2[0];

          // Update Acceptance counter
          if (A_Pseudo(j, l) != a) {

            // Increase AccptA
            AccptA = AccptA + 1;

            // Update Trace values
            Trace3 = Output2[1];
            Trace4 = Output2[2];
            Trace5 = Output2[3];
            Trace6 = Output2[4];

            // Update logdet and (I - A)^(-1)
            logdet = Output2[5];
            InvMat = Rcpp::as<arma::mat>(Output2[6]);


          }

          // Update A_Pseudo, A and Gamma
          A_Pseudo(j, l) = a;

          A(j, l) = (std::abs(A_Pseudo(j, l)) > tA) * A_Pseudo(j, l);

          Gamma(j, l) = (std::abs(A_Pseudo(j, l)) > tA) * 1;

        }

      }

    }

    // Propose tA_new
    double tA_new = Sample_tn(tA, t_sd, 0, t0);

    // Create A_new
    arma::mat A_new = A_Pseudo % (arma::abs(A_Pseudo) > tA_new);

    // Calculate Difference
    double Diff_A = LL(A_new, B, S_YY, S_YX, S_XX, Sigma_Inv, p, n) - LL(A, B, S_YY, S_YX, S_XX, Sigma_Inv, p, n) + log(tn_pdf(tA, tA_new, t_sd, 0, t0)) - log(tn_pdf(tA_new, tA, t_sd, 0, t0));

    // Compare Diff with log of a random number from uniform(0, 1)
    if (Diff_A > log(Rcpp::runif(1, 0, 1)(0))) {

      // Update A, tA and Accpt_tA
      A = A_new;

      tA = tA_new;

      Accpt_tA = Accpt_tA + 1;

    }


    // Store posterior samples
    if((i > nBurnin) && (i % Thin == 0)){

      A_Pst.slice(Itr) = A;
      A0_Pst.slice(Itr) = A_Pseudo;
      B_Pst.slice(Itr) = B;
      B0_Pst.slice(Itr) = B_Pseudo;
      Gamma_Pst.slice(Itr) = Gamma;
      Tau_Pst.slice(Itr) = Tau;
      Phi_Pst.slice(Itr) = Phi;
      Eta_Pst.slice(Itr) = Eta;
      tA_Pst(Itr) = tA;
      tB_Pst(Itr) = tB;
      Sigma_Pst.slice(Itr) = 1 / Sigma_Inv.t();
      LL_Pst(Itr) = LL(A, B, S_YY, S_YX, S_XX, Sigma_Inv, p, n);

      // Increase Itr by 1
      Itr = Itr + 1;

    }


  }

  // Calculate estimates based on posterior samples
  arma::mat A_Est = mean(A_Pst, 2);
  arma::mat B_Est = mean(B_Pst, 2);
  arma::mat A0_Est = mean(A0_Pst, 2);
  arma::mat B0_Est = mean(B0_Pst, 2);
  arma::mat Gamma_Est = mean(Gamma_Pst, 2);
  arma::mat Tau_Est = mean(Tau_Pst, 2);
  arma::mat Phi_Est = mean(Phi_Pst, 2);
  arma::mat Eta_Est = mean(Eta_Pst, 2);
  double tA_Est = mean(tA_Pst);
  double tB_Est = mean(tB_Pst);
  arma::mat Sigma_Est = mean(Sigma_Pst, 2);

  // Construct the graph structures
  arma::umat logicalGraph_A = (Gamma_Est > 0.5);
  arma::umat logicalGraph_B = (Phi_Est > 0.5);
  arma::mat zA_Est = arma::conv_to<arma::mat>::from(logicalGraph_A);
  arma::mat zB_Est = arma::conv_to<arma::mat>::from(logicalGraph_B);


  // Return outputs
  return Rcpp::List::create(Rcpp::Named("A_Est") = A_Est, Rcpp::Named("B_Est") = B_Est,
                            Rcpp::Named("zA_Est") = zA_Est, Rcpp::Named("zB_Est") = zB_Est,
                            Rcpp::Named("A0_Est") = A0_Est, Rcpp::Named("B0_Est") = B0_Est,
                            Rcpp::Named("Gamma_Est") = Gamma_Est, Rcpp::Named("Tau_Est") = Tau_Est,
                            Rcpp::Named("Phi_Est") = Phi_Est, Rcpp::Named("Eta_Est") = Eta_Est,
                            Rcpp::Named("tA_Est") = tA_Est, Rcpp::Named("tB_Est") = tB_Est,
                            Rcpp::Named("Sigma_Est") = Sigma_Est,
                            Rcpp::Named("AccptA") = AccptA / (p * (p - 1) * nIter) * 100, Rcpp::Named("AccptB") = AccptB / (arma::accu(D) * nIter) * 100,
                            Rcpp::Named("Accpt_tA") = Accpt_tA / (nIter) * 100, Rcpp::Named("Accpt_tB") = Accpt_tB / (nIter) * 100,
                            Rcpp::Named("LL_Pst") = LL_Pst);



}






// Do MCMC sampling with Spike and Slab Prior
// [[Rcpp::export]]
Rcpp::List RGM_SpikeSlab(const arma::mat& S_YY, const arma::mat& S_YX, const arma::mat& S_XX, const arma::mat& D, double n, int nIter, int nBurnin, int Thin, double a_rho = 3.0, double b_rho = 1.0, double nu_1 = 0.001, double a_psi = 0.5, double b_psi = 0.5, double nu_2 = 0.0001, double a_sigma = 0.01, double b_sigma = 0.01, double Prop_VarA = 0.01, double Prop_VarB = 0.01){


  // Calculate number of nodes from S_YY matrix
  int p = S_YY.n_cols;

  // Calculate number of columns of S_XX
  int k = S_XX.n_cols;

  // Initialize matrix A, B
  arma::mat A = arma::zeros(p, p);
  arma::mat B = arma::zeros(p, k);

  // Initialize Sigma_Inv
  arma::colvec Sigma_Inv = Rcpp::rgamma(p, a_sigma, 1 / b_sigma);

  // Initialize Rho, Psi, Gamma, Phi, Tau and Eta matrix
  arma::mat Rho = arma::zeros(p, p);
  arma::mat Psi = arma::zeros(p, k);
  arma::mat Gamma = arma::ones(p, p);
  arma::mat Phi = arma::ones(p, k);
  arma::mat Tau = arma::ones(p, p);
  arma::mat Eta = arma::ones(p, k);

  // Make the diagonals of Gamma and Tau matrix to be 0
  Gamma.diag().zeros();
  Tau.diag().zeros();

  // Make Phi[i, j] = 0 and Eta[i, j] = 0 if D[i, j] = 0
  Phi = Phi % D;
  Eta = Eta % D;

  // Initialize acceptance counter
  double AccptA = 0;
  double AccptB = 0;

  // Calculate number of posterior samples
  int nPst = std::floor((nIter - nBurnin) / Thin);

  // Initialize Itr to index the posterior samples
  int Itr = 0;

  // Initialize posterior arrays and matrices
  arma::cube A_Pst = arma::zeros(p, p, nPst);
  arma::cube B_Pst = arma::zeros(p, k, nPst);
  arma::cube Gamma_Pst = arma::zeros(p, p, nPst);
  arma::cube Tau_Pst = arma::zeros(p, p, nPst);
  arma::cube Rho_Pst = arma::zeros(p, p, nPst);
  arma::cube Phi_Pst = arma::zeros(p, k, nPst);
  arma::cube Eta_Pst = arma::zeros(p, k, nPst);
  arma::cube Psi_Pst = arma::zeros(p, k, nPst);
  arma::cube Sigma_Pst = arma::zeros(1, p, nPst);

  // Initialize LogLikelihood vector
  arma::colvec LL_Pst = arma::zeros(nPst);


  // Run a loop to do MCMC sampling
  for(int i = 1; i <= nIter; i++){

    // Update B
    // Calculate I_p - A
    arma::mat MultMat = arma::eye(p, p) - A;

    // Calculate trace values
    double Trace1 = -2 * n * arma::trace(S_YX * B.t() * arma::diagmat(Sigma_Inv) * MultMat);
    double Trace2 = n * arma::trace(S_XX * B.t() * arma::diagmat(Sigma_Inv) * B);

    // Update Psi, Eta, Phi and b
    for (int j = 0; j < p; j++) {

      for (int l = 0; l < k; l++) {

        // Don't update if the corresponding D entry is 0
        if (D(j, l) != 0) {

          // Sample Psi
          Psi(j, l) = Sample_Psi(Phi(j, l), a_psi, b_psi);

          // Sample Eta
          Eta(j, l) = Sample_Eta(B(j, l), Phi(j, l), Eta(j, l), nu_2);

          // Sample Phi
          Phi(j, l) = Sample_Phi(B(j, l), Eta(j, l), Psi(j, l), nu_2);

          // Sample b
          Rcpp::List Output1 = Sample_B(S_YX, S_XX, B, B, j, l, Sigma_Inv, MultMat, n, Phi(j, l), Eta(j, l), nu_2, Prop_VarB, -1, Trace1, Trace2);
          double b = Output1[0];

          // Update acceptance counter
          if (B(j, l) != b) {

            // Increase AccptB
            AccptB = AccptB + 1;

            // Update Trace Values
            Trace1 = Output1[1];
            Trace2 = Output1[2];

          }

          // Update B
          B(j, l) = b;

        }

      }

    }



    ////////////////////
    // Update Sigma
    for (int j = 0; j < p; j++) {

      // Calculate Sum
      double z_sum = n * arma::accu(MultMat.row(j) * S_YY * MultMat.row(j).t()) - 2 * n * arma::accu(MultMat.row(j) * S_YX * B.row(j).t()) + n * arma::accu(B.row(j) * S_XX * B.row(j).t());

      // Sample Sigma_Inv
      Sigma_Inv(j) = 1 / Sample_Sigma(n, z_sum, a_sigma, b_sigma);

    }



    ////////////////////
    // Update A
    // Calculate trace values
    double Trace3 = - n * arma::trace(S_YY * A.t() * arma::diagmat(Sigma_Inv));
    double Trace4 = - n * arma::trace(S_YY * arma::diagmat(Sigma_Inv) * A);
    double Trace5 = n * arma::trace(S_YY * A.t() * arma::diagmat(Sigma_Inv) * A);
    double Trace6 = 2 * n * arma::trace(S_YX * B.t() * arma::diagmat(Sigma_Inv) * A);

    // Calculate logdet and (I - A)^(-1)
    double logdet = real(arma::log_det(MultMat));
    arma::mat InvMat = arma::inv(MultMat);

    // Update Rho, Tau, Gamma and a
    for (int j = 0; j < p; j++) {

      for (int l = 0; l < p; l++) {

        // Don't update the diagonal entries
        if (l != j) {

          // Sample Rho
          Rho(j, l) = Sample_Rho(Gamma(j, l), a_rho, b_rho);

          // Sample Tau
          Tau(j, l) = Sample_Tau(A(j, l), Gamma(j, l), Tau(j, l), nu_1);

          // Sample Gamma
          Gamma(j, l) = Sample_Gamma(A(j, l), Tau(j, l), Rho(j, l), nu_1);

          // Sample a
          Rcpp::List Output2 = Sample_A(S_YY, S_YX, A, A, j, l, Sigma_Inv, n, p, B, Gamma(j, l), Tau(j, l), nu_1, Prop_VarA, -1, Trace3, Trace4, Trace5, Trace6, InvMat, logdet);
          double a = Output2[0];

          // Update acceptance counter
          if (A(j, l) != a) {

            // Increase AccptA
            AccptA = AccptA + 1;

            // Update trace values
            Trace3 = Output2[1];
            Trace4 = Output2[2];
            Trace5 = Output2[3];
            Trace6 = Output2[4];

            // Update logdet and (I - A)^(-1)
            logdet = Output2[5];
            InvMat = Rcpp::as<arma::mat>(Output2[6]);

          }

          // Update A
          A(j, l) = a;

        }

      }

    }


    // Store posterior samples
    if((i > nBurnin) && (i % Thin == 0)){

      A_Pst.slice(Itr) = A;
      B_Pst.slice(Itr) = B;
      Gamma_Pst.slice(Itr) = Gamma;
      Tau_Pst.slice(Itr) = Tau;
      Rho_Pst.slice(Itr) = Rho;
      Phi_Pst.slice(Itr) = Phi;
      Eta_Pst.slice(Itr) = Eta;
      Psi_Pst.slice(Itr) = Psi;
      Sigma_Pst.slice(Itr) = 1 / Sigma_Inv.t();
      LL_Pst(Itr) = LL(A, B, S_YY, S_YX, S_XX, Sigma_Inv, p, n);

      // Increase Itr by 1
      Itr = Itr + 1;

    }


  }

  // Calculate estimates based on posterior samples
  arma::mat A_Est = mean(A_Pst, 2);
  arma::mat B_Est = mean(B_Pst, 2);
  arma::mat Gamma_Est = mean(Gamma_Pst, 2);
  arma::mat Tau_Est = mean(Tau_Pst, 2);
  arma::mat Rho_Est = mean(Rho_Pst, 2);
  arma::mat Phi_Est = mean(Phi_Pst, 2);
  arma::mat Eta_Est = mean(Eta_Pst, 2);
  arma::mat Psi_Est = mean(Psi_Pst, 2);
  arma::mat Sigma_Est = mean(Sigma_Pst, 2);

  // Construct the graph structures
  arma::umat logicalGraph_A = (Gamma_Est > 0.5);
  arma::umat logicalGraph_B = (Phi_Est > 0.5);
  arma::mat zA_Est = arma::conv_to<arma::mat>::from(logicalGraph_A);
  arma::mat zB_Est = arma::conv_to<arma::mat>::from(logicalGraph_B);

  // Return outputs
  return Rcpp::List::create(Rcpp::Named("A_Est") = A_Est, Rcpp::Named("B_Est") = B_Est,
                            Rcpp::Named("zA_Est") = zA_Est, Rcpp::Named("zB_Est") = zB_Est,
                            Rcpp::Named("Gamma_Est") = Gamma_Est, Rcpp::Named("Tau_Est") = Tau_Est,
                            Rcpp::Named("Rho_Est") = Rho_Est, Rcpp::Named("Phi_Est") = Phi_Est,
                            Rcpp::Named("Eta_Est") = Eta_Est, Rcpp::Named("Psi_Est") = Psi_Est,
                            Rcpp::Named("Sigma_Est") = Sigma_Est,
                            Rcpp::Named("AccptA") = AccptA / (p * (p - 1) * nIter) * 100, Rcpp::Named("AccptB") = AccptB / (arma::accu(D) * nIter) * 100,
                            Rcpp::Named("LL_Pst") = LL_Pst);


}
