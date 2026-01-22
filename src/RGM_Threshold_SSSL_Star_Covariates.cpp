#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// For Inverse Wishart distribution
#include <RcppDist.h>
// [[Rcpp::depends(RcppDist)]]

// Include the GIGrvg library for the GIG distribution
#include <GIGrvg.h>
// [[Rcpp::depends(GIGrvg)]]



//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

// ------------------------------
// Helpers: file-local only
// ------------------------------
namespace {

// Sample Tau
double Sample_Tau(double a, double gamma, double tau, double nu_1) {

  // Sample Epsilon based on old tau
  double Epsilon = 1 / Rcpp::rgamma(1, 1, 1 / (1 + 1 / tau))(0);

  // Initialize Tau
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




// Function to update the Sigma matrix column by column when Sigma is full with spike and slab prior on it
arma::mat Sample_Sigma_Full(const arma::mat& Sigma, const arma::mat& S, const arma::mat& Z, double lambda,
                            double v0_square, double v1_square, double n) {

  // Get number of nodes
  int p = Sigma.n_cols;

  // Create a new sigma matrix
  arma::mat New_Sigma = Sigma;

  // Create matrix V based on Z
  arma::mat V = Z * v1_square + (1 - Z) * v0_square;

  // Set diagonal entries of V to 0
  V.diag().zeros();

  // --- get GIGrvg::rgig safely (no need to attach GIGrvg) ---
  Rcpp::Environment GIGrvg = Rcpp::Environment::namespace_env("GIGrvg");
  Rcpp::Function rgig = GIGrvg["rgig"];

  // Run a loop for each column of Sigma
  for (int i = 0; i < p; i++) {

    // Column indices apart from i
    arma::uvec not_col_index = arma::find(arma::linspace<arma::uvec>(1, p, p) != static_cast<unsigned int>(i + 1));


    // Create Submatrices from Sigma
    arma::mat Sigma11 = New_Sigma.submat(not_col_index, not_col_index);
    arma::colvec sigma12 = New_Sigma.submat(not_col_index, arma::uvec{static_cast<unsigned int>(i)});
    double sigma22 = New_Sigma(i, i);

    // Create Submatrices from S
    arma::mat S11 = S.submat(not_col_index, not_col_index);
    arma::colvec s12 = S.submat(not_col_index, arma::uvec{static_cast<unsigned int>(i)});
    double s22 = S(i, i);


    // Create Submatrices from V
    arma::mat V11 = V.submat(not_col_index, not_col_index);
    arma::colvec v12 = V.submat(not_col_index, arma::uvec{static_cast<unsigned int>(i)});

    // Calculate Sigmma_11 inverse
    arma::mat Sigma11_Inv = arma::inv(Sigma11);


    // Calculate the intermediate result: sigma12.t() * arma::inv(Sigma11) * sigma12
    arma::mat intermediate_result = sigma12.t() * Sigma11_Inv * sigma12;

    // Extract the double value from the matrix
    double intermediate_value = intermediate_result(0, 0);

    // Calculate v
    double v = sigma22 - intermediate_value;

    // Calculate diagonal matrix and B
    arma::mat diag_v12_inv = arma::diagmat(1 / v12);
    arma::mat B = (Sigma11_Inv * S11 * Sigma11_Inv) / v + lambda * Sigma11_Inv;
    arma::mat C = arma::inv(B + diag_v12_inv);

    // Make C symmetric
    C = 0.5 * (C + C.t());

    // Sample u from multivariate normal distribution
    arma::colvec u = arma::mvnrnd((C * Sigma11_Inv * s12) / v, C);

    // Sample v from GIG distribution
    Rcpp::NumericVector vv = rgig(Named("n") = 1, Named("lambda") = 1 - n / 2,
                                  Named("chi") = (arma::as_scalar(u.t() * Sigma11_Inv * S11 * Sigma11_Inv * u) -
                                  2 * arma::as_scalar(s12.t() * Sigma11_Inv * u) + s22), Named("psi") = lambda);
    v = vv[0];

    // Update New_Sigma
    New_Sigma.submat(not_col_index, arma::uvec{static_cast<unsigned int>(i)}) = u;
    New_Sigma.submat(arma::uvec{static_cast<unsigned int>(i)}, not_col_index) = u.t();
    New_Sigma(i, i) = v + arma::as_scalar(u.t() * Sigma11_Inv * u);

  }

  // Return New_Sigma
  return New_Sigma;

}


// Function to update the latent variable matrix Z
arma::mat Sample_Z(const arma::mat& Sigma, double pi, double v0_square, double v1_square) {

  // Number of nodes
  int p = Sigma.n_cols;

  // Initiate Z
  arma::mat Z = arma::ones(p, p);

  // Run a loop to update Z
  for (int i = 0; i < p - 1; ++i) {

    for (int j = i + 1; j < p; ++j) {

      // Get the corresponding Sigma Value
      double sigma_ij = Sigma(i, j);

      // Calculate the probability
      double prob = (1 / sqrt(v1_square) * exp(-0.5 * (sigma_ij * sigma_ij / v1_square)) * pi) /
        (1 / sqrt(v1_square) * exp(-0.5 * (sigma_ij * sigma_ij / v1_square)) * pi
           + 1 / sqrt(v0_square) * exp(-0.5 * (sigma_ij * sigma_ij / v0_square)) * (1 - pi));

      // Check if prob is NaN and assign 0 for it
      if (std::isnan(prob)) {
        prob = 0;  // Set prob to 0 if NaN is encountered
      }

      //double prob = R::dnorm(sigma_ij, 0, sqrt(v1_square), FALSE) * pi /
      //(R::dnorm(sigma_ij, 0, sqrt(v1_square), FALSE) * pi + R::dnorm(sigma_ij, 0, sqrt(v0_square), FALSE) * (1 - pi));

      // Sample from bernoulli distribution
      Z(i, j) = Rcpp::rbinom(1, 1, prob)(0);

      Z(j, i) = Z(i, j);

    }

  }

  // Return Z
  return Z;

}





// Sample random number from truncated normal
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
double tn_pdf(double x, double mu, double sigma, double a, double b) {

  // Calculate alpha and beta
  double alpha = (a - mu) / sigma;
  double beta = (b - mu) / sigma;

  // Calculate density
  double d = arma::normpdf(x, mu, sigma) / (arma::normcdf(beta, 0.0, 1.0) - arma::normcdf(alpha, 0.0, 1.0));

  // Return d
  return d;

}



// Calculate log-likelihood for the model when Y and U are there and sigma is full
double LL_Full_Star(const arma::mat& A, const arma::mat& C, const arma::mat& S_YY, const arma::mat& S_YU,
                    const arma::mat& S_UU, const arma::mat& Sigma_Inv, double p, double N) {

  // Calculate (I_p - A)
  const arma::mat& Mult_Mat = arma::eye(p, p) - A;

  // Calculate Sum
  double Sum = N * arma::trace(S_YY * Mult_Mat.t() * Sigma_Inv * Mult_Mat)
                + N * arma::trace(S_UU * C.t() * Sigma_Inv * C)
                - 2 * N * arma::trace(S_YU * C.t() * Sigma_Inv * Mult_Mat);

  // Calculate log-likelihood
  double LL = N * real(arma::log_det(Mult_Mat)) + N / 2 * real(arma::log_det(Sigma_Inv)) - Sum / 2 - N / 2 * log(2 * arma::datum::pi);

  // Return log-likelihood
  return LL;

}


// Calculate target value for a particular A for the model when Y and U are there
double Target_A_Star(double a, double N, double gamma, double tau, double nu_1, double Trace4, double Trace5, double Trace6, double Trace8, double logdet) {

  // Calculate Sum term inside exponential in likelihood
  double Sum = Trace4 + Trace5 + Trace6 + Trace8;

  // Calculate Target value
  double Target = N * logdet - Sum / 2 - gamma * (a * a / (2 * tau)) - (1 - gamma) * (0.5 * log(nu_1) + a * a / (2 * nu_1 * tau));

  // Return Target
  return Target;

}



// Sample a particular entry of matrix A for the model when Y and U are there and Sigma is a full matrix
Rcpp::List Sample_A_Full_Star(const arma::mat& S_YY, const arma::mat& S_YU, const arma::mat& A, const arma::mat& A_Pseudo, double i, double j,
                         const arma::mat& Sigma_Inv, double N, double p, const arma::mat& C, double gamma, double tau, double nu_1,
                         double prop_var1, double tA, double Trace4, double Trace5, double Trace6, double Trace8, arma::mat& InvMat, double logdet) {

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
  double Trace4_New = Trace4 - N * ((fabs(a_new) > tA) * a_new - (fabs(a) > tA) * a) * arma::trace(Sigma_Inv.row(i) * S_YY.col(j));
  double Trace5_New = Trace5 - N * ((fabs(a_new) > tA) * a_new - (fabs(a) > tA) * a) * arma::trace(S_YY.row(j) * Sigma_Inv.col(i));
  double Trace6_New = Trace6 + N * ((fabs(a_new) > tA) * a_new - (fabs(a) > tA) * a) * arma::trace(S_YY.row(j) * A.t() * Sigma_Inv.col(i)
                                                                                                     + S_YY.col(j).t() * A_new.t() * Sigma_Inv.row(i).t());
  double Trace8_New = Trace8 + 2 * N * ((fabs(a_new) > tA) * a_new - (fabs(a) > tA) * a) * arma::trace(S_YU.row(j) * C.t() * Sigma_Inv.col(i));

  // Calculate target values with a and a_new
  double Target1 = Target_A_Star(a_new, N, gamma, tau, nu_1, Trace4_New, Trace5_New, Trace6_New, Trace8_New, logdet_new);
  double Target2 = Target_A_Star(a, N, gamma, tau, nu_1, Trace4, Trace5, Trace6, Trace8, logdet);

  // Calculate r i.e. the differnce between two target values
  double r = Target1 - Target2;

  // Sample u from Uniform(0, 1)
  double u = Rcpp::runif(1, 0, 1)(0);

  // Compare u and r
  if (r >= log(u)) {

    // Update a, trace values, logdet and InvMat
    a = a_new;

    Trace4 = Trace4_New;
    Trace5 = Trace5_New;
    Trace6 = Trace6_New;
    Trace8 = Trace8_New;

    logdet = logdet_new;

    InvMat = InvMat - (((fabs(a) > tA) * a - (fabs(a_new) > tA) * a_new) / (1 + ((fabs(a) > tA) * a -
      (fabs(a_new) > tA) * a_new) * InvMat(j, i))) * (InvMat.col(i) * InvMat.row(j));

  }

  // Return a, trace values, logdet and InvMat
  return Rcpp::List::create(Rcpp::Named("a") = a, Rcpp::Named("Trace4") = Trace4,
                            Rcpp::Named("Trace5") = Trace5, Rcpp::Named("Trace6") = Trace6,
                            Rcpp::Named("Trace8") = Trace8,
                            Rcpp::Named("logdet") = logdet, Rcpp::Named("InvMat") = InvMat);


}



// Sample from matrix variate normal
arma::mat rmatrixvarnorm(const arma::mat& M, const arma::mat& U, const arma::mat& V) {
  Rcpp::RNGScope scope; // Ensures R's RNG is used
  int n = M.n_rows;
  int p = M.n_cols;

  // Cholesky decomposition of U (row covariance) and V (column covariance)
  arma::mat U_chol = arma::chol(U);
  arma::mat V_chol = arma::chol(V);

  // Generate matrix Z of independent standard normal random variables
  arma::mat Z = arma::randn(n, p); // Z ~ N(0, I)

  // Compute the sample from matrix-variate normal distribution: X = M + U_chol * Z * V_chol.t()
  arma::mat X = M + U_chol.t() * Z * V_chol;

  return X;
}

} // end anonymous namespace

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////





// Do MCMC sampling with threshold prior for the model when Y and U are there and Sigma is full with spike and slab prior on it
// [[Rcpp::export]]
Rcpp::List RGM_Threshold_SSSL_Star_Covariates(const arma::mat& S_YY, const arma::mat& S_YU, const arma::mat& S_UU,
                                         double n, int nIter, int nBurnin,
                                         int Thin, double nu_1 = 0.0001,
                                         double Prop_VarA = 0.01, double TAU = 2){


  // Calculate number of nodes from S_YY matrix
  int p = S_YY.n_cols;

  // Calculate number of columns of S_UU
  int l = S_UU.n_cols;

  // Take b_Tilde = p + 1 (b_Tilde should be greater than p-1; p+1 is taken to ensure mean is not infinity)
  int b_Tilde = p + 1;

  // Initialize A, C and A_Pseudo matrices
  arma::mat A = arma::zeros(p, p);
  arma::mat C = arma::zeros(p, l);
  arma::mat A_Pseudo = arma::zeros(p, p);

  // Initialize Sigma
  arma::mat Sigma = riwish(b_Tilde, arma::eye(p, p));

  // Initialize Gamma and Tau matrices
  arma::mat Gamma = arma::ones(p, p);
  arma::mat Tau = arma::ones(p, p);

  // Initialize Z
  arma::mat Z = arma::ones(p, p);

  // Make the diagonals of Gamma, Tau and Z matrices to be 0
  Gamma.diag().zeros();
  Tau.diag().zeros();
  Z.diag().zeros();

  // Initialize tA, t0 and t_sd
  double tA = 0;
  double t0 = 0.1;
  double t_sd = 0.1;

  // Initialize acceptance counter
  double AccptA = 0;
  double Accpt_tA = 0;

  // Calculate number of posterior samples
  int nPst = std::floor((nIter - nBurnin) / Thin);

  // Initiate Itr to index the posterior samples
  int Itr = 0;

  // Initialize posterior arrays and matrices
  arma::cube A_Pst = arma::zeros(p, p, nPst);
  arma::cube A0_Pst = arma::zeros(p, p, nPst);
  arma::cube C_Pst = arma::zeros(p, l, nPst);
  arma::cube Gamma_Pst = arma::zeros(p, p, nPst);
  arma::cube Tau_Pst = arma::zeros(p, p, nPst);
  arma::colvec tA_Pst = arma::zeros(nPst);
  arma::cube Sigma_Pst = arma::zeros(p, p, nPst);
  arma::cube Z_Pst = arma::zeros(p, p, nPst);

  // Initialize LogLikelihood vector
  arma::colvec LL_Pst = arma::zeros(nPst);

  // Calculate Sigma_Inv
  arma::mat Sigma_Inv = arma::inv(Sigma);

  // Parameters for spike and slab on Sigma (v0 >= 0.01, v1 = hv0, h <= 1000, lambda = 1/5/10, Pi = Inclusion probability)
  double v0_square = 0.005;
  double v1_square = 5;
  double Pi = 0.5;
  double lambda = 1;

  // Caculate the sub matrix to sample C
  arma::mat SubMat_C = arma::inv(n * S_UU + arma::eye(l, l) / TAU);

  // Run a loop to do MCMC sampling
  for(int i = 1; i <= nIter; i++){

    // Calculate I_p - A
    arma::mat MultMat = arma::eye(p, p) - A;


    ////////////////////
    // Update C
    arma::mat Mean_C = (n * MultMat * S_YU) * SubMat_C;
    C = rmatrixvarnorm(Mean_C, Sigma, SubMat_C);


    ////////////////////
    // Change here based on C
    arma::mat S = n * MultMat * S_YY * MultMat.t()
      + n * C * S_UU * C.t() - n * MultMat * S_YU * C.t() - n * C * S_YU.t() * MultMat.t()
      + C * C.t() / TAU;

      // Update Sigma
      Sigma = Sample_Sigma_Full(Sigma, S, Z, lambda, v0_square, v1_square, n);

      // Calculate Sigma_Inv
      Sigma_Inv = arma::inv(Sigma);

      // Sample Z
      Z = Sample_Z(Sigma, Pi, v0_square, v1_square);


      ////////////////////
      // Update A
      // Calculate trace values
      double Trace4 = - n * arma::trace(S_YY * A.t() * Sigma_Inv);
      double Trace5 = - n * arma::trace(S_YY * Sigma_Inv * A);
      double Trace6 = n * arma::trace(S_YY * A.t() * Sigma_Inv * A);
      double Trace8 = 2 * n * arma::trace(S_YU * C.t() * Sigma_Inv * A);

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
            Rcpp::List Output2 = Sample_A_Full_Star(S_YY, S_YU, A, A_Pseudo, j, l, Sigma_Inv, n, p, C, 1, Tau(j, l), nu_1,
                                               Prop_VarA, tA, Trace4, Trace5, Trace6, Trace8, InvMat, logdet);


            double a = Output2[0];


            // Update Acceptance counter
            if (A_Pseudo(j, l) != a) {

              // Increase AccptA
              AccptA = AccptA + 1;

              // Update Trace values
              Trace4 = Output2[1];
              Trace5 = Output2[2];
              Trace6 = Output2[3];
              Trace8 = Output2[4];

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
      double Diff_A = LL_Full_Star(A_new, C, S_YY, S_YU, S_UU, Sigma_Inv, p, n) - LL_Full_Star(A, C, S_YY, S_YU, S_UU, Sigma_Inv, p, n)
        + log(tn_pdf(tA, tA_new, t_sd, 0, t0)) - log(tn_pdf(tA_new, tA, t_sd, 0, t0));


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
        C_Pst.slice(Itr) = C;
        Gamma_Pst.slice(Itr) = Gamma;
        Tau_Pst.slice(Itr) = Tau;
        tA_Pst(Itr) = tA;
        Sigma_Pst.slice(Itr) = Sigma;
        Z_Pst.slice(Itr) = Z;
        LL_Pst(Itr) = LL_Full_Star(A, C, S_YY, S_YU, S_UU, Sigma_Inv, p, n);

        // Increase Itr by 1
        Itr = Itr + 1;

      }

      // Check if `i` is divisible by 100 and print progress
      //if (i % 100 == 0) {
      //Rcpp::Rcout << "Iterations " << i << " done." << std::endl;
      //}


  }

  // Calculate estimates based on posterior samples
  arma::mat A_Est = mean(A_Pst, 2);
  arma::mat A0_Est = mean(A0_Pst, 2);
  arma::mat C_Est = mean(C_Pst, 2);
  arma::mat Gamma_Est = mean(Gamma_Pst, 2);
  arma::mat Tau_Est = mean(Tau_Pst, 2);
  double tA_Est = mean(tA_Pst);
  arma::mat Sigma_Est = mean(Sigma_Pst, 2);
  arma::mat Z_Est = mean(Z_Pst, 2);

  // Construct the graph structures
  arma::umat logicalGraph_A = (Gamma_Est > 0.5);
  arma::mat zA_Est = arma::conv_to<arma::mat>::from(logicalGraph_A);


  // Return outputs
  return Rcpp::List::create(Rcpp::Named("A_Est") = A_Est,
                            Rcpp::Named("C_Est") = C_Est,
                            Rcpp::Named("zA_Est") = zA_Est,
                            Rcpp::Named("A0_Est") = A0_Est,
                            Rcpp::Named("Gamma_Est") = Gamma_Est, Rcpp::Named("Tau_Est") = Tau_Est,
                            Rcpp::Named("tA_Est") = tA_Est,
                            Rcpp::Named("Sigma_Est") = Sigma_Est, Rcpp::Named("Z_Est") = Z_Est,
                            Rcpp::Named("AccptA") = AccptA / (p * (p - 1) * nIter) * 100,
                            Rcpp::Named("Accpt_tA") = Accpt_tA / (nIter) * 100,
                            Rcpp::Named("LL_Pst") = LL_Pst, Rcpp::Named("Gamma_Pst") = Gamma_Pst);



}
