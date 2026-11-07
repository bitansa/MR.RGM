//////////////////////////////////////////////////////////////////////////////////////////////////////
// rgm_mcmc.h  --  MCMC building blocks shared by every RGM_* sampler
//
// 1. RWAdapt: adaptive random-walk Metropolis-within-Gibbs proposal scales, one per entry (Roberts and Rosenthal,
//    2009, "Examples of adaptive MCMC", J. Comput. Graph. Statist. 18, 349-367). During the burn-in only, after every
//    batch of `batch` iterations the log proposal standard deviation of each entry is increased by delta if the
//    entry's acceptance rate in that batch exceeded `target` (0.44, the optimal rate for one-dimensional random-walk
//    updates) and decreased by delta otherwise, delta = min(0.5, batch_number^(-1/2)) (Roberts and Rosenthal use min(0.01, .); the
//    larger early steps let the scales travel the orders of magnitude between small and large n within the burn-in). The scales are frozen at the end
//    of the burn-in, so the retained draws come from a fixed (non-adaptive) Markov kernel with the posterior as its
//    stationary distribution; no diminishing-adaptation argument is needed.
// 2. Initial proposal variances: 2.38^2 / (approximate posterior precision of the entry), the precision taken from
//    the quadratic term of the log-likelihood at the starting value (Gelman, Roberts and Gilks, 1996).
// 3. Gibbs_B_Row: exact Gibbs draw of the allowed entries of one row of B (jointly for up to 100 allowed instruments,
//    entry by entry within the row beyond that). Given A, Sigma (and C), the log-likelihood
//    is quadratic in B and the (spike-and-slab or slab) prior of every entry is normal given its indicator, so the
//    full conditional of b_j,L (L = allowed instruments of response j) is N(P^{-1} h, P^{-1}) with
//        P = n Omega_jj S_XX[L, L] + diag(1 / v_l),     v_l = eta_jl (phi_jl = 1) or nu_2 eta_jl (phi_jl = 0),
//        h = n [ (Omega T)_j,L - (Omega B S_XX)_j,L + Omega_jj (B_j S_XX)_L ],
//    where Omega = Sigma^{-1}, T = (I - A) S_YX - C S_XU^T (C = 0 without covariates), B_j = current row j of B.
// 4. Starting values: B from least squares of each response on its allowed instruments (A = 0), Sigma from the
//    residual second-moment matrix at these values.
//////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef MR_RGM_MCMC_H
#define MR_RGM_MCMC_H

#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>

namespace rgm {

// --------------------------------------------------------------------------------------------------
// 1. Adaptive random-walk proposal scales
// --------------------------------------------------------------------------------------------------
class RWAdapt {
public:
  arma::mat logsd;      // log proposal standard deviation of every entry
  arma::mat acc;        // acceptances of every entry in the current batch
  arma::umat mask;      // entries that are updated (1) or not (0)
  int batch, nbatch, nBurnin;
  double target;

  RWAdapt() : batch(50), nbatch(0), nBurnin(0), target(0.44) {}

  RWAdapt(const arma::mat& var0, const arma::umat& mask_, int nBurnin_, int batch_ = 50, double target_ = 0.44)
    : logsd(0.5 * arma::log(var0)), acc(arma::zeros(var0.n_rows, var0.n_cols)), mask(mask_),
      batch(batch_), nbatch(0), nBurnin(nBurnin_), target(target_) {}

  // Proposal variance of entry (i, j)
  double var(int i, int j) const { return std::exp(2.0 * logsd(i, j)); }

  // Record whether the proposal of entry (i, j) was accepted (only matters during the burn-in)
  void record(int i, int j, bool accepted) { if (accepted) acc(i, j) += 1.0; }

  // Call once at the end of iteration `iter` (1, ..., nIter): adapts after every batch during the burn-in
  void end_iteration(int iter) {
    if (iter > nBurnin || iter % batch != 0) return;
    nbatch++;
    double delta = std::min(0.5, 1.0 / std::sqrt(static_cast<double>(nbatch)));
    for (arma::uword c = 0; c < logsd.n_cols; c++) {
      for (arma::uword r = 0; r < logsd.n_rows; r++) {
        if (mask(r, c)) {
          logsd(r, c) += (acc(r, c) / batch > target) ? delta : -delta;
          // keep the scale within [1e-7, 10] (standard deviation)
          logsd(r, c) = std::min(std::max(logsd(r, c), std::log(1e-7)), std::log(10.0));
        }
      }
    }
    acc.zeros();
  }

  // Final (frozen) proposal variances; 0 for entries that are not updated
  arma::mat variances() const {
    arma::mat V = arma::exp(2.0 * logsd);
    V.elem(arma::find(mask == 0)).zeros();
    return V;
  }
};

// --------------------------------------------------------------------------------------------------
// 2. Initial proposal variances
//    Prop_Var > 0: that value for every entry (user-supplied); Prop_Var <= 0: automatic,
//    2.38^2 / (n * Omega_ii * S_jj + 1) for entry (i, j), with S = S_YY (entries of A) or S_XX (entries of B);
//    the "+ 1" is the order of the prior precision of a slab effect on standardized data.
// --------------------------------------------------------------------------------------------------
inline arma::mat Initial_PropVar(double Prop_Var, const arma::colvec& Omega_diag, const arma::colvec& S_diag, double n) {
  arma::mat V(Omega_diag.n_elem, S_diag.n_elem);
  for (arma::uword j = 0; j < S_diag.n_elem; j++) {
    for (arma::uword i = 0; i < Omega_diag.n_elem; i++) {
      V(i, j) = (Prop_Var > 0) ? Prop_Var : 5.6644 / (n * Omega_diag(i) * S_diag(j) + 1.0);
    }
  }
  return V;
}

// Entries of A that are updated: every off-diagonal entry
inline arma::umat OffDiag_Mask(int p) {
  arma::umat M = arma::ones<arma::umat>(p, p);
  M.diag().zeros();
  return M;
}

// --------------------------------------------------------------------------------------------------
// 3. Exact Gibbs draw of the allowed entries of row j of B (see the header)
//    T = (I - A) S_YX - C S_XU^T (p x k); Omega = Sigma^{-1} (p x p; diagonal for the diagonal error model);
//    prior variance of entry (j, l): Eta(j, l) if Phi(j, l) == 1, nu_2 * Eta(j, l) otherwise.
//    The standard normal draws use R's generator (R::norm_rand), so results follow set.seed().
// --------------------------------------------------------------------------------------------------
inline void Gibbs_B_Row(arma::mat& B, int j, const arma::mat& T, const arma::mat& S_XX, const arma::mat& Omega,
                        const arma::mat& D, const arma::mat& Phi, const arma::mat& Eta, double nu_2, double n,
                        arma::uword max_block = 100) {

  arma::uvec L = arma::find(D.row(j).t() != 0);
  if (L.n_elem == 0) return;

  // h = n [ (Omega T)_j,L - (Omega B S_XX)_j,L + Omega_jj (B_j S_XX)_L ]
  arma::rowvec OT  = Omega.row(j) * T;
  arma::rowvec OBS = Omega.row(j) * B * S_XX;
  double Ojj = Omega(j, j);

  // Prior variances of the allowed entries
  arma::colvec v(L.n_elem);
  for (arma::uword m = 0; m < L.n_elem; m++) {
    v(m) = (Phi(j, L(m)) == 1) ? Eta(j, L(m)) : nu_2 * Eta(j, L(m));
    v(m) = std::max(v(m), 1e-300);           // guard against an underflowed slab variance (1 / v would be infinite)
  }

  if (L.n_elem > max_block) {
    // Many allowed instruments (e.g. D = 1 with large k): exact entry-by-entry Gibbs within the row, cost O(|L| k)
    // instead of the O(|L|^3) block draw. Scalar full conditional of b_jl: precision n Omega_jj S_ll + 1 / v_l,
    // mean n [ (Omega T)_jl - (Omega B S_XX)_jl + Omega_jj b_jl S_ll ] / precision; (Omega B S_XX)_j,. is updated
    // after every change of b_jl by Omega_jj (b_new - b_old) S_XX[l, .].
    for (arma::uword m = 0; m < L.n_elem; m++) {
      arma::uword l = L(m);
      double b_old = B(j, l);
      double prec = n * Ojj * S_XX(l, l) + 1.0 / v(m);
      double h = n * (OT(l) - OBS(l) + Ojj * b_old * S_XX(l, l));
      double b_new = h / prec + R::norm_rand() / std::sqrt(prec);
      OBS += Ojj * (b_new - b_old) * S_XX.row(l);
      B(j, l) = b_new;
    }
    return;
  }

  // Block draw of the whole row (exact joint full conditional of b_j,L)
  arma::rowvec own = Ojj * (B.row(j) * S_XX);
  arma::rowvec hr  = n * (OT - OBS + own);
  arma::colvec h   = hr.elem(L);

  // P = n Omega_jj S_XX[L, L] + diag(1 / v)
  arma::mat P = n * Ojj * S_XX.submat(L, L);
  P.diag() += 1.0 / v;
  P = 0.5 * (P + P.t());

  // b = P^{-1} h + R^{-1} z,  P = R^T R (Cholesky, R upper triangular), z ~ N(0, I)
  arma::mat R;
  if (!arma::chol(R, P)) {
    P.diag() += 1e-10 * std::max(1.0, arma::max(P.diag()));
    if (!arma::chol(R, P)) Rcpp::stop("Gibbs update of B: the full-conditional precision matrix is not positive definite.");
  }
  arma::colvec z(L.n_elem);
  for (arma::uword m = 0; m < L.n_elem; m++) z(m) = R::norm_rand();
  // Triangular solves by substitution (solve_opts::fast: no condition-number check). P is often badly scaled -- an
  // entry in the spike has prior precision 1 / (nu_2 eta), orders of magnitude above the data precision -- which makes
  // its reciprocal condition number tiny although the triangular solve itself is backward stable; without `fast`
  // Armadillo would warn and switch to an approximate (least-squares) solution.
  arma::colvec w  = arma::solve(arma::trimatl(R.t()), h, arma::solve_opts::fast);
  arma::colvec mu = arma::solve(arma::trimatu(R), w, arma::solve_opts::fast);
  arma::colvec b  = mu + arma::solve(arma::trimatu(R), z, arma::solve_opts::fast);

  for (arma::uword m = 0; m < L.n_elem; m++) B(j, L(m)) = b(m);
}

// --------------------------------------------------------------------------------------------------
// 4. Starting values
// --------------------------------------------------------------------------------------------------
// B: least squares of each response on its allowed instruments (rows of D), with A = 0
inline arma::mat Start_B(const arma::mat& S_YX, const arma::mat& S_XX, const arma::mat& D) {
  int p = S_YX.n_rows, k = S_XX.n_cols;
  arma::mat B = arma::zeros(p, k);
  for (int j = 0; j < p; j++) {
    arma::uvec L = arma::find(D.row(j).t() != 0);
    if (L.n_elem == 0) continue;
    arma::mat G = S_XX.submat(L, L);
    G.diag() += 1e-8 * std::max(1.0, arma::mean(G.diag()));
    arma::rowvec r = S_YX.row(j);
    arma::colvec s = r.elem(L);
    arma::colvec b;
    if (!arma::solve(b, G, s, arma::solve_opts::no_approx)) b = arma::pinv(G) * s;
    for (arma::uword m = 0; m < L.n_elem; m++) B(j, L(m)) = b(m);
  }
  return B;
}

// Sigma: residual second-moment matrix S_YY - B S_YX^T - S_YX B^T + B S_XX B^T (A = 0), with its eigenvalues
// floored at 1e-3 times the largest so that it is positive definite
inline arma::mat Start_Sigma(const arma::mat& R0) {
  arma::mat R = 0.5 * (R0 + R0.t());
  arma::vec ev; arma::mat V;
  arma::eig_sym(ev, V, R);
  double floor_ = 1e-3 * std::max(ev.max(), 1e-8);
  ev = arma::clamp(ev, floor_, arma::datum::inf);
  return V * arma::diagmat(ev) * V.t();
}
inline arma::mat Start_Sigma(const arma::mat& S_YY, const arma::mat& S_YX, const arma::mat& S_XX, const arma::mat& B) {
  return Start_Sigma(S_YY - B * S_YX.t() - S_YX * B.t() + B * S_XX * B.t());
}

} // namespace rgm

#endif
