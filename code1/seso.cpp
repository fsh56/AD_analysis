#include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Implementation of my_rinvgamma
NumericVector my_rinvgamma(int n, double shape, double rate) {
  Function ff("rinvgamma", Environment::namespace_env("invgamma"));
  NumericVector res = ff(n, Named("shape") = shape, _["rate"] = rate);
  return res;
}

// Implementation of my_rinvwishart
// [[Rcpp::export]]
arma::mat my_rinvwishart(double nu, arma::mat S) {
  return arma::iwishrnd(S, nu);
}

// Helper function for multivariate normal sampling with robust Cholesky
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  sigma = 0.5 * (sigma + sigma.t());
  arma::mat L;
  bool success = arma::chol(L, sigma, "lower");
  if (!success) {
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, sigma);
    eigval = arma::clamp(eigval, 1e-8, arma::datum::inf);
    L = eigvec * arma::diagmat(arma::sqrt(eigval));
  }
  return arma::repmat(mu, 1, n).t() + Y * L.t();
}

// [[Rcpp::export]]
Rcpp::NumericVector my_rdirichlet(int n, Rcpp::NumericVector alpha) {
  if (n != 1) Rcpp::stop("my_rdirichlet currently supports n = 1 only.");
  int K = alpha.size();
  Rcpp::NumericVector draws(K);
  double sum = 0.0;
  for (int k = 0; k < K; ++k) {
    draws[k] = R::rgamma(alpha[k], 1.0);
    sum += draws[k];
  }
  // Normalize
  for (int k = 0; k < K; ++k) draws[k] /= sum;
  return draws;
}

// [[Rcpp::export]]
NumericVector gibbs_seso_uhp_only(int niter,
                                  NumericVector Gamma_hat,
                                  NumericVector gamma_hat,
                                  NumericVector s_hat_Gamma,
                                  NumericVector s_hat_gamma) {
  
  NumericVector s2_hat_Gamma = pow(s_hat_Gamma,2);
  NumericVector s2_hat_gamma = pow(s_hat_gamma,2);
  // define K
  int K = gamma_hat.size();
  // store beta
  NumericVector beta_tk(niter);
  // MCMC trace
  NumericMatrix theta_tk(niter, K);
  NumericMatrix gamma_tk(niter, K);
  NumericVector sigma2_gamma_tk(niter);
  NumericVector sigma2_theta_tk(niter);
  
  // setup priors based on input data
  double a_gamma = std::max(2.0, K / 4.0);
  double a_theta = a_gamma;
  double b_gamma = std::max(1e-3, (var(Gamma_hat) - pow(mean(s_hat_Gamma), 2)))*(a_gamma-1);
  double b_theta = std::max(1e-5, (var(gamma_hat) - pow(mean(s_hat_gamma), 2)))*(a_theta-1);
  
  // initialize starting values
  double beta_cur = 0.0;
  NumericVector theta_cur(K);
  NumericVector gamma_cur(K);
  double sigma2_gamma_cur = 1.0;
  double sigma2_theta_cur = 1.0;
  
  for (int k = 0; k < K; k++) {
    theta_cur[k] = rnorm(1, 0, 0.1)[0];
    gamma_cur[k] = rnorm(1, 0, 0.1)[0];
  }
  beta_tk[0] = beta_cur;
  theta_tk(0, _) = theta_cur;
  gamma_tk(0, _) = gamma_cur;
  sigma2_gamma_tk[0] = sigma2_gamma_cur;
  sigma2_theta_tk[0] = sigma2_theta_cur;
  
  // Gibbs
  for (int iter = 0; iter < (niter-1); iter++) {
    // update gamma_k
    for (int k = 0; k < K; k++) {
      double Ak_gamma = beta_cur * beta_cur / s2_hat_Gamma[k] + 1.0 / s2_hat_gamma[k] + 1.0 / sigma2_gamma_cur;
      double Bk_gamma = beta_cur * (Gamma_hat[k] - theta_cur[k]) / s2_hat_Gamma[k] + gamma_hat[k] / s2_hat_gamma[k];
      gamma_cur[k] = rnorm(1, Bk_gamma / Ak_gamma, sqrt(1.0 / Ak_gamma))[0];
    }
    gamma_tk(iter + 1, _) = gamma_cur;
    
    // update theta_k
    for (int k = 0; k < K; k++) {
      double Ak_theta = 1.0 / s2_hat_Gamma[k] + 1.0 / sigma2_theta_cur;
      double Bk_theta = (Gamma_hat[k] - beta_cur * gamma_cur[k]) / s2_hat_Gamma[k];
      theta_cur[k] = rnorm(1, Bk_theta / Ak_theta, sqrt(1.0 / Ak_theta))[0];
    }
    theta_tk(iter + 1, _) = theta_cur;
    
    // update sigma2_gamma and sigma2_theta
    sigma2_gamma_cur = my_rinvgamma(1, a_gamma + K / 2.0, b_gamma + 0.5 * sum(pow(gamma_cur, 2)))[0];
    sigma2_theta_cur = my_rinvgamma(1, a_theta + K / 2.0, b_theta + 0.5 * sum(pow(theta_cur, 2)))[0];
    
    // add bounds to prevent numerical issues
    if (sigma2_gamma_cur < 1e-10) sigma2_gamma_cur = 1e-10;
    if (sigma2_gamma_cur > 100) sigma2_gamma_cur = 10;
    if (sigma2_theta_cur < 1e-10) sigma2_theta_cur = 1e-10;
    if (sigma2_theta_cur > 100) sigma2_theta_cur = 10;
    
    sigma2_gamma_tk[iter + 1] = sigma2_gamma_cur;
    sigma2_theta_tk[iter + 1] = sigma2_theta_cur;
    
    // update beta
    vec U_beta = Gamma_hat - theta_cur;
    vec W_beta = gamma_cur;
    mat Omega_hat_Gamma(K, K, fill::zeros);
    for (int ii = 0; ii < K; ii++) {
      Omega_hat_Gamma(ii, ii) = 1.0 / s2_hat_Gamma[ii];
    }
    
    mat A_beta_mat = trans(W_beta) * Omega_hat_Gamma * W_beta;
    double A_beta = std::max(A_beta_mat(0,0), 1e-8);
    mat temp = trans(W_beta) * Omega_hat_Gamma * U_beta;
    double mu_beta = temp(0,0) / A_beta;
    beta_cur = rnorm(1, mu_beta, sqrt(1.0 / A_beta))[0];
    beta_tk[iter + 1] = beta_cur;
  }
  
  return beta_tk;
}

