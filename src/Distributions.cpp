#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

const double log2pi = std::log(2.0 * M_PI);

arma::mat rmvnorm_mat(int iN, arma::vec vMu, arma::mat mSigma) {
  int incols = mSigma.n_cols;
  arma::mat Y = arma::randn(iN, incols);
  return arma::repmat(vMu, 1, iN).t() + Y * chol(mSigma);
}
