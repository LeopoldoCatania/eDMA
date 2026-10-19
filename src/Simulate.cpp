#include <RcppArmadillo.h>
#include <cmath>
// [[Rcpp::export]]
Rcpp::List SimulateDLM(double iT, arma::mat mX, arma::vec vBeta0, arma::mat mW,
                      double dV, double dPhi) {
  if (!std::isfinite(iT) || iT<1 || iT!=std::floor(iT) || iT!=mX.n_rows || !mX.n_cols ||
      vBeta0.n_elem!=mX.n_cols || mW.n_rows!=mX.n_cols || mW.n_cols!=mX.n_cols)
    Rcpp::stop("iT must be a positive integer matching the design rows, with compatible coefficient/covariance dimensions.");
  if (!mX.is_finite() || !vBeta0.is_finite() || !mW.is_finite() || !std::isfinite(dV) || dV<0 || !std::isfinite(dPhi))
    Rcpp::stop("Simulation inputs must be finite and dV must be nonnegative.");
  double tol=1e-12*std::max(1.0,arma::norm(mW,"inf"));
  if (arma::norm(mW-mW.t(),"inf")>tol) Rcpp::stop("mW must be symmetric.");
  arma::vec eigen; arma::mat vectors;
  if (!arma::eig_sym(eigen,vectors,0.5*(mW+mW.t())) || eigen.min() < -tol)
    Rcpp::stop("mW must be positive semidefinite.");
  arma::mat root=vectors*arma::diagmat(arma::sqrt(arma::clamp(eigen,0.0,arma::datum::inf)));
  arma::uword T=mX.n_rows,K=mX.n_cols;
  arma::mat coefficients(K,T); arma::vec y(T);
  coefficients.col(0)=vBeta0;
  for (arma::uword t=0;t<T;++t) {
    if (t) coefficients.col(t)=dPhi*coefficients.col(t-1)+root*arma::randn<arma::vec>(K);
    y(t)=arma::dot(mX.row(t),coefficients.col(t))+std::sqrt(dV)*R::rnorm(0,1);
  }
  return Rcpp::List::create(Rcpp::Named("vY")=y,Rcpp::Named("mBeta")=coefficients);
}
