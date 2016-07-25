#include <RcppArmadillo.h>
#include "Distributions.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
List SimulateDLM(int iT, arma::mat mX, arma::vec vBeta0, arma::mat mW,
                 double dV, double dPhi){

  int i;
  int iK = mX.n_cols;

  arma::mat mBeta(iK,iT);
  arma::vec vY(iT);

  mBeta.col(0) = vBeta0;
  vY(0) = as_scalar(mX.row(0) * mBeta.col(0)) + pow(dV,0.5) * Rf_rnorm(0.0,1.0);

  arma::vec vMu_foo(iK);vMu_foo.zeros();

  for(i = 1;i<iT;i++){
    mBeta.col(i) = dPhi*mBeta.col(i-1) + arma::trans(rmvnorm_mat(1, vMu_foo, mW));
    vY(i)        = as_scalar(mX.row(i) * mBeta.col(i)) + pow(dV,0.5) * Rf_rnorm(0.0,1.0);
  }

  List lOut;

  lOut["vY"]    = vY;
  lOut["mBeta"] = mBeta;

  return lOut;

}

