#include "EstimateCore.h"
// [[Rcpp::export]]
Rcpp::List funcEstimate_Eff(arma::vec vY, arma::mat mX, arma::vec vDelta, double dAlpha,
                          arma::vec vKeep, double dBeta, bool bZellnerPrior = true,
                          double dG = 50.0, double dS0 = 1.0, double dN0 = 2.0) {
  return estimate_core(vY,mX,vDelta,dAlpha,vKeep,dBeta,bZellnerPrior,dG,1,dS0,dN0);
}
