#include "EstimateCore.h"
// [[Rcpp::export]]
Rcpp::List funcEstimate_Eff_par(arma::vec vY, arma::mat mX, arma::vec vDelta, double dAlpha,
                              arma::vec vKeep, double dBeta, bool bZellnerPrior = true,
                              double dG = 50.0, double iCores = 2,
                              double dS0 = 1.0, double dN0 = 2.0) {
  if (!std::isfinite(iCores) || iCores < 1 || iCores != std::floor(iCores) || iCores > 65536)
    Rcpp::stop("iCores must be a positive integer no greater than 65536.");
#ifndef _OPENMP
  Rcpp::warning("OpenMP is unavailable; using sequential execution.");
#endif
  return estimate_core(vY,mX,vDelta,dAlpha,vKeep,dBeta,bZellnerPrior,dG,static_cast<int>(iCores),dS0,dN0);
}
