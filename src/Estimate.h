#ifndef ESTIMATE_H
#define ESTIMATE_H
Rcpp::List funcEstimate_Eff(arma::vec vY, arma::mat mX, arma::vec vDelta, double dAlpha, arma::vec vKeep, double dBeta, bool bZellnerPrior=true, double dG=50.0);
#endif
