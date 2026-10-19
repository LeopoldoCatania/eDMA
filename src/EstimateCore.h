#ifndef EDMA_ESTIMATE_CORE_H
#define EDMA_ESTIMATE_CORE_H
#include <RcppArmadillo.h>
Rcpp::List estimate_core(const arma::vec&, const arma::mat&, const arma::vec&, double,
                        const arma::vec&, double, bool, double, int, double, double);
#endif
