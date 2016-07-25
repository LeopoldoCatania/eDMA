#ifndef COMBINATIONS_H
#define COMBINATIONS_H
Rcpp::List PowerSet(int iK);
arma::field<arma::uvec> PowerSet2(int iK);
arma::field<arma::uvec> PowerSet2_f(int iK, arma::vec vKeep);
#endif
