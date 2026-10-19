#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include "Combinations.h"

// Enumerate optional bits only, in the original ascending bit-mask order.
arma::field<arma::uvec> PowerSet2_f(int K, arma::vec keep) {
  if (K < 1) Rcpp::stop("The design must contain at least one column.");
  std::vector<bool> required(K, false);
  bool none = keep.n_elem == 1 && keep(0) == -9999;
  if (!none) {
    if (keep.is_empty()) Rcpp::stop("vKeep must be nonempty or NULL.");
    for (double k : keep) {
      if (!std::isfinite(k) || k != std::floor(k) || k < 0 || k >= K || required[k])
        Rcpp::stop("Invalid or duplicated vKeep index.");
      required[k] = true;
    }
  }
  std::vector<unsigned> optional;
  for (int k=0; k<K; ++k) if (!required[k]) optional.push_back(k);
  if (optional.size() > 20) Rcpp::stop("At most 20 optional design columns are supported; use vKeep to reduce the model space.");
  size_t count = size_t(1) << optional.size();
  arma::field<arma::uvec> out(count - (none ? 1 : 0));
  size_t at=0;
  for (size_t bits=(none ? 1 : 0); bits<count; ++bits) {
    std::vector<arma::uword> columns;
    size_t b=0;
    for (int k=0; k<K; ++k) {
      if (required[k]) columns.push_back(k);
      else { if (bits & (size_t(1)<<b)) columns.push_back(k); ++b; }
    }
    out(at++) = arma::uvec(columns);
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List PowerSet(double iK) {
  if (!std::isfinite(iK) || iK != std::floor(iK) || iK < 0 || iK > 20)
    Rcpp::stop("iK must be an integer between 0 and 20.");
  size_t count = size_t(1) << static_cast<unsigned>(iK);
  Rcpp::List out(count);
  for (size_t bits=0; bits<count; ++bits) {
    std::vector<double> columns;
    for (unsigned k=0; k<iK; ++k) if (bits & (size_t(1)<<k)) columns.push_back(k);
    out[bits] = Rcpp::wrap(columns);
  }
  return out;
}
