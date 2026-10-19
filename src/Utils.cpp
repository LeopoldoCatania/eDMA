#include <RcppArmadillo.h>
#include <limits>

using namespace Rcpp;
using namespace arma;

int MaxFinder(arma::vec vX){

  int iN = vX.size();
  int i;
  int iMax = 0;
  double dX_max = vX(0);

  for(i=1;i<iN;i++){
    if(vX(i)>dX_max){
      dX_max = vX(i);
      iMax = i;
    }
  }
  return iMax;
}

int int_pow(int base, int exp) {
  if (base < 0 || exp < 0) Rcpp::stop("Integer power requires nonnegative arguments.");
  int result = 1;
  for (int i=0; i<exp; ++i) {
    if (base && result > std::numeric_limits<int>::max()/base)
      Rcpp::stop("Integer power overflow.");
    result *= base;
  }
  return result;
}
