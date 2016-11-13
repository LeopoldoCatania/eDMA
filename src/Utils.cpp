#include <RcppArmadillo.h>
#if defined(_OPENMP)
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]
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
