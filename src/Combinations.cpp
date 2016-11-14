#include <RcppArmadillo.h>
#include "Utils.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List PowerSet(int iK) {

  int i,j;
  NumericVector els(iK);
  for(i=0;i<iK;i++)
    els[i]=i;

  int pwrset_card = int_pow(2,iK);
  List out(pwrset_card);
  out[0] = NumericVector::create();
  NumericVector tmp;
  int counter = 0;
  for (i=0; i < iK; ++i) {
    int cnt2 = counter;            // capture counter state
    for (j =0; j <= counter; ++j) {
      cnt2++;                   // capture counter j steps
      tmp = as<NumericVector>(out[j]);
      tmp.push_back(els[i]);
      out[cnt2] = Rcpp::as<arma::uvec>(wrap(tmp));
    }
    counter = cnt2;             // update counter state
  }

  return out;
}

field<uvec> PowerSet2(int iK) {

  int i,j;
  NumericVector els(iK);
  for(i=0;i<iK;i++)
    els[i]=i;

  int pwrset_card = int_pow(2,iK);
  NumericVector tmp;
  field<uvec> fOut(pwrset_card);

  int counter = 0;
  for (i=0; i < iK; ++i) {
    int cnt2 = counter;            // capture counter state
    for (j =0; j <= counter; ++j) {
      cnt2++;                   // capture counter + j steps
      tmp = as<NumericVector>(wrap(fOut(j)));
      tmp.push_back(els[i]);
      fOut(cnt2) = Rcpp::as<arma::uvec>(wrap(tmp));
    }
    counter = cnt2;             // update counter state
  }

  fOut = fOut.rows(1,pwrset_card-1);

  return fOut;
}

field<uvec> PowerSet2_withkeep(int iK, arma::vec vKeep) {

  int i,j,k;
  NumericVector els(iK);
  for(i=0;i<iK;i++)
    els[i]=i;

  int pwrset_card = int_pow(2,iK);
  NumericVector tmp;
  field<uvec> fOut(pwrset_card);

  int vKeep_size = vKeep.size();
  double idummy=0.0;

  int cnt, counter = 0;
  for (i=0; i < iK; ++i) {
    cnt = counter;            // capture counter state
    for (j =0; j <= counter; ++j) {
      tmp = as<NumericVector>(wrap(fOut(j)));
      tmp.push_back(els[i]);
      cnt++;                   // capture counter + j steps
      fOut(cnt) = Rcpp::as<arma::uvec>(wrap(tmp));
    }
    counter = cnt;             // update counter state
  }

  arma::vec vKeep_Index(pwrset_card);
  arma::uvec tmp2;
  int Size2Keep = 0;

  for(j=1;j<pwrset_card;j++){
    idummy = 0;
    tmp2   = fOut(j);
    for(k=0;k<vKeep_size;k++){
      idummy += sum(tmp2==vKeep(k));
    }
    if(idummy==vKeep_size){
      vKeep_Index(Size2Keep) = j;
      Size2Keep++;
    }
  }

  field<uvec> fOut2(Size2Keep);

  for(j=0;j<Size2Keep;j++){
    fOut2(j) = fOut(vKeep_Index(j));
  }

  return fOut2;
}

field<uvec> PowerSet2_f(int iK, arma::vec vKeep){
  field<uvec> fOut;

  if(vKeep(0) == -9999){
    fOut = PowerSet2(iK);
  }else{
    fOut = PowerSet2_withkeep(iK, vKeep);
  }

  return fOut;
}

