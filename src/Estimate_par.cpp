#include <RcppArmadillo.h>
#include "Utils.h"
#include "Combinations.h"
#include "Estimate.h"
#if defined(_OPENMP)
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

using namespace Rcpp;
using namespace arma;


//[[Rcpp::export]]
List funcEstimate_Eff_par(arma::vec vY, arma::mat mX, arma::vec vDelta, double dAlpha, arma::vec vKeep, bool bZellnerPrior=true, double dG=50.0,int iCores=2){
  List lOut;
#if defined(_OPENMP)
  omp_set_num_threads(iCores);

  // memory allocation
  int k,i,j,d, iK = mX.n_cols, iT = vY.size(), iD = vDelta.size();

  field<uvec> fuv_Foo = PowerSet2_f(iK,vKeep);

  int iM = fuv_Foo.size();
  int iCeil = ceil(0.1 * iM);

  // DLM memory
  arma::vec  vTheta0;
  arma::mat  mVtheta0;

  field<cube> fcC(iM);
  arma::mat   mR;
  field<cube> fcmm(iM);
  arma::vec   vA;
  arma::vec   vf(iM);
  arma::cube  cS_DLM(2,iM,iD);

  arma::vec  vLpdf_DLM(iM);
  double  dmFmRmF_DLM=0.0;

  double  dn=0.0;
  double  dQ=0.0;
  double  de=0.0;

  arma::uvec vFoo;
  arma::uvec vInt(1);

  int iK_foo;

  for(j=0;j<iM;j++){
    vFoo = fuv_Foo(j);
    iK_foo     = vFoo.size();
    fcC(j)     = zeros<cube>(iK_foo,iK_foo,iD);
    fcmm(j)    = zeros<cube>(iK_foo,2,iD);
  }
  // DMA memory
  arma::mat mpm_DMA(iD,iM); mpm_DMA.fill(1.0/iM);
  //
  double dsumprob_DMA,dsumx_DMA,dsum_DMA;
  //
  // DMA-Delta memory
  arma::cube ampmt    = zeros(2,iM,iD);
  //
  arma::mat  mmm     = zeros(iD,iK);
  arma::mat  mcov    = zeros(iD,iK);
  arma::mat  mypred  = zeros(iT,iD);
  arma::mat  mLpdf   = zeros(iT,iD);
  arma::mat  mmhat   = zeros(iT,iK);
  arma::mat  mcovhat = zeros(iT,iK);
  arma::vec  vS      = zeros(iD);
  arma::vec  vFmRmF  = zeros(iD);
  arma::vec  vyhati  = zeros(iD);

  arma::rowvec vpm_DMAD(iD); vpm_DMAD.fill(1.0/iD);
  arma::rowvec vpmt_DMAD(iD);
  arma::mat mpmt_DMAD = zeros(iT,iD);mpmt_DMAD.row(0) = vpm_DMAD;

  arma::vec vpmt2_DMAD(iM);vpmt2_DMAD.fill(1.0/(iM));

  arma::mat mincpmt=zeros(iT,iK);
  arma::vec vsize(iT);vsize.zeros();

  double dsumx_DMAD,dsum_DMAD,dsumprob_DMAD,dSum_foo=0.0,dSum_foo2=0.0;
  arma::mat mFoo_DMAD;

  int iC=0;
  double idummy=0.0;

  arma::uvec vg(iM); vg.zeros();
  arma::uvec vg2_foo;
  arma::rowvec vBaz;

  //DMA-DELTA output variables
  arma::vec vyhat(iT);//
  arma::vec vLpdfhat(iT);//
  arma::vec vShat(iT);//
  arma::vec vdeltahat(iT);//

  arma::vec vobs(iT);
  arma::vec vcoeff(iT);
  arma::vec vmod(iT);
  arma::vec vtvp(iT);
  arma::vec vtotal(iT);

  vdeltahat(0) = accu(vDelta.t() % mpmt_DMAD.row(0));

  // DMS Memory and output variables
  arma::vec vyhat_DMS(iT);//
  arma::vec vLpdfhat_DMS(iT);//
  arma::vec vsize_DMS(iT);vsize_DMS.zeros();
  arma::vec vhighmp_DMS(iT);vhighmp_DMS.zeros(); // highest model probabilitites
  arma::vec vhighmpTop01_DMS(iT);vhighmpTop01_DMS.zeros(); // highest model probabilitites

  arma::vec vf_DMS(iD); // best point forecast conditionally on delta_i
  arma::vec vLpdf_DMS(iD); // best log density forecast conditionally on delta_i
  arma::mat mMax_DMS_ModelGivenDelta(iD, iT); // argmax_i P(M_i | delta_j, F_t)
  arma::vec vMax_DMS_Delta(iT); // argmax_j P(delta_j| F_t)
  arma::vec vMax_DMS_Model(iT); // argmax_j P(M_j| F_t)

  mMax_DMS_ModelGivenDelta.zeros();
  vMax_DMS_Delta.zeros();
  vMax_DMS_Model.zeros();

  // initialisation DLM
  //
  arma::mat invFoo;
  arma::mat mX_foo;
  double dS0;

  dn          = 2.0;
  for(d = 0;d<iD;d++){
    ampmt.slice(d).row(0).fill(1.0/iM);
#pragma omp parallel for default(none)                                                  \
    shared(iM,d,mpmt_DMAD,ampmt,fuv_Foo,mX,iT,vY,bZellnerPrior,dG,vf,fcmm,dn,cS_DLM,fcC)\
      private(mX_foo,j,vFoo,iK_foo,vTheta0,invFoo,dS0,mVtheta0,dQ,de,vA,mR)
      for(j = 0;j<iM;j++){
        vFoo    = fuv_Foo(j);
        mX_foo  = mX.cols(vFoo);
        iK_foo  = vFoo.size();
        vTheta0 = zeros(iK_foo);
        if(bZellnerPrior){
          invFoo   = arma::inv(mX_foo.t()*mX_foo);
          dS0      = (1.0/(iT-1.0))*as_scalar(vY.t()*(eye(iT,iT)-mX_foo*invFoo*mX_foo.t())*vY);
          mVtheta0 = dG*dS0*invFoo;
        }else{
          mVtheta0 = dG*eye(iK_foo,iK_foo);
        }
        mR            = mVtheta0;
        vf(j)         = as_scalar(mX_foo.row(0) * vTheta0);
        dQ            = as_scalar(mX_foo.row(0) * mR * mX_foo.row(0).t());
        de            = vY(0)-vf(j);
        vA            = mR*mX_foo.row(0).t() / dQ  ;
        fcmm(j).slice(d).col(0) = vTheta0 + vA  * de  ;
        cS_DLM(0,j,d)    = (pow(vY(0),2.0)+( pow(de ,2.0))/dQ )/dn;
        fcC(j).slice(d)  = mVtheta0;
      }
  }

  // size and inclusion probabilitites
#pragma omp parallel for default(none) reduction(+: dSum_foo)\
  shared(iM,fuv_Foo,vpmt2_DMAD)                              \
    private(iK_foo,vFoo)
    for(j=0;j<iM;j++){
      vFoo = fuv_Foo(j);
      iK_foo = vFoo.size();
      dSum_foo += vpmt2_DMAD(j)*iK_foo;
    }
    vsize(0) =dSum_foo;

  for(k=0;k<iK;k++){
    vg.zeros();
    iC = 0.0;
    idummy=0.0;
    //
    for(j=0;j<iM;j++){
      vFoo = fuv_Foo(j);
      idummy = accu(vFoo==k);
      if(idummy==1){
        vg(iC) = j;
        iC++ ;
      }
    }
    vg2_foo = vg.subvec(0,iC-1);
    mincpmt(0,k) += accu(vpmt2_DMAD.elem(vg2_foo));
  }

  arma::rowvec vX_t(iK);
  arma::rowvec vX_t_foo;
  //
  arma::mat mFoo1(iK,iM),mFoo2(iK,iM);
  //
  for(i=1;i<iT;i++){
    vpmt2_DMAD.zeros();
    vFmRmF.zeros();
    vX_t = mX.row(i);
    // Do DLM for each j at time t
    dn         += 1.0;
    for(d = 0;d<iD;d++){
      mFoo1.zeros();
      mFoo2.zeros();
      dSum_foo  = 0.0;
      dSum_foo2 = 0.0;
#pragma omp parallel for default(none) reduction(+:dSum_foo2,dSum_foo)                       \
      shared(fuv_Foo,iM,vX_t,vDelta,d,fcC,vf,fcmm,cS_DLM,vY,i,dn,vLpdf_DLM,ampmt,mFoo1,mFoo2)\
        private(j,vFoo,vX_t_foo,mR,dQ,de,vA,dmFmRmF_DLM,iK_foo,k)
        for(j=0;j<iM;j++){
          vFoo = fuv_Foo(j);
          iK_foo = vFoo.size();
          vX_t_foo                  = vX_t.elem(vFoo).t();
          mR                        = (1.0/vDelta(d))*fcC(j).slice(d) ;
          vf(j)                     = as_scalar(vX_t_foo * fcmm(j).slice(d).col(0));
          dQ                        = as_scalar(vX_t_foo * mR * vX_t_foo.t()) + cS_DLM(0,j,d);
          de                        = vY(i) - vf(j) ;
          vA                        = mR*vX_t_foo.t()/dQ  ;
          fcmm(j).slice(d).col(1)   = fcmm(j).slice(d).col(0) + vA * de ;
          cS_DLM(1,j,d)             = cS_DLM(0,j,d)+(cS_DLM(0,j,d)/ (dn-1) )*(( pow(de ,2.0)/dQ  ) - 1.0);
          fcC(j).slice(d)           = cS_DLM(1,j,d)*( mR-vA*vA.t()*dQ  )/cS_DLM(0,j,d);
          vLpdf_DLM(j)              = Rf_dt(de /pow(dQ  ,0.5),dn,1) - 0.5*log(dQ );
          dmFmRmF_DLM               = as_scalar(vX_t_foo*mR*vX_t_foo.t());
          //
          dSum_foo              += ampmt(0,j,d) * dmFmRmF_DLM;//
          dSum_foo2             += ampmt(0,j,d) * exp(vLpdf_DLM(j));
          //
          for(k=0;k<iK_foo;k++){
            mFoo1(vFoo(k),j) = fcmm(j).at(k,1,d);
            mFoo2(vFoo(k),j) = fcC(j).at(k,k,d);
          }
          //store
          fcmm(j).slice(d).col(0) = fcmm(j).slice(d).col(1);
        }

        vFmRmF(d) += dSum_foo;
      mLpdf(i,d) = log(dSum_foo2);

      // Do DMA at time t
      dsumx_DMA             = accu(pow(mpm_DMA.row(d),dAlpha));
      mpm_DMA.row(d)        = pow(mpm_DMA.row(d),dAlpha)/dsumx_DMA;
      dsum_DMA              = accu(exp(vLpdf_DLM.t()) % mpm_DMA.row(d));
      ampmt.slice(d).row(1) = (exp(vLpdf_DLM.t()) % mpm_DMA.row(d))/dsum_DMA;
      dsumprob_DMA          = accu(ampmt.slice(d).row(1));
      mpm_DMA.row(d)        = ampmt.slice(d).row(1)/dsumprob_DMA;
      //
      mmm.row(d)            = arma::trans(sum(mFoo1 % repmat(ampmt.slice(d).row(0),iK,1),1));
      mcov.row(d)           = arma::trans(sum(mFoo2 % repmat(ampmt.slice(d).row(0),iK,1),1));
      mypred(i,d)           = accu(ampmt.slice(d).row(0) % vf.t());
      vS(d)                 = accu(ampmt.slice(d).row(0) % cS_DLM.slice(d).row(1));
      //
      vyhati(d)             = accu(  pow(vf-mypred(i,d),2.0) % ampmt.slice(d).row(1).t());

      // Do DMS conditionally on delta
      vf_DMS(d)    = vf(mMax_DMS_ModelGivenDelta(d, i-1));
      vLpdf_DMS(d) = vLpdf_DLM(mMax_DMS_ModelGivenDelta(d, i-1));
      mMax_DMS_ModelGivenDelta(d, i)  = MaxFinder(arma::trans(mpm_DMA.row(d)));

      //store
      cS_DLM.slice(d).row(0) = cS_DLM.slice(d).row(1);
      ampmt.slice(d).row(0)  = ampmt.slice(d).row(1);
    }
    // Do DMA-D at time t
    dsumx_DMAD       = accu(pow(vpm_DMAD,dAlpha));
    vpm_DMAD         = pow(vpm_DMAD,dAlpha)/dsumx_DMAD;
    dsum_DMAD        = accu(exp(mLpdf.row(i)) % vpm_DMAD);
    vpmt_DMAD        = (exp(mLpdf.row(i)) % vpm_DMAD)/dsum_DMAD;
    mpmt_DMAD.row(i) = vpmt_DMAD;
    //
    mmhat.row(i)   = arma::trans(sum(mmm.t() % repmat(mpmt_DMAD.row(i-1),iK,1),1));
    mcovhat.row(i) = arma::trans(sum(mcov.t() % repmat(mpmt_DMAD.row(i-1),iK,1),1));
    //
    dsumprob_DMAD = accu(vpmt_DMAD);
    vpm_DMAD      = vpmt_DMAD/dsumprob_DMAD;

    for(d = 0;d<iD;d++){
      for(j=0;j<iM;j++){
        vpmt2_DMAD(j) += ampmt(1,j,d) * mpmt_DMAD(i,d);
      }
    }
    // size and inclusion probabilities
    dSum_foo = 0.0;
#pragma omp parallel for default(none) reduction(+: dSum_foo)\
    shared(iM,fuv_Foo,vpmt2_DMAD)                            \
      private(iK_foo,vFoo)
      for(j=0;j<iM;j++){
        vFoo = fuv_Foo(j);
        dSum_foo += vpmt2_DMAD(j)*vFoo.size();
      }
      vsize(i) = dSum_foo;

    for(k=0;k<iK;k++){
      vg.zeros();
      iC = 0.0;
      idummy=0.0;
      //
      for(j=0;j<iM;j++){
        vFoo = fuv_Foo(j);
        idummy = accu(vFoo==k);
        if(idummy==1){
          vg(iC) = j;
          iC++ ;
        }
      }
      vg2_foo = vg.subvec(0,iC-1);
      mincpmt(i,k) += accu(vpmt2_DMAD.elem(vg2_foo));
    }

    // Do DMS at time t

    // vyhat_DMS(i)    = accu(vf_DMS % vpm_DMAD.t());
    // vLpdfhat_DMS(i) = log(accu(exp(vLpdf_DMS) % vpm_DMAD.t()));
    vyhat_DMS(i)    = vf_DMS(vMax_DMS_Delta(i - 1));
    vLpdfhat_DMS(i) = vLpdf_DMS(vMax_DMS_Delta(i - 1));
    vMax_DMS_Delta(i) = MaxFinder(vpm_DMAD.t());

    // find max over P(M_i|F_t)
    vMax_DMS_Model(i) = MaxFinder(vpmt2_DMAD);
    vhighmp_DMS(i)  = vpmt2_DMAD(vMax_DMS_Model(i));
    vpmt2_DMAD = sort(vpmt2_DMAD, "descend");
    vhighmpTop01_DMS(i) =  accu(vpmt2_DMAD.subvec(0, iCeil - 1));
    vFoo            = fuv_Foo(vMax_DMS_Model(i));
    vsize_DMS(i)    = vFoo.size();

    // store output quantities
    vyhat(i)     = accu(mypred.row(i) % mpmt_DMAD.row(i));
    vLpdfhat(i)  = log(accu(exp(mLpdf.row(i)) % mpmt_DMAD.row(i)));
    vShat(i)     = accu(vS.t() % mpmt_DMAD.row(i));
    vdeltahat(i) = accu(vDelta.t() % mpmt_DMAD.row(i));

    vobs(i)      = accu(vS.t() % mpmt_DMAD.row(i));
    vcoeff(i)    = accu(vFmRmF.t() % mpmt_DMAD.row(i));
    vmod(i)      = accu(vyhati.t() % mpmt_DMAD.row(i));
    vtvp(i)      = accu(pow(mypred.row(i) - vyhat(i),2.0));
    vtotal(i)    = vobs(i) + vcoeff(i) + vmod(i) + vtvp(i);

  }

  lOut["mincpmt"] = mincpmt;
  lOut["vsize"] = vsize;
  lOut["mmhat"] = mmhat;
  lOut["mcovhat"] = mcovhat;
  lOut["mpmt"] = mpmt_DMAD;
  lOut["vyhat"] = vyhat;
  lOut["veps"] = vY - vyhat;
  lOut["vLpdfhat"] = vLpdfhat;
  lOut["mLpdf"] = mLpdf;
  lOut["vShat"] = vShat;
  lOut["vdeltahat"] = vdeltahat;
  // variance
  arma::mat mvdec(iT, 5);

  mvdec.col(0) = vtotal;
  mvdec.col(1) = vobs;
  mvdec.col(2) = vcoeff;
  mvdec.col(3) = vmod;
  mvdec.col(4) = vtvp;

  lOut["mvdec"] = mvdec;
  lOut["vobs"] = vobs;
  lOut["vcoeff"] = vcoeff;
  lOut["vmod"] = vmod;
  lOut["vtvp"] = vtvp;
  lOut["vtotal"] = vtotal;
  lOut["iM"]  = iM;
  // DMS output
  lOut["vsize_DMS"]    = vsize_DMS;
  lOut["vyhat_DMS"]    = vyhat_DMS;
  lOut["vLpdfhat_DMS"] = vLpdfhat_DMS;
  lOut["vhighmp_DMS"]  = vhighmp_DMS;
  lOut["vhighmpTop01_DMS"]  = vhighmp_DMS;
  lOut["veps_DMS"]     = vY - vyhat_DMS;
  //
#else
  Rf_warning("OpenMP is not available. Sequential processing is used.");
  lOut = funcEstimate_Eff(vY, mX, vDelta, dAlpha, vKeep, bZellnerPrior, dG);
#endif
  return lOut;
}
