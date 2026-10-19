#include "EstimateCore.h"
#include "Combinations.h"
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif
using arma::mat;
using arma::vec;
using arma::uword;

namespace {
uword first_max(const vec& z) {
  uword best=0;
  for (uword i=1;i<z.n_elem;++i) if (z(i)>z(best)) best=i;
  return best;
}
double lse(const vec& z) {
  double m=z.max();
  if (!std::isfinite(m)) Rcpp::stop("Invalid probability normalizer.");
  return m + std::log(arma::accu(arma::exp(z-m)));
}
// No R API calls: safe in OpenMP workers. Avoid squaring a large standardized error.
double log_student(double e, double q, double nu) {
  double z = e == 0 ? -std::numeric_limits<double>::infinity() :
    2*std::log(std::abs(e))-std::log(q)-std::log(nu);
  double tail = z > 0 ? z+std::log1p(std::exp(-z)) : std::log1p(std::exp(z));
  return std::lgamma((nu+1)/2)-std::lgamma(nu/2)-0.5*(std::log(nu)+std::log(arma::datum::pi)+std::log(q))-(nu+1)*tail/2;
}
struct State { vec m; mat C; double S; };
}

Rcpp::List estimate_core(const vec& y, const mat& X, const vec& delta, double alpha,
                        const vec& keep, double beta, bool zellner, double g, int cores,
                        double S0, double n0) {
  const uword T=y.n_elem, K=X.n_cols, D=delta.n_elem;
  if (!T || X.n_rows!=T || !K || !X.is_finite()) Rcpp::stop("Invalid response/design dimensions or nonfinite design.");
  bool future=std::isnan(y(T-1));
  uword N=T-(future ? 1 : 0);
  if (!N || !y.head(N).is_finite()) Rcpp::stop("Responses must be finite, except for an optional terminal NA.");
  if (!D || !delta.is_finite() || arma::any(delta<=0) || arma::any(delta>1) ||
      !std::isfinite(alpha) || alpha<=0 || alpha>1 || !std::isfinite(beta) || beta<=0 || beta>1 ||
      !std::isfinite(g) || g<=0 || !std::isfinite(S0) || S0<=0 || !std::isfinite(n0) || n0<=0)
    Rcpp::stop("Discounts must be in (0,1]; dG, dS0 and dN0 must be finite and positive.");
  arma::field<arma::uvec> models=PowerSet2_f(K,keep);
  const uword M=models.n_elem;
  // Bound resident allocations before creating model/delta covariance arrays.
  long double cells=static_cast<long double>(T)*(3*K+3*D+24)+static_cast<long double>(M)*D*10;
  for (uword j=0;j<M;++j) { long double p=models(j).n_elem; cells+=D*(p*p+p)+p; }
  if (cells*sizeof(double)>1024.0L*1024*1024) Rcpp::stop("Estimated working memory exceeds 1 GiB; reduce optional predictors or the delta grid.");
  std::vector<State> states(M*D);
  for (uword j=0;j<M;++j) {
    uword p=models(j).n_elem;
    mat C=g*arma::eye(p,p);
    if (zellner) {
      if (N<2 || N<p) Rcpp::stop("Zellner prior requires at least two observations and a full-column-rank design.");
      mat selected=X.cols(models(j));
      mat Z=selected.head_rows(N);
      mat Q,R;
      if (!arma::qr_econ(Q,R,Z) || arma::rcond(R)<1e-12)
        Rcpp::stop("Rank-deficient or ill-conditioned design in Zellner prior (model %d).",static_cast<int>(j+1));
      vec residual=y.head(N)-Q*(Q.t()*y.head(N));
      double s=std::max(arma::dot(residual,residual)/(N-1),S0*1e-12);
      mat inverseR=arma::solve(arma::trimatu(R),arma::eye(p,p));
      C=g*s*(inverseR*inverseR.t());
      if (!C.is_finite()) Rcpp::stop("Nonfinite Zellner covariance.");
    }
    for (uword d=0;d<D;++d) states[d*M+j]={vec(p,arma::fill::zeros),C,S0};
  }
  mat loga(M,D); loga.fill(-std::log(double(M)));
  vec logb(D); logb.fill(-std::log(double(D)));
  mat inc(T,K,arma::fill::zeros), mean(T,K,arma::fill::zeros), cov(T,K,arma::fill::zeros);
  mat postdelta(T,D), preddelta(T,D), loglikdelta(T,D), vdec(T,5);
  vec forecast(T), errors(T), loglik(T), Shat(T), deltahat(T), size(T), df(T);
  vec dmsforecast(T), dmsloglik(T), dmserrors(T), dmssize(T), maxprob(T), topmass(T);
  mat f(M,D), obs(M,D), coeff(M,D), lp(M,D);
  double n=n0;
  for (uword t=0;t<T;++t) {
    Rcpp::checkUserInterrupt();
    const bool observed=t<N;
    const double nu=beta*n+1;
    df(t)=nu;
    mat la=alpha*loga;
    for (uword d=0;d<D;++d) la.col(d)-=lse(la.col(d));
    vec lb=alpha*logb; lb-=lse(lb);
    mat a=arma::exp(la); vec b=arma::exp(lb);
    preddelta.row(t)=b.t();
    uword dstar=first_max(lb), jstar=first_max(la.col(dstar));
    const vec fullx=X.row(t).t();
    std::vector<std::string> failures(M*D);
    const int workers=static_cast<int>(std::min<uword>(cores,M*D));
#ifndef _OPENMP
    (void)workers;
#endif
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(workers) if(workers>1)
#endif
    for (long long h=0;h<static_cast<long long>(M*D);++h) {
      try {
        uword d=h/M,j=h%M;
        State& s=states[h];
        vec x=fullx.elem(models(j));
        mat R=s.C/(t==0 ? 1.0 : delta(d));
        vec Rx=R*x;
        double c=arma::dot(x,Rx), mu=arma::dot(x,s.m), q=c+s.S;
        if (!std::isfinite(q) || q<=0 || !std::isfinite(mu) || c < -1e-10*std::max(1.0,q))
          throw std::runtime_error("invalid predictive mean or covariance");
        f(j,d)=mu; obs(j,d)=s.S; coeff(j,d)=std::max(0.0,c);
        lp(j,d)=NA_REAL;
        if (observed) {
          double e=y(t)-mu;
          lp(j,d)=log_student(e,q,nu);
          vec gain=Rx/q;
          vec updated=s.m+gain*e;
          // Rank-one plug-in covariance update, O(p^2) per component.
          mat C=R-gain*Rx.t();
          double standardized=e/std::sqrt(q);
          double S=s.S*((nu-1)/nu)+(s.S/nu*standardized)*standardized;
          if (!std::isfinite(lp(j,d)) || !updated.is_finite() || !C.is_finite() || !std::isfinite(S) || S<=0)
            throw std::runtime_error("nonfinite state update; rescale the data or increase the prior scale");
          s.m=updated; s.C=arma::symmatu(C); s.S=S;
        }
      } catch (const std::exception& e) { failures[h]=e.what(); }
      catch (...) { failures[h]="unknown numerical failure"; }
    }
    for (uword h=0;h<M*D;++h) if (!failures[h].empty())
      Rcpp::stop("Time %d, model %d, delta %d: %s",int(t+1),int(h%M+1),int(h/M+1),failures[h].c_str());
    vec mud=arma::sum(a%f,0).t();
    double mu=arma::dot(b,mud), vo=0,vc=0,vm=0,vt=0;
    for (uword d=0;d<D;++d) {
      vo+=b(d)*arma::dot(a.col(d),obs.col(d));
      vc+=b(d)*arma::dot(a.col(d),coeff.col(d));
      vm+=b(d)*arma::dot(a.col(d),arma::square(f.col(d)-mud(d)));
      vt+=b(d)*std::pow(mud(d)-mu,2);
    }
    // Student-t scale is not its variance. Zero coefficient variance remains zero.
    if (nu>2) { vo*=nu/(nu-2); vc*=nu/(nu-2); }
    else { vo=arma::datum::inf; if (vc>0) vc=arma::datum::inf; }
    vdec.row(t)=arma::rowvec({vo+vc+vm+vt,vo,vc,vm,vt});
    forecast(t)=mu; errors(t)=observed ? y(t)-mu : NA_REAL;
    dmsforecast(t)=f(jstar,dstar); dmsloglik(t)=lp(jstar,dstar);
    dmserrors(t)=observed ? y(t)-dmsforecast(t) : NA_REAL;
    if (observed) {
      vec ld(D);
      for (uword d=0;d<D;++d) {
        ld(d)=lse(la.col(d)+lp.col(d));
        loga.col(d)=la.col(d)+lp.col(d)-ld(d);
      }
      loglik(t)=lse(lb+ld); logb=lb+ld-loglik(t);
      loglikdelta.row(t)=ld.t(); n=nu;
    } else {
      loglik(t)=NA_REAL; loglikdelta.row(t).fill(NA_REAL);
      // Prediction only: no response, state, degrees-of-freedom or posterior-weight update.
    }
    mat ap=arma::exp(loga); vec bp=arma::exp(logb), marginal=ap*bp;
    postdelta.row(t)=bp.t(); deltahat(t)=arma::dot(bp,delta); Shat(t)=0; size(t)=0;
    for (uword d=0;d<D;++d) for (uword j=0;j<M;++j) {
      double w=bp(d)*ap(j,d); const State& s=states[d*M+j];
      Shat(t)+=w*s.S;
      for (uword k=0;k<models(j).n_elem;++k) {
        uword col=models(j)(k);
        mean(t,col)+=w*s.m(k);
        inc(t,col)+=w;
      }
    }
    // Mixture variance includes the mass at zero for absent coefficients.
    for (uword d=0;d<D;++d) for (uword j=0;j<M;++j) {
      double w=bp(d)*ap(j,d); const State& s=states[d*M+j];
      vec padded(K,arma::fill::zeros); padded.elem(models(j))=s.m;
      cov.row(t)+=w*arma::square(padded.t()-mean.row(t));
      for (uword k=0;k<models(j).n_elem;++k) cov(t,models(j)(k))+=w*s.C(k,k);
    }
    for (uword j=0;j<M;++j) size(t)+=marginal(j)*models(j).n_elem;
    uword best=first_max(marginal); dmssize(t)=models(best).n_elem; maxprob(t)=marginal(best);
    vec sorted=arma::sort(marginal,"descend");
    topmass(t)=arma::accu(sorted.head(static_cast<uword>(std::ceil(0.1*M))));
  }
  Rcpp::List out;
  out["mincpmt"]=inc; out["vsize"]=size; out["mmhat"]=mean; out["mcovhat"]=cov;
  out["mpmt"]=postdelta; out["mpmt_pred"]=preddelta; out["vdf"]=df;
  out["vyhat"]=forecast; out["veps"]=errors; out["vLpdfhat"]=loglik; out["mLpdf"]=loglikdelta;
  out["vShat"]=Shat; out["vdeltahat"]=deltahat; out["mvdec"]=vdec;
  out["vtotal"]=vec(vdec.col(0)); out["vobs"]=vec(vdec.col(1)); out["vcoeff"]=vec(vdec.col(2));
  out["vmod"]=vec(vdec.col(3)); out["vtvp"]=vec(vdec.col(4));
  out["vsize_DMS"]=dmssize; out["vyhat_DMS"]=dmsforecast; out["vLpdfhat_DMS"]=dmsloglik;
  out["vhighmp_DMS"]=maxprob; out["vhighmpTop01_DMS"]=topmass; out["veps_DMS"]=dmserrors;
  out["iM"]=static_cast<double>(M);
  return out;
}
