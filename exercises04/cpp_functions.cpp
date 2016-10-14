#include <RcppArmadillo.h>
//#include <cmath.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

//[[Rcpp::export]]
List AdaGrad_biggish_cpp(vector<unsigned int> y,
                         vector<unsigned int> X_dim,
                         vector<unsigned int> X_p,
                         vector<unsigned int> X_j,
                         vector<double> X_x,
                         double alpha0,
                         vector<double> beta0,
                         double stepsize,
                         double ada_eps){
  //Rcpp::Rcout << "alpha is " << alpha0 << std::endl;
  int n=X_dim[0];
  int p=X_dim[1];
  vector<double> unit_negloglik(n);
  double diag_G_const=ada_eps;
  vector<double> diag_G(p);
  for(int k=0;k<p;++k){diag_G[k]=ada_eps;}
  for(int i=0;i<n;++i){
    int j_start=X_p[i];
    int j_end=X_p[i+1]-1;
    int looping=j_end-j_start+1;
    vector<double> active_Xs(looping);
    vector<double> values_Xs(looping); 
    double Xtbeta=0;
    for(int k=0;k<looping;++k){
      active_Xs[k]=X_j[j_start+k];
      values_Xs[k]=X_x[j_start+k];
      Xtbeta+=values_Xs[k]*beta0[active_Xs[k]];}
    Xtbeta+=alpha0;
    double value_y=y[i];
    double expbeta=1+exp(-Xtbeta);
    unit_negloglik[i]=(1-value_y)*Xtbeta+log(expbeta);
    double grad_const=1/expbeta-value_y;
    diag_G_const+=grad_const*grad_const;
    alpha0-=stepsize/sqrt(diag_G_const)*grad_const;
    vector<double> gradient(looping);
    for(int k=0;k<looping;++k){
      gradient[k]=values_Xs[k]*grad_const;
      diag_G[active_Xs[k]]+=gradient[k]*gradient[k];
      beta0[active_Xs[k]]-=stepsize/sqrt(diag_G[active_Xs[k]])*gradient[k];}}
  return Rcpp::List::create(
    _["alphahat"]=alpha0,
    _["betahat"]=beta0,
    _["unit_negloglik"]=unit_negloglik);
}
