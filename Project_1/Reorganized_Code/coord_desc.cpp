#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double myabs(double x){
  if(x < 0){
    x = -1*x;
  }
  return x;
}

// [[Rcpp::export]]
double lasso(double x, double y){
  double result = 0;
  if(myabs(x) <= y){
    result = 0;
  } else{
    result = x - y*arma::sign(x);
  }
  return result;
}

// [[Rcpp::export]]

arma::rowvec subcolvec(arma::mat A, int a){
  arma::rowvec c(A.n_rows-1);
  for(int k = 0; k < A.n_rows; ++k){
    if((k != a) & (k < a)){
      c(k) = A(k,a);
    }
    if((k != a) & (k > a)){
      c(k-1) = A(k,a);
    }
  }
  return c;
}

// [[Rcpp::export]]

arma::vec subrowvec(arma::mat A, int a, int b){
  arma::vec c(A.n_cols-1);
  for(int k = 0; k < A.n_cols; ++k){
    if((k != b) & (k < b)){
      c(k) = A(a,k);
    }
    if((k != b) & (k > b)){
      c(k-1) = A(a,k);
    }
  }
  return c;
}

// [[Rcpp::export]]
arma::mat coordesc(int n, arma::mat AB, arma::vec AB0, arma::mat XY, arma::mat EW, arma::mat condvar, arma::mat Phi12, double tuningp, arma::mat weights) {
  int a = AB.n_rows, b = AB.n_cols, iter = 0;
  double delta = 0, theta = 0, omega = 0, lambdaAB = 0, abbar = 0, phi12val = 0;
  arma::vec inter (n), XYrow (n), AB0row (n);
  arma::rowvec EWrow (n);
  arma::mat secmom = EW*EW.t() + n*condvar;
  for(int i = 0; i < a; ++i){
    XYrow = XY.col(i);
    AB0row = arma::vec(n, arma::fill::value(AB0(i)));
    inter = XYrow-AB0row;
    phi12val = Phi12(i,i);
    for(int j = 0; j < b; ++j){
      if(i > j){
        EWrow = EW.row(j);
        delta = secmom(j,j);
        theta = arma::dot(EWrow,inter);
        omega = arma::dot(subcolvec(secmom,j),subrowvec(AB, i, j));
        abbar = ((theta-omega))/(delta);
        if(weights(i,j) == 0){
          lambdaAB = (tuningp*phi12val)/(delta*(1/pow(1,-5)));
        } 
        else {
          lambdaAB = tuningp*phi12val/(delta*myabs(weights(i,j)));
        }
        AB(i,j) = lasso(abbar, lambdaAB);
      }
    }
  }
  return AB;
}

