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
  int d = A.n_rows-1;
  for(int k = 0; k < d; ++k){
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
  int d = A.n_cols-1;
  for(int k = 0; k < d; ++k){
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
arma::mat coordesc_constr(int n, arma::mat loadest, arma::mat loadfix, arma::mat XY, arma::mat Eest, arma::mat Efix, arma::mat condvarest, arma::mat condvarcross, arma::mat Phi12, double tuningp, arma::mat weights) {
  int a = loadest.n_rows, b = loadest.n_cols, c = loadfix.n_cols;
  double epsilon = 0, theta = 0, omega = 0, tau = 0, lambdaAB = 0, abbar = 0, phi12val = 0;
  arma::colvec XYrow (n);
  arma::rowvec Eestrow (n);
  arma::mat secmomest = Eest*Eest.t() + n*condvarest, secmomcross = Eest*Efix.t() + n*condvarcross;
  for(int i = 0; i < a; ++i){
    XYrow = XY.col(i);
    phi12val = Phi12(i,i);
    for(int j = 0; j < b; ++j){
      if(i > j){
        Eestrow = Eest.row(j);
        epsilon = secmomest(j,j);
        theta = arma::dot(Eestrow,XYrow);
        omega = arma::dot(subcolvec(secmomest,j),subrowvec(loadest, i, j));
        for(int k = 0; k < c; k++){
          tau += secmomcross(k,j)*loadfix(i,k);
        }
        abbar = (theta-omega-tau)/(epsilon);
        if(weights(i,j) == 0){
          lambdaAB = (tuningp*phi12val)/(epsilon*(1/pow(1,-5)));
        } 
        else {
          lambdaAB = tuningp*phi12val/(epsilon*myabs(weights(i,j)));
        }
        loadest(i,j) = lasso(abbar, lambdaAB);
      }
    }
  }
  return loadest;
}

// [[Rcpp::export]]
arma::mat coordesc_full(int n, arma::mat loadest, arma::mat loadfix, arma::mat XY, arma::mat Eest, arma::mat Efix, arma::mat condvarest, arma::mat condvarcross, arma::mat Phi12, double tuningp, arma::mat weights) {
  int a = loadest.n_rows, b = loadest.n_cols, c = loadfix.n_cols;
  double epsilon = 0, theta = 0, omega = 0, tau = 0, lambdaAB = 0, abbar = 0, phi12val = 0;
  arma::colvec XYrow (n);
  arma::rowvec Eestrow (n);
  arma::mat secmomest = Eest*Eest.t() + n*condvarest, secmomcross = Efix*Eest.t() + n*condvarcross.t();
  for(int i = 0; i < a; ++i){
    XYrow = XY.col(i);
    phi12val = Phi12(i,i);
    for(int j = 0; j < b; ++j){
        Eestrow = Eest.row(j);
        epsilon = secmomest(j,j);
        theta = arma::dot(Eestrow,XYrow);
        omega = arma::dot(subcolvec(secmomest,j),subrowvec(loadest, i, j));
        for(int k = 0; k < c; k++){
          tau += secmomcross(k,j)*loadfix(i,k);
        }
        abbar = (theta-omega-tau)/(epsilon);
        if(weights(i,j) == 0){
          lambdaAB = (tuningp*phi12val)/(epsilon*(1/pow(1,-5)));
        } 
        else {
          lambdaAB = tuningp*phi12val/(epsilon*myabs(weights(i,j)));
        }
        loadest(i,j) = lasso(abbar, lambdaAB);
    }
  }
  return loadest;
}
