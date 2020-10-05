#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::mat mmult(arma::mat A, arma::mat B) {
  arma::mat res = A * B;
  return res;
}

// [[Rcpp::export]]
arma::mat minv(arma::mat M) {
  arma::mat res = inv(M);
  return res;
}
