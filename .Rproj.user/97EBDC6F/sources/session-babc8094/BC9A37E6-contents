#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

// Target Distribution ---------------------------------------------------------
// [[Rcpp::export]]
double tg_dist_ind(double x, double a, double b){
  return std::exp(-std::pow(std::abs(x), a)/b);
}

// [[Rcpp::export]]
arma::vec tg_dist(arma::vec x, double a, double b){
  return arma::exp(-arma::pow(arma::abs(x), a)/b);
}

// Important Sampling ----------------------------------------------------------
// [[Rcpp::export]]
arma::mat isRcpp(unsigned int n, double a, double b, double df){
  
  arma::mat result(n, 5, arma::fill::zeros);
  Rcpp::NumericVector xRcpp = Rcpp::rt(n, df); // Envelope distribution: t with df = 3
  Rcpp::NumericVector fxRcpp = Rcpp::dt(xRcpp, df);
  
  arma::vec x = Rcpp::as<arma::vec>(Rcpp::wrap(xRcpp));
  arma::vec fx = Rcpp::as<arma::vec>(Rcpp::wrap(fxRcpp));
  
  result.col(0) = x;
  result.col(1) = fx;
  result.col(2) = tg_dist(x, a, b);
  result.col(3) = tg_dist(x, a, b)/fx;
  
  double sumW = arma::accu(tg_dist(x, a, b)/fx);
  result.col(4) = result.col(3)/sumW;
  
  return result;
  
}

// Rejection Sampling ----------------------------------------------------------
// [[Rcpp::export]]
arma::vec rsRcpp(unsigned int n, double alp, double a, double b, double df){
  
  arma::vec result(n, arma::fill::zeros);
  
  unsigned int counter = 0;
  double y = 0;
  double u = 0;
  double e = 0;
  
  while(counter < n){
    y = R::rt(df);
    u = R::runif(0.0, 1.0);
    e = R::dt(y, df, 0)/alp;
    if(u <= tg_dist_ind(y, a, b)/e){
      result.row(counter).fill(y);
      counter += 1;
    }
  }
  
  return result;
  
}

// -----------------------------------------------------------------------------




