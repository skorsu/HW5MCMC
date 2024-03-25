#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

// Update Algorithms: ----------------------------------------------------------
// [[Rcpp::export]]
int update_theta(int theta_old, double lambda1, double lambda2, arma::vec dat){
  
  // Propose a new theta from prior
  Rcpp::IntegerVector thetapvec = Rcpp::sample(111, 1, false);
  int theta_propose = thetapvec[0];
  
  // Calculate the acceptance probability in the log-scale
  double lprob = 0.0;
  lprob += (arma::accu(dat.head_rows(theta_propose)) * std::log(lambda1));
  lprob += (arma::accu(dat.head_rows(112 - theta_propose)) * std::log(lambda2));
  lprob -= (-theta_propose * (lambda1 - lambda2));
  lprob -= (arma::accu(dat.head_rows(theta_old)) * std::log(lambda1));
  lprob -= (arma::accu(dat.head_rows(112 - theta_old)) * std::log(lambda2));
  lprob += (-theta_old * (lambda1 - lambda2));
  
  // Decide
  int theta_new = theta_old;
  double logU = std::log(R::runif(0.0, 1.0));
  if(logU < lprob){
    theta_new = theta_propose;
  }
  
  return theta_new;
  
}

// [[Rcpp::export]]
double update_l1hn(double l1_old, int theta, double s2_1, arma::vec dat){
  
  // Propose a new theta from proposal distribution
  double l1_propose = R::rgamma(l1_old, 1);
  
  // Calculate the acceptance probability in the log-scale
  double lprob = 0.0;
  lprob += (arma::accu(dat.head_rows(theta)) * std::log(l1_propose));
  lprob -= (l1_propose * theta);
  lprob -= (std::pow(l1_propose, 2.0)/(2 * s2_1));
  lprob += R::dgamma(l1_old, l1_propose, 1, 1);
  lprob -= (arma::accu(dat.head_rows(theta)) * std::log(l1_old));
  lprob += (l1_old * theta);
  lprob += (std::pow(l1_old, 2.0)/(2 * s2_1));
  lprob -= R::dgamma(l1_propose, l1_old, 1, 1);
  
  // Decide
  double l1_new = l1_old;
  double logU = std::log(R::runif(0.0, 1.0));
  if(logU < lprob){
    l1_new = l1_propose;
  }
  
  return l1_new;
  
}

// [[Rcpp::export]]
double update_l2hn(double l2_old, int theta, double s2_2, arma::vec dat){
  
  // Propose a new theta from proposal distribution
  double l2_propose = R::rgamma(l2_old, 1);
  
  // Calculate the acceptance probability in the log-scale
  double lprob = 0.0;
  lprob += (arma::accu(dat.head_rows(112 - theta)) * std::log(l2_propose));
  lprob -= (l2_propose * (112 - theta));
  lprob -= (std::pow(l2_propose, 2.0)/(2 * s2_2));
  lprob += R::dgamma(l2_old, l2_propose, 1, 1);
  lprob -= (arma::accu(dat.head_rows(112 - theta)) * std::log(l2_old));
  lprob += (l2_old * (112 - theta));
  lprob += (std::pow(l2_old, 2.0)/(2 * s2_2));
  lprob -= R::dgamma(l2_propose, l2_old, 1, 1);
  
  // Decide
  double l2_new = l2_old;
  double logU = std::log(R::runif(0.0, 1.0));
  if(logU < lprob){
    l2_new = l2_propose;
  }
  
  return l2_new;
  
}

// Gibbs Sampler: --------------------------------------------------------------
// [[Rcpp::export]]
arma::mat gibbsGamma(int iter, arma::vec dat){
  
  arma::mat result(iter, 4, arma::fill::zeros);
  
  // Initialize the parameters from prior
  double ap = R::rgamma(10, 1/10);
  double l1 = R::rgamma(3, 1/ap);
  double l2 = R::rgamma(3, 1/ap);
  Rcpp::IntegerVector thetavec = Rcpp::sample(111, 1, false);
  int tt = thetavec[0];
  
  int tt_new = 0;
  double l1_new = 0.0;
  double l2_new = 0.0;
  double ap_new = 0.0;
  
  for(int i = 0; i < iter; ++i){
    tt_new = update_theta(tt, l1, l2, dat);
    l1_new = R::rgamma(3 + arma::accu(dat.head_rows(tt_new)), 1/(ap + tt_new));
    l2_new = R::rgamma(3 + arma::accu(dat.head_rows(112 - tt_new)), 1/(ap + 112 - tt_new));
    ap_new = R::rgamma(19, 1/(10 + l1_new + l2_new));
    
    result.row(i).col(0).fill(tt_new);
    result.row(i).col(1).fill(l1_new);
    result.row(i).col(2).fill(l2_new);
    result.row(i).col(3).fill(ap_new);
    
    tt = tt_new;
    l1 = l1_new;
    l2 = l2_new;
    ap = ap_new;
    
  }
  
  return result;
  
}

// [[Rcpp::export]]
arma::mat gibbsHalfN(int iter, double s2_1, double s2_2, arma::vec dat){
  
  arma::mat result(iter, 3, arma::fill::zeros);
  
  // Initialize the parameters from prior
  double l1 = R::rgamma(3, 1);
  double l2 = R::rgamma(3, 1);
  Rcpp::IntegerVector thetavec = Rcpp::sample(111, 1, false);
  int tt = thetavec[0];
  
  int tt_new = 0;
  double l1_new = 0.0;
  double l2_new = 0.0;
  
  for(int i = 0; i < iter; ++i){
    tt_new = update_theta(tt, l1, l2, dat);
    l1_new = update_l1hn(l1, tt_new, s2_1, dat);
    l2_new = update_l2hn(l2, tt_new, s2_2, dat);
    
    result.row(i).col(0).fill(tt_new);
    result.row(i).col(1).fill(l1_new);
    result.row(i).col(2).fill(l2_new);
    
    tt = tt_new;
    l1 = l1_new;
    l2 = l2_new;
    
  }
  
  return result;
  
}

// -----------------------------------------------------------------------------




