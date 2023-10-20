#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>
#include <string>
#include <list>

#include <RcppDist.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends( RcppDist , RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export]]


NumericMatrix procrustes_cpp(NumericMatrix X, NumericMatrix Y){
  
  
  arma::mat X_aux = as<arma::mat>(X);
  arma::mat Y_aux = as<arma::mat>(Y);
  arma::mat XY = X_aux.t()*Y_aux;
  
  arma::mat U;
  arma::vec s;
  arma::mat V;
  
  arma::svd(U,s,V,XY);
  
  arma::mat A = V*U.t();
  arma::mat result = Y_aux*A;
  
  
  return wrap(result);
}


// [[Rcpp::export]]
NumericVector cseq(double first, double last, double by){
  
  int n = (last - first)/by + 1;
  NumericVector result(n);
  result[0] = first;
  
  
  for(int i = 1; i< n ; i++) {
    
    result[i] = result[i-1] + by;
    
  }
  
  return result;
}



// [[Rcpp::export]]

NumericMatrix procrustes_preprocessing(NumericMatrix ref_position_c, NumericVector ref_position_c_mean , NumericVector result_x ,  NumericVector result_y ,   int N, int iterations ){
  
  NumericMatrix RESULT(N*iterations, 2);
  NumericVector v_x(N*iterations);
  NumericVector v_y(N*iterations);
  int idx = 0;
  
  for(int ite = 0; ite < iterations; ite++) {
    
    
    // double prct = (ite+1)/iterations;
    // printProgress(prct);
    
    NumericVector index =  cseq(0 + idx, (N-1) +idx ,  1);
    
    NumericVector position_x = result_x[index];
    NumericVector position_y = result_y[index];
    
    position_x =  position_x - mean(position_x);
    position_y =  position_y - mean(position_y);
    
    NumericMatrix position(N,2) ;
    
    position(_,0) = position_x;
    position(_,1) = position_y;
    
    position = procrustes_cpp(ref_position_c , position);
    
    position_x =  position(_,0);
    position_y = position(_,1);
    
    position_x = position_x + ref_position_c_mean[0];
    position_y = position_y + ref_position_c_mean[1];
    
    v_x[index] = position_x;
    v_y[index] = position_y;
    
    idx = idx + N;
    Rcout << "The iteration is : " << ite << "\n";
    
  }
  
  
  RESULT(_,0) = v_x;
  RESULT(_,1) = v_y;
  
  return RESULT;
  
}