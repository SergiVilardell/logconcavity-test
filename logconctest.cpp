#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int logconctest_cpp(NumericVector x, NumericVector x_g) {
  
  int n = x.size() -1;
  NumericVector out(n);
  int i;  

  for(i = 0; i < n-1; i++){
    
   if( (x_g[i+1]-x_g[i] == 0) || (x_g[i+2]-x_g[i+1] == 0)){continue;}
   
   out[i] = (log(x[i+1])-log(x[i]))/(x_g[i+1]-x_g[i]) ;
   out[i+1] = (log(x[i+2])-log(x[i+1]))/(x_g[i+2]-x_g[i+1]);
   if((out[i+1] - out[i])>0){return 0;}
  }
  return 1;
}