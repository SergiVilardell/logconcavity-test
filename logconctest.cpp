#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int logconctest_cpp(NumericVector x) {
  
  int n = x.size() -1;
  NumericVector out(n);
  int i;  

  for(i = 0; i < n-1; i++){
   out[i] = log(x[i])-log(x[i+1]);
   out[i + 1] = log(x[i+1])-log(x[i+2]);
   if(out[i]-out[i+1]>0){return 0;}
  }
  return 1;
}