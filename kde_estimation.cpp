#include <Rcpp.h>

// [[Rcpp::plugins(openmp)]]
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector kde_cpp(NumericVector s, NumericVector x, double h) {
  
  int n = s.size();
  int m = x.size();
  NumericVector out(m);
  double sum;
  int j;
  int i;  
  #pragma omp parallel for private(i, j, sum) num_threads( 3 )
  for(j = 0; j < m; j++){
    sum = 0;
    for(i = 0; i < n; i++){
      sum += exp(-(pow((x[j]-s[i])/h,2))/2)/sqrt(2*M_PI);
    }
    out[j] = sum/(n*h);
  }
  return out;
}
