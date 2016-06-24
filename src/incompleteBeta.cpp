// Compute the distribution of the largest root in a double Wishart setting
// This is based on Chiani 2016
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector incompleteBeta(NumericVector x, double alpha, double beta) {
    int n = x.size();
  NumericVector res(n);

  for (int i =0; i < n; i++) {
      res[i] = R::pbeta(x[i], alpha, beta, 1, 0) * exp(R::lgammafn(alpha) + R::lgammafn(beta) - R::lgammafn(alpha + beta));
  }
  return res;
}
