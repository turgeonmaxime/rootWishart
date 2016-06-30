// Compute the distribution of the largest root in a double Wishart setting
// This is based on Chiani 2016
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double incompleteBeta_C(double x, double alpha, double beta) {
  //   int n = x.size();
  // NumericVector res(n);
  //
  // for (int i =0; i < n; i++) {
  //     res[i] = R::pbeta(x[i], alpha, beta, 1, 0) * exp(R::lgammafn(alpha) + R::lgammafn(beta) - R::lgammafn(alpha + beta));
  // }

  double res = R::pbeta(x, alpha, beta, 1, 0) * exp(R::lgammafn(alpha) + R::lgammafn(beta) - R::lgammafn(alpha + beta));
  return res;
}

// [[Rcpp::export]]
double mgamma_C(double x, int m, bool logar) {

    double res;

    if (logar) {
        res = 0.25*m*(m-1) * log(M_PI);
        for(int i = 0; i < m; i++) {
            res += R::lgammafn(x - 0.5*i);
        }
    } else {
        res = pow(M_PI, 0.25*m/(m-1));
        for(int i = 0; i < m; i++) {
            res *= R::gammafn(x - 0.5*i);
        }
    }

    return res;
}
