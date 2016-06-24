// Compute the distribution of the largest root in a double Wishart setting
// This is based on Chiani 2016

#include <Rmath.h>

double incompleteBeta(double x, double alpha, double beta) {
  double res;
  res = pbeta(x, alpha, beta, 1, 0) * exp(lgammafn(alpha) + lgammafn(beta) - lgammafn(alpha + beta));
  return res;
}
