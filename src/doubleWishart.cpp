#include <Rcpp.h>
#include <Eigen/Dense>
#include <Eigen/LU>
#include "incompleteBeta.h"

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
double doubleWishart_C(double x, int s, double m, double n) {
    // Initialize vector
    VectorXd b(s);
    int d = s + (s % 2);

    // Initialize matrix with zeroes on diagonal
    MatrixXd A(d, d);
    for (int i = 0; i < d; i++) {
        A(i, i) = 0;
    }


    if (s != d) {
        // Fill in extra column
        for (int i = 0; i < s; i++) {
            A(i, s) = incompleteBeta_C(x, m + i + 1, n + 1);
            A(s, i) = - A(i, s);
        }
    }


    if (s != 1) {
        for (int i = 0; i < s; i++) {
            b(i) = 0.5 * pow(incompleteBeta_C(x, m + i + 1, n + 1), 2);

            for (int j = i; j < (s-1); j++) {
                b(j+1) = ((m+j+1)*b(j) - incompleteBeta_C(x, 2*m+i+j+2, 2*n+2))/(m+j+n+2);
                A(i,j+1) = incompleteBeta_C(x, m+i+1, n+1) * incompleteBeta_C(x, m+j+2, n+1) -
                    2*b(j+1);
                A(j+1, i) = - A(i, j+1);
            }
        }
    }


    // Compute scaling constant
    double C1 = 0;
    for (int i = 1; i <= s; i++) {
        C1 += lgamma(0.5*(i+2*m+2*n+s+2)) - lgamma(0.5*i) -lgamma(0.5*(i+2*m+1)) - lgamma(0.5*(i+2*n+1));
    }
    double C = pow(M_PI, 0.5 * s) * exp(C1);

    double det;
    if(d > 4) {
        det = A.fullPivLu().determinant();
    } else {
        det = A.determinant();
    }

    double result = C * sqrt(det);
    return result;

}

