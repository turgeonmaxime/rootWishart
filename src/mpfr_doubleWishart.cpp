#include <Rcpp.h>
#include <Eigen/Dense>
#include <Eigen/LU>
#include "utility.h"

// [[Rcpp::depends(BH)]]

#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

using namespace boost::multiprecision;
using namespace Rcpp;
using namespace Eigen;
using boost::math::constants::pi;

typedef number<cpp_dec_float<100>> mp_float;

// Declare matrix and vector types with multi-precision scalar type
typedef Eigen::Matrix<mp_float,Dynamic,Dynamic>  MatrixXmp;
typedef Eigen::Matrix<mp_float,Dynamic,1>        VectorXmp;


// [[Rcpp::export]]
double doubleWishart_mp(double x, int s, double m, double n) {
    mp_float xx(x);
    // Initialize vector
    VectorXmp b(s);
    int d = s + (s % 2);

    // Initialize matrix with zeroes on diagonal
    MatrixXmp A(d, d);
    for (int i = 0; i < d; i++) {
        A(i, i) = 0;
    }

    if (s != d) {
        // Fill in extra column
        for (int i = 0; i < s; i++) {
            A(i, s) = boost::math::beta(m + i + 1, n + 1, xx);
            A(s, i) = - A(i, s);
        }
    }


    if (s != 1) {
        for (int i = 0; i < s; i++) {
            b(i) = 0.5 * pow(boost::math::beta(m + i + 1, n + 1, xx), 2);

            for (int j = i; j < (s-1); j++) {
                b(j+1) = ((m+j+1)*b(j) - boost::math::beta(2*m+i+j+2, 2*n+2, xx))/(m+j+n+2);
                A(i,j+1) = boost::math::beta(m+i+1, n+1, xx) * boost::math::beta(m+j+2, n+1, xx) -
                    2*b(j+1);
                A(j+1, i) = - A(i, j+1);
            }
        }
    }


    // Compute scaling constant
    mp_float C1 = 0;
    for (int i = 1; i <= s; i++) {
        C1 += boost::math::lgamma(0.5*(i+2*m+2*n+s+2)) - boost::math::lgamma(0.5*i) -
            boost::math::lgamma(0.5*(i+2*m+1)) - boost::math::lgamma(0.5*(i+2*n+1));
    }
    mp_float C = pow(pi<mp_float>(), 0.5 * s) * exp(C1);

    mp_float det;
    if(d > 4) {
        det = A.fullPivLu().determinant();
    } else {
        det = A.determinant();
    }

    mp_float result = C * sqrt(det);
    return result.convert_to<double>();

}

