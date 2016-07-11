#include <Rcpp.h>
#include <Eigen/Dense>
#include <Eigen/LU>

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
// typedef Eigen::Matrix<mp_float,Dynamic,Dynamic>  MatrixXmp;
// typedef Eigen::Matrix<mp_float,Dynamic,1>        VectorXmp;


template <class T>
T doubleWishart(T xx, int s, T mm, T nn) {
    // Initialize vector
    Eigen::Matrix<T,Dynamic,1> b(s);
    int d = s + (s % 2);

    // Initialize matrix with zeroes on diagonal
    Eigen::Matrix<T,Dynamic,Dynamic> A(d, d);
    for (int i = 0; i < d; i++) {
        A(i, i) = 0;
    }

    if (s != d) {
        // Fill in extra column
        for (int i = 0; i < s; i++) {
            A(i, s) = boost::math::beta(mm + i + 1, nn + 1, xx);
            A(s, i) = - A(i, s);
        }
    }


    if (s != 1) {
        for (int i = 0; i < s; i++) {
            b(i) = 0.5 * pow(boost::math::beta(mm + i + 1, nn + 1, xx), 2);

            for (int j = i; j < (s-1); j++) {
                b(j+1) = ((mm+j+1)*b(j) - boost::math::beta(2*mm+i+j+2, 2*nn+2, xx))/(mm+j+nn+2);
                A(i,j+1) = boost::math::beta(mm+i+1, nn+1, xx) * boost::math::beta(mm+j+2, nn+1, xx) -
                    2*b(j+1);
                A(j+1, i) = - A(i, j+1);
            }
        }
    }


    // Compute scaling constant
    T C1 = 0;
    for (int i = 1; i <= s; i++) {
        C1 += boost::math::lgamma(0.5*(i+2*mm+2*nn+s+2)) - boost::math::lgamma(0.5*i) -
            boost::math::lgamma(0.5*(i+2*mm+1)) - boost::math::lgamma(0.5*(i+2*nn+1));
    }
    T C = pow(pi<T>(), 0.5 * s) * exp(C1);

    T det;
    if(d > 4) {
        det = A.fullPivLu().determinant();
    } else {
        det = A.determinant();
    }

    T result = C * sqrt(det);
    return result;
}

// [[Rcpp::export]]
double doubleWishart_raw(double x, int s, double m, double n, bool mp) {
    double result;
    if (mp) {
        mp_float xx(x);
        mp_float mm(m);
        mp_float nn(n);

        result = doubleWishart(xx, s, mm, nn).convert_to<double>();
    } else {
        result = doubleWishart(x, s, m, n);
    }

    return result;
}

