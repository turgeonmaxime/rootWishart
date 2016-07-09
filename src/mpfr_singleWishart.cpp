#include <Rcpp.h>
#include <Eigen/Dense>
#include <Eigen/LU>
#include "mpfr_pgamma.h"
#include "MPRealSupport.h"

using namespace mpfr;
using namespace Eigen;

// [[Rcpp::export]]
double singleWishart_mpfr(double xx, int n_min, int n_max) {
    // set precision to 256 bits
    mpreal::set_default_prec(256);
    mpreal x = xx;

    const mpreal pi = mpfr::const_pi();
    mpreal b;
    int n_mat = n_min + (n_min % 2);
    mpreal alpha = 0.5*(n_max - n_min - 1);

    // Initialize matrix with zeroes on diagonal
    Matrix<mpreal,Dynamic,Dynamic> A(n_mat, n_mat);
    for (int i = 0; i < n_mat; i++) {
        A(i, i) = 0;
    }

    // Pre-compute some values
    Matrix<mpreal,Dynamic,1> p(n_min);
    Matrix<mpreal,Dynamic,1> g(n_min);
    Matrix<mpreal,Dynamic,1> q(2*n_min - 2);

    if (x < 1) {
        for (int i = 0; i < n_min; i++) {
            p(i) = pgamma_smallx(0.5*x, alpha + i + 1, TRUE, FALSE);
            g(i) = gamma(alpha + i + 1);
            q(i) = pow(0.5, 2*alpha + i + 2) * gamma(2*alpha+i+2)*pgamma_smallx(x, 2*alpha + i + 2, TRUE, FALSE);
        }
        for (int i = n_min; i < (2*n_min - 2); i++) {
            q(i) = pow(0.5, 2*alpha + i + 2) * gamma(2*alpha+i+2)*pgamma_smallx(x, 2*alpha + i + 2, TRUE, FALSE);
        }
    } else {
        create_mpfr_p(p, x, n_max, n_min);
        create_mpfr_q(q, x, n_max, n_min);
        for (int i = 0; i < n_min; i++) {
            g(i) = gamma(alpha + i + 1);
        }
    }

    // std::cout << p << std::endl;
    // std::cout << g << std::endl;
    // std::cout << q << std::endl;

    if (n_min != n_mat) {
        // Fill in extra column
        Matrix<mpreal,Dynamic,1> alpha_vec(n_min);
        mpreal alpha_last = 0.5*(n_max + n_min + 1);
        for (int i = 0; i < n_min; i++) {
            alpha_vec(i) = alpha + i + 1;
            A(i, n_min) = pow(2, -alpha_last)*p(i)/gamma(alpha_last);
            A(n_min, i) = - A(i, n_min);
        }
    }

    // std::cout << A << std::endl;

    for (int i = 0; i < n_min; i++) {
        b = 0.5*p(i)*p(i);
        for(int j = i; j < (n_min - 1); j++) {
            b -= q(i+j)/(g(i)*g(j+1));
            A(i, j+1) = p(i)*p(j+1) - 2*b;
            A(j+1, i) = -A(i, j+1);
        }
    }

    // std::cout << A << std::endl;

    // Compute constant
    mpreal K1 = pow(pi, 0.5*n_min*n_min);
    K1 /= pow(2, 0.5*n_min*n_max) * mgamma_mpfr(mpreal(0.5*n_max), n_min) * mgamma_mpfr(mpreal(0.5*n_min), n_min);
    mpreal K2 = pow(2, alpha*n_mat+0.5*n_mat*(n_mat+1));
    for (int k = 0; k < n_mat; k++) {
        K2 *= gamma(alpha + k + 1);
    }

    // Compute Pfaffian
    mpreal det;
    if(n_mat > 4) {
        det = A.fullPivLu().determinant();
    } else {
        det = A.determinant();
    }

    mpreal result = K1 * K2 * sqrt(det);
    return result.toDouble();
}
