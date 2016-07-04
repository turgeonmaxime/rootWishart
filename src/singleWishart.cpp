#include <Rcpp.h>
#include <Eigen/Dense>
#include <Eigen/LU>
#include "utility.h"

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
double singleWishart_C(double x, int n_min, int n_max) {
    double b;
    int n_mat = n_min + (n_min % 2);
    double alpha = 0.5*(n_max - n_min - 1);

    // Initialize matrix with zeroes on diagonal
    MatrixXd A(n_mat, n_mat);
    for (int i = 0; i < n_mat; i++) {
        A(i, i) = 0;
    }

    // Pre-compute some values
    VectorXd p(n_min);
    VectorXd g(n_min);
    VectorXd q(2*n_min - 2);

    if (n_min != n_mat) {
        // Fill in extra column
        VectorXd alpha_vec(n_min);
        double alpha_last = 0.5*(n_max + n_min + 1);
        for (int i = 0; i < n_min; i++) {
            alpha_vec(i) = alpha + i + 1;
            A(i, n_min) = pow(2, -alpha_last)*R::pgamma(0.5*x, alpha_vec(i), 1.0, 1, 0)/R::gammafn(alpha_last);
            A(n_min, i) = - A(i, n_min);
        }
    }

    for (int i = 0; i < n_min; i++) {
        p(i) = R::pgamma(0.5*x, alpha + i + 1, 1.0, 1, 0);
        g(i) = R::gammafn(alpha + i + 1);
    }
    for (int i = 0; i < (2*n_min - 2); i++) {
        q(i) = pow(0.5, 2*alpha + i + 2) * R::gammafn(2*alpha+i+2)*R::pgamma(x, 2*alpha + i + 2, 1.0, 1, 0);
    }

    for (int i = 0; i < n_min; i++) {
        b = 0.5*p(i)*p(i);
        for(int j = i; j < (n_min - 1); j++) {
            b -= q(i+j)/(g(i)*g(j+1));
            A(i, j+1) = p(i)*p(j+1) - 2*b;
            A(j+1, i) = -A(i, j+1);
        }
    }

    // Compute constant
    double K1 = pow(M_PI, 0.5*n_min*n_min);
    K1 /= pow(2, 0.5*n_min*n_max)*mgamma_C(0.5*n_max, n_min, FALSE)*mgamma_C(0.5*n_min, n_min, FALSE);
    double K2 = pow(2, alpha*n_mat+0.5*n_mat*(n_mat+1));
    for (int k = 0; k < n_mat; k++) {
        K2 *= R::gammafn(alpha + k + 1);
    }

    // Compute Pfaffian
    double det;
    if(n_mat > 4) {
        det = A.fullPivLu().determinant();
    } else {
        det = A.determinant();
    }

    double result = K1 * K2 * sqrt(det);
    return result;
}
