#include <Eigen/Dense>
#include "mpfr_pgamma.h"
#include "MPRealSupport.h"
#include <Rcpp.h>

using namespace mpfr;
using namespace Eigen;

void create_mpfr_p(Matrix<mpreal,Dynamic,1>& p, mpreal x, int n_max, int n_min) {
    mpreal::set_default_prec(256);
    const mpreal pi = mpfr::const_pi();

    int n_diff = n_max - n_min;

    if(n_diff % 2) { /* ODD */
        // std::cout << "ODD" << std::endl;
        int alph = 0.5*(n_diff - 1);
        if (n_diff == 1) {
            // std::cout << "DIFF1 " << alph<< std::endl;
            p(0) = 1 - exp(-0.5*x);
        } else {
            // std::cout << "DIFFNOT1" << std::endl;
            mpreal sum = "0.0";
            for (int i = 0; i <= alph; i++) {
                sum += exp(-0.5*x)*pow(0.5*x, i)/mpfr::gamma(i+1);
            }
            p(0) = 1 - sum;
        }

        for (int i = 1; i < n_min; i++) {
            // std::cout << exp(-0.5*x) << " " << pow(0.5*x, alph + i) << " " << gamma(alph + i) << std::endl;
            p(i) = p(i-1) - exp(-0.5*x)*pow(0.5*x, alph + i)/mpfr::gamma(alph + i + 1);
        }
    } else { /* EVEN */
        // std::cout << "EVEN" << std::endl;
        int n = floor(0.5*(n_diff - 1));
        if (n_max == n_min) {
            p(0) = sqrt(pi) * erf(sqrt(0.5*x))/mpfr::gamma(0.5);
        } else {
            mpreal sum = "0.0";
            for (int i = 0; i <= n; i++) {
                sum += exp(-0.5*x)*pow(0.5*x, i+0.5)/((i + 0.5)*mpfr::gamma(i+0.5));
            }
            p(0) = erf(sqrt(0.5*x)) - sum;
        }

        for (int i = 1; i < n_min; i++) {
            p(i) = p(i-1) - exp(-0.5*x)*pow(0.5*x, n + i + 0.5)/((n + i + 0.5)*mpfr::gamma(n + i + 0.5));
        }
    }

}

void create_mpfr_q(Matrix<mpreal,Dynamic,1>& q, mpreal x, int n_max, int n_min) {
    int n_diff = n_max - n_min;

    if (n_max == n_min) {
        q(0) = 1 - exp(-x);
    } else {
        mpreal sum = "0.0";
        for (int i = 0; i <= n_diff; i++) {
            sum += exp(-x)*pow(x, i)/mpfr::gamma(i+1);
        }
        q(0) = 1 - sum;
    }

    for (int i = 1; i < (2*n_min - 2); i++) {
        q(i) = q(i-1) - exp(-x)*pow(x, n_diff + i)/mpfr::gamma(n_diff + i + 1);
        q(i-1) *= pow(0.5, n_diff + i) * mpfr::gamma(n_diff + i); // THINK ABOUT REPLACING pow() BY BIT SHIFTING
    }
    q(2*n_min - 3) *= pow(0.5, n_max+n_min-2) * mpfr::gamma(n_max+n_min-2); // THINK ABOUT REPLACING pow() BY BIT SHIFTING
}

// Abramowitz and Stegun 6.5.29 [right]
mpreal pgamma_smallx (mpreal x, mpreal alph, bool lower_tail, bool log_p) {
    // set precision to 256 bits
    mpreal::set_default_prec(256);
    const mpreal EPSILON  = std::numeric_limits<mpreal>::epsilon();

    double n = 0;
    mpreal sum = 0, c = alph, term;

    /* Relative to 6.5.29 all terms have been multiplied by alph
    and the first, thus being 1, is omitted. */

    do {
    n++;
    c *= -x / n;
    term = c / (alph + n);
    sum += term;
    } while (fabs (term) > EPSILON * fabs (sum));

    if (lower_tail) {
    mpreal f1 = log_p ?  log1p(sum) : 1 + sum;
    mpreal f2;
    if (alph > 1) {
        f2 = mpfr_dpois (alph, x, log_p);
        f2 = log_p ? f2 + x : f2 * exp (x);
    } else if (log_p)
        f2 = alph * log (x) - lngamma(alph+1);
    else
        f2 = pow (x, alph) / gamma(alph + 1);
    return log_p ? f1 + f2 : f1 * f2;
    } else {
    mpreal lf2 = alph * log (x) - lngamma(alph + 1);
    if (log_p)
        return log(1 - exp(log1p(sum) + lf2));
        // return R_Log1_Exp (log1p (sum) + lf2);
    else {
        mpreal f1m1 = sum;
        mpreal f2m1 = expm1(lf2);
        return -(f1m1 + f2m1 + f1m1 * f2m1);
    }
    }
}

mpreal mpfr_dpois(mpreal x, mpreal lambda, bool give_log) {
    if(give_log) {
        return(-lambda + x*log(lambda) -lgamma(x+1));
    } else {
        return(exp(-lambda + x*log(lambda) -lgamma(x+1)));
    }
}

mpreal mgamma_mpfr(mpreal x, int m) {
    const mpreal pi = mpfr::const_pi();
    mpreal res;
    res = pow(pi, 0.25*m*(m-1));
    for(int i = 0; i < m; i++) {
        res *= gamma(x - 0.5*i);
    }
    return res;
}
