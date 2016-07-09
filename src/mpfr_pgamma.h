#ifndef MPFRPGAMMA_H
#define MPFRPGAMMA_H
#include "mpreal.h"

using mpfr::mpreal;

void create_mpfr_p(Eigen::Matrix<mpreal, Eigen::Dynamic,1>&, mpreal, int, int);
void create_mpfr_q(Eigen::Matrix<mpreal, Eigen::Dynamic,1>&, mpreal, int, int);

mpreal pgamma_smallx(mpreal, mpreal, bool, bool) ;
mpreal mpfr_dpois(mpreal, mpreal, bool);

mpreal mgamma_mpfr(mpreal, int);

#endif
