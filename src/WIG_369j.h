//https://www.gnu.org/software/gsl/doc/html/specfunc.html?highlight=3j#coupling-coefficients
#ifndef _WIGNER_H
#define _WIGNER_H
#include <gsl/gsl_sf_coupling.h>
#include <cmath>

double WIG_3j_1(int j1, int j2, int j3, int m1, int m2, int m3)
double WIG_3j_2(int two_j1, int two_j2, int two_j3, int two_m1, int two_m2,
  int two_m3)

double WIG_cg_1(int j1, int m1, int j1, int m2, int J, int M)
double WIG_cg_2(int two_j1, int two_m1, int two_j1, int two_m2, int two_J,
  int two_M);


#endif