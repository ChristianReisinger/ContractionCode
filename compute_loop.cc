#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <getopt.h>

#include "fields.hh"
#include "geometry.hh"
#include "linear_algebra.hh"
#include "compute_loop.hh"

using namespace std;

// Computes <chi^dagger, gamma_5 gamma phi>

void compute_loop(complex * CF, int per[], int sig[],
		  int per_5[], int sig_5[],
		  double chi[], double phi[],
		  const unsigned long int T, const unsigned long int L) {

  complex tmp;
  double psi0[24];
  int perm[24], sign[24];
  gamma_eq_gamma_ti_gamma(perm, sign, per, per_5, sig, sig_5);

  for(unsigned long int t = 0; t < T; t++) {
    tmp.re = 0.;
    tmp.im = 0.;
    for(unsigned long int i = 0; i <L*L*L; i++) {
      int v = gsi( t*L*L*L + i );
      fv_eq_gamma_ti_fv(psi0, &phi[ v ], perm, sign);
      co_eq_fv_dag_ti_fv(&tmp, &chi[ v ], psi0);
      CF[t].re += tmp.re;
      CF[t].im += tmp.im;
    }
  }
  return;
}
