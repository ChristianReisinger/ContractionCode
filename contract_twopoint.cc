#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdlib.h>
// #include <getopt.h>

#include "fields.hh"
#include "geometry.hh"
#include "linear_algebra.hh"
#include "contract_twopoint.hh"

using namespace std;

// Performs contraction of 
// <gamma_nu gamma_5 chi^dagger gamma_5 gamma_mu phi>
// for all t
//
// this routine is for gamma_5 diagonal only!
void contract_twopoint(complex * CF, const int gamma_source, const int gamma_sink,
		       double * chi[], double * phi[],
		       const unsigned int T, const unsigned int L, 
		       const unsigned int ts) {
  complex tmp;
  double psi0[24], psi1[24];
  int tt=0;
  int isimag = (gamma_permutation[gamma_source][0])%2;
  int p[4];

  for(int mu = 0; mu < 4; mu++) {
    p[mu] = gamma_permutation[gamma_source][6*mu]/6;
  }
  for(unsigned int t = 0; t < T; t++) {
    tmp.re = 0.;
    tmp.im = 0.;
    
    for(unsigned int i = 0; i < L*L*L; i++) {
      unsigned long int v = gsi( get_index(ts+t, 0, 0, 0, T, L) + i );
      for(unsigned int mu = 0; mu < 4; mu++) {
	// here comes the gamma_sink and gamma_5 at sink
	fv_eq_gamma_ti_fv(psi0, gamma_sink, &phi[ mu ][ v ]);
	fv_eq_gamma_ti_fv(psi1, 5, psi0);
	
	// spinor index p[mu] takes care of gamma_source at source
	//co_eq_fv_dag_ti_fv(&tmp, psi1, &chi[ p ][ v ]);
	co_eq_fv_dag_ti_fv(&tmp, &chi[ p[mu] ][ v ], psi1);
	
	// sign at source (also gamma_5) and possible factor of i
	if( !isimag ) {
	  CF[tt].re += gamma_sign[5][6*mu]*gamma_sign[gamma_source][6*mu]*tmp.re;
	  CF[tt].im += gamma_sign[5][6*mu]*gamma_sign[gamma_source][6*mu]*tmp.im;
	}
	else {
	  CF[tt].re += gamma_sign[5][6*mu]*gamma_sign[gamma_source][6*mu]*tmp.im;
	  CF[tt].im += gamma_sign[5][6*mu]*gamma_sign[gamma_source][6*mu]*tmp.re;
	}
      }
    }
    tt++;
  }
  return;
}

void contract_twopoint(complex * CF, const int gamma_source, const int gamma_sink,
		       double * chi[], double * phi[],
		       const unsigned int T, const unsigned int L, 
		       const unsigned int ts, 
		       int * const perm[], int * const sign[]) {
  complex tmp;
  double psi0[24], psi1[24];
  int tt=0;
  int isimag = (perm[gamma_source][0])%2;
  int p[4];

  for(unsigned int mu = 0; mu < 4; mu++) {
    p[mu] = perm[gamma_source][6*mu]/6;
  }
  for(unsigned int t = 0; t < T; t++) {
    tmp.re = 0.;
    tmp.im = 0.;
    
    for(unsigned int i = 0; i < L*L*L; i++) {
      unsigned long int v = gsi( get_index(ts+t, 0, 0, 0, T, L) + i );
      for(unsigned int mu = 0; mu < 4; mu++) {
	
	// here comes the gamma_sink and gamma_5 at sink
	fv_eq_gamma_ti_fv(psi0, &phi[ mu ][ v ], perm[gamma_sink], sign[gamma_sink]);
	fv_eq_gamma_ti_fv(psi1, psi0, perm[5], sign[5]);
	
	// spinor index p[mu] takes care of gamma_source at source
	//co_eq_fv_dag_ti_fv(&tmp, psi1, &chi[ p ][ v ]);
	co_eq_fv_dag_ti_fv(&tmp, &chi[ p[mu] ][ v ], psi1);
	
	// sign at source (also gamma_5) and possible factor of i
	if( !isimag ) {
	  CF[tt].re += sign[5][6*mu]*sign[gamma_source][6*mu]*tmp.re;
	  CF[tt].im += sign[5][6*mu]*sign[gamma_source][6*mu]*tmp.im;
	}
	else {
	  CF[tt].re += sign[5][6*mu]*sign[gamma_source][6*mu]*tmp.im;
	  CF[tt].im += sign[5][6*mu]*sign[gamma_source][6*mu]*tmp.re;
	}
      }
    }
    tt++;
  }
  return;
}


// next version is for general gamma_5
// should compute every gamma matrix combination at source and sink 
// using gamma_5 hermiticity
// this version can also cope with 12 components instead of only four

void contract_twopoint(complex * CF, int const per_source[], int const sig_source[],
		       int per_sink[], int sig_sink[],
		       int per_5[], int sig_5[],
		       double * chi[], double * phi[],
		       const unsigned int T, const unsigned int L, 
		       const unsigned int ts, const unsigned int n_c) {
  complex tmp;
  double psi0[24];
  unsigned int tt=0;
  int isimag = (per_source[0])%2;
  int p[4];
  int psource[24], psink[24], ssource[24], ssink[24];

  // compute gamma_source gamma_5
  gamma_eq_gamma_ti_gamma(psource, ssource, per_source, per_5, sig_source, sig_5);
  // and gamma_5 gamma_sink
  gamma_eq_gamma_ti_gamma(psink, ssink, per_5, per_sink, sig_5, sig_sink);

  isimag = psource[0]%2;
  for(unsigned int mu = 0; mu < 4; mu++) {
    p[mu] = psource[6*mu]/6;
  }

  if(n_c == 1) {
    for(unsigned int t = 0; t < T; t++) {
      tmp.re = 0.;
      tmp.im = 0.;
      
      for(unsigned int i = 0; i < L*L*L; i++) {
	unsigned long int v = gsi( get_index(ts+t, 0, 0, 0, T, L) + i );
	for(unsigned int mu = 0; mu < 4; mu++) {
	  
	  // here comes the gamma_sink * gamma_5 at sink
	  fv_eq_gamma_ti_fv(psi0, &phi[ mu ][ v ], psink, ssink);
	  
	  // spinor index p[mu] takes care of gamma_source at source
	  co_eq_fv_dag_ti_fv(&tmp, &chi[ p[mu] ][ v ], psi0);
	  
	  // sign at source and possible factor of i
	  if( !isimag ) {
	    CF[tt].re += ssource[6*mu]*tmp.re;
	    CF[tt].im += ssource[6*mu]*tmp.im;
	  }
	  else {
	    CF[tt].re += ssource[6*mu]*tmp.im;
	    CF[tt].im += ssource[6*mu]*tmp.re;
	  }
	}
      }
      tt++;
    }
  }
  else {
    for(unsigned int t = 0; t < T; t++) {
      tmp.re = 0.;
      tmp.im = 0.;
      
      for(unsigned int i = 0; i < L*L*L; i++) {
	unsigned long int v = gsi( get_index(ts+t, 0, 0, 0, T, L) + i );
	// spin indices
	for(unsigned int mu = 0; mu < 4; mu++) {
	  // colour indices
	  for(unsigned int c = 0; c < n_c; c++) {
	    
	    // here comes the gamma_sink * gamma_5 at sink
	    fv_eq_gamma_ti_fv(psi0, &phi[ mu*n_c + c ][ v ], psink, ssink);
	    
	    // spinor index p[mu] takes care of gamma_source at source
	    co_eq_fv_dag_ti_fv(&tmp, &chi[ p[mu]*n_c + c ][ v ], psi0);

	    // sign at source and possible factor of i
	    if( !isimag ) {
	      CF[tt].re += ssource[6*mu]*tmp.re;
	      CF[tt].im += ssource[6*mu]*tmp.im;
	    }
	    else {
	      CF[tt].re += ssource[6*mu]*tmp.im;
	      CF[tt].im += ssource[6*mu]*tmp.re;
	    }
	  }
	}
      }
      tt++;
    }
  }
  return;
}
