// ********************



// Q_phi.cc

// Author: Marc Wagner
// Date: September 2007
// Carsten Urbach


// ********************



#include "Q_phi.hh"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.hh"
#include "linear_algebra.hh"



// ********************



// Computes xi = Q phi, where Q is the tm Dirac operator.
// Q is written as Q = m_0 + i mu gamma_5 + ...
// not in the hopping parameter representation
//
// if diag = true only the hopping part is computed and
// the result is then in the hopping parameter normalisation
// note that there might be a relative "-" compared to the 
// HMC hopping matrix convention

void Q_phi(double *xi, double *phi, double *gauge_field, int T, int L, double kappa, double mu, bool diag)
{
  int it, ix, iy, iz;
  double SU3_1[18];
  double spinor1[24], spinor2[24];


  double _1_2_kappa = 0.5 / kappa;

  for(it = 0; it < T; it++)
    {
      fprintf(stderr, "void Q_phi(...   -->   it = %2d ...\n", it);

      double angle_B = M_PI / (double)T;

      complex B;
      B.re = cos(angle_B);
      B.im = sin(angle_B);

      for(ix = 0; ix < L; ix++)
	{
	  for(iy = 0; iy < L; iy++)
	    {
	      for(iz = 0; iz < L; iz++)
		{
		  int index_s = gsi(get_index(it, ix, iy, iz, T, L));

		  double *xi_ = xi + index_s;
		  double *phi_;
		  double *U_;


		  fv_eq_zero(xi_);


		  // Negative t-direction.

		  phi_ = phi + gsi(get_index(it-1, ix, iy, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 0, phi_);
		  fv_mi(spinor1);

		  fv_mi_eq_fv(spinor1, phi_);

		  U_ = gauge_field + ggi(get_index(it-1, ix, iy, iz, T, L), 0);

		  cm_eq_cm_ti_co(SU3_1, U_, &B);
		  fv_pl_eq_cm_dag_ti_fv(xi_, SU3_1, spinor1);


		  // Positive t-direction.

		  phi_ = phi + gsi(get_index(it+1, ix, iy, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 0, phi_);

		  fv_mi_eq_fv(spinor1, phi_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 0);

		  cm_eq_cm_ti_co(SU3_1, U_, &B);
		  fv_pl_eq_cm_ti_fv(xi_, SU3_1, spinor1);


		  // Negative x-direction.

		  phi_ = phi + gsi(get_index(it, ix-1, iy, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 1, phi_);
		  fv_mi(spinor1);

		  fv_mi_eq_fv(spinor1, phi_);

		  U_ = gauge_field + ggi(get_index(it, ix-1, iy, iz, T, L), 1);

		  fv_pl_eq_cm_dag_ti_fv(xi_, U_, spinor1);


		  // Positive x-direction.

		  phi_ = phi + gsi(get_index(it, ix+1, iy, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 1, phi_);

		  fv_mi_eq_fv(spinor1, phi_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 1);

		  fv_pl_eq_cm_ti_fv(xi_, U_, spinor1);


		  // Negative y-direction.

		  phi_ = phi + gsi(get_index(it, ix, iy-1, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 2, phi_);
		  fv_mi(spinor1);

		  fv_mi_eq_fv(spinor1, phi_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy-1, iz, T, L), 2);

		  fv_pl_eq_cm_dag_ti_fv(xi_, U_, spinor1);


		  // Positive y-direction.

		  phi_ = phi + gsi(get_index(it, ix, iy+1, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 2, phi_);

		  fv_mi_eq_fv(spinor1, phi_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 2);

		  fv_pl_eq_cm_ti_fv(xi_, U_, spinor1);


		  // Negative z-direction.

		  phi_ = phi + gsi(get_index(it, ix, iy, iz-1, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 3, phi_);
		  fv_mi(spinor1);

		  fv_mi_eq_fv(spinor1, phi_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy, iz-1, T, L), 3);

		  fv_pl_eq_cm_dag_ti_fv(xi_, U_, spinor1);


		  // Positive z-direction.

		  phi_ = phi + gsi(get_index(it, ix, iy, iz+1, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 3, phi_);

		  fv_mi_eq_fv(spinor1, phi_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 3);

		  fv_pl_eq_cm_ti_fv(xi_, U_, spinor1);


		  if(diag) {
		    // Multiplication with 1/2.
		    
		    fv_ti_eq_re(xi_, 0.5);
		    
		    // Diagonal elements.

		    phi_ = phi + index_s;
		    
		    fv_eq_fv_ti_re(spinor1, phi_, _1_2_kappa);
		    fv_pl_eq_fv(xi_, spinor1);
		    
		    fv_eq_gamma_ti_fv(spinor1, 5, phi_);
		    fv_eq_fv_ti_im(spinor2, spinor1, mu);
		    fv_pl_eq_fv(xi_, spinor2);
		  }
		  else {
		    fv_ti_eq_re(xi_, kappa);
		  }
		}
	    }
	}
    }
}



// ********************



void Q_phi_timeslice(double *xi, double *phi, double *gauge_field, int T, int L, double kappa, double mu, int timeslice, int T_B)
{
  int it, ix, iy, iz;
  double SU3_1[18];
  double spinor1[24], spinor2[24];


  double _1_2_kappa = 0.5 / kappa;

  it = timeslice;

  fprintf(stderr, "void Q_phi(...   -->   it = %2d ...\n", it);

  double angle_B = M_PI / (double)T_B;

  complex B;
  B.re = cos(angle_B);
  B.im = sin(angle_B);

  for(ix = 0; ix < L; ix++)
    {
      for(iy = 0; iy < L; iy++)
	{
	  for(iz = 0; iz < L; iz++)
	    {
	      int index_s = gsi(get_index(it, ix, iy, iz, T, L));

	      double *xi_ = xi + gsi(get_index(0, ix, iy, iz, T, L));
	      double *phi_;
	      double *U_;


	      fv_eq_zero(xi_);


	      // Negative t-direction.

	      phi_ = phi + gsi(get_index(it-1, ix, iy, iz, T, L));

	      fv_eq_gamma_ti_fv(spinor1, 0, phi_);
	      fv_mi(spinor1);

	      fv_mi_eq_fv(spinor1, phi_);

	      U_ = gauge_field + ggi(get_index(it-1, ix, iy, iz, T, L), 0);

	      cm_eq_cm_ti_co(SU3_1, U_, &B);
	      fv_pl_eq_cm_dag_ti_fv(xi_, SU3_1, spinor1);


	      // Positive t-direction.

	      phi_ = phi + gsi(get_index(it+1, ix, iy, iz, T, L));

	      fv_eq_gamma_ti_fv(spinor1, 0, phi_);

	      fv_mi_eq_fv(spinor1, phi_);

	      U_ = gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 0);

	      cm_eq_cm_ti_co(SU3_1, U_, &B);
	      fv_pl_eq_cm_ti_fv(xi_, SU3_1, spinor1);


	      // Negative x-direction.

	      phi_ = phi + gsi(get_index(it, ix-1, iy, iz, T, L));

	      fv_eq_gamma_ti_fv(spinor1, 1, phi_);
	      fv_mi(spinor1);

	      fv_mi_eq_fv(spinor1, phi_);

	      U_ = gauge_field + ggi(get_index(it, ix-1, iy, iz, T, L), 1);

	      fv_pl_eq_cm_dag_ti_fv(xi_, U_, spinor1);


	      // Positive x-direction.

	      phi_ = phi + gsi(get_index(it, ix+1, iy, iz, T, L));

	      fv_eq_gamma_ti_fv(spinor1, 1, phi_);

	      fv_mi_eq_fv(spinor1, phi_);

	      U_ = gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 1);

	      fv_pl_eq_cm_ti_fv(xi_, U_, spinor1);


	      // Negative y-direction.

	      phi_ = phi + gsi(get_index(it, ix, iy-1, iz, T, L));

	      fv_eq_gamma_ti_fv(spinor1, 2, phi_);
	      fv_mi(spinor1);

	      fv_mi_eq_fv(spinor1, phi_);

	      U_ = gauge_field + ggi(get_index(it, ix, iy-1, iz, T, L), 2);

	      fv_pl_eq_cm_dag_ti_fv(xi_, U_, spinor1);


	      // Positive y-direction.

	      phi_ = phi + gsi(get_index(it, ix, iy+1, iz, T, L));

	      fv_eq_gamma_ti_fv(spinor1, 2, phi_);

	      fv_mi_eq_fv(spinor1, phi_);

	      U_ = gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 2);

	      fv_pl_eq_cm_ti_fv(xi_, U_, spinor1);


	      // Negative z-direction.

	      phi_ = phi + gsi(get_index(it, ix, iy, iz-1, T, L));

	      fv_eq_gamma_ti_fv(spinor1, 3, phi_);
	      fv_mi(spinor1);

	      fv_mi_eq_fv(spinor1, phi_);

	      U_ = gauge_field + ggi(get_index(it, ix, iy, iz-1, T, L), 3);

	      fv_pl_eq_cm_dag_ti_fv(xi_, U_, spinor1);


	      // Positive z-direction.

	      phi_ = phi + gsi(get_index(it, ix, iy, iz+1, T, L));

	      fv_eq_gamma_ti_fv(spinor1, 3, phi_);

	      fv_mi_eq_fv(spinor1, phi_);

	      U_ = gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 3);

	      fv_pl_eq_cm_ti_fv(xi_, U_, spinor1);


	      // Multiplication with 1/2.

	      fv_ti_eq_re(xi_, 0.5);


	      // Diagonal elements.

	      phi_ = phi + index_s;
		    
	      fv_eq_fv_ti_re(spinor1, phi_, _1_2_kappa);
	      fv_pl_eq_fv(xi_, spinor1);
		    
	      fv_eq_gamma_ti_fv(spinor1, 5, phi_);
	      fv_eq_fv_ti_im(spinor2, spinor1, mu);
	      fv_pl_eq_fv(xi_, spinor2);
	    }
	}
    }
}


void B_phi(double * xi, double * phi, const int T, const int L, 
	   const double kappa, const double mu, bool dagger) {

  int it, ix, iy, iz;
  double spinor1[24], spinor2[24];
  double nrm = 1./(1. + 4.*kappa*kappa*mu*mu);

  for(it = 0; it < T; it++) {
    fprintf(stderr, "void B_phi(...   -->   it = %2d ...\n", it);

    for(ix = 0; ix < L; ix++) {
      for(iy = 0; iy < L; iy++) {
	for(iz = 0; iz < L; iz++) {
	  int index_s = gsi(get_index(it, ix, iy, iz, T, L));

	  double *xi_ = xi + index_s;
	  double *phi_;

	  phi_ = phi + index_s;
		    
	  fv_eq_fv_ti_re(xi_, phi_, nrm);
	  
	  fv_eq_gamma_ti_fv(spinor1, 5, phi_);
	  if(dagger) fv_eq_fv_ti_im(spinor2, spinor1, 2*kappa*mu*nrm);
	  else fv_eq_fv_ti_im(spinor2, spinor1, -2*kappa*mu*nrm);
	  fv_pl_eq_fv(xi_, spinor2);
	}
      }
    }
  }
  return;
}


// ********************



// Computes (xi_c,xi_s) = Q_h (phi_c,phi_s), where Q_h is the heavy tm Dirac
// operator.

void Q_h_phi(double *xi_c, double *xi_s, double *phi_c, double *phi_s, double *gauge_field, int T, int L, double kappa, double mu_sigma, double mu_delta, bool diag)
{
  c_Q_h_phi(xi_c, phi_c, phi_s, gauge_field, T, L, kappa, mu_sigma, mu_delta, diag);
  s_Q_h_phi(xi_s, phi_c, phi_s, gauge_field, T, L, kappa, mu_sigma, mu_delta, diag);
}



// ********************



// Computes (xi_c,...) = Q_h (phi_c,phi_s), where Q_h is the heavy tm Dirac
// operator.

void c_Q_h_phi(double *xi_c, double *phi_c, double *phi_s, double *gauge_field, int T, int L, double kappa, double mu_sigma, double mu_delta, bool diag)
{
  int it, ix, iy, iz;
  double SU3_1[18];
  double spinor1[24], spinor2[24];


  double _1_2_kappa = 0.5 / kappa;

  for(it = 0; it < T; it++)
    {
      fprintf(stderr, "void Q_h_phi(...   -->   it = %2d ...\n", it);

      double angle_B = M_PI / (double)T;

      complex B;
      B.re = cos(angle_B);
      B.im = sin(angle_B);

      for(ix = 0; ix < L; ix++)
	{
	  for(iy = 0; iy < L; iy++)
	    {
	      for(iz = 0; iz < L; iz++)
		{
		  int index_s = gsi(get_index(it, ix, iy, iz, T, L));

		  double *xi_c_ = xi_c + index_s;
		  double *phi_c_, *phi_s_;
		  double *U_;


		  fv_eq_zero(xi_c_);


		  // Negative t-direction.

		  phi_c_ = phi_c + gsi(get_index(it-1, ix, iy, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 0, phi_c_);
		  fv_mi(spinor1);

		  fv_mi_eq_fv(spinor1, phi_c_);

		  U_ = gauge_field + ggi(get_index(it-1, ix, iy, iz, T, L), 0);

		  cm_eq_cm_ti_co(SU3_1, U_, &B);
		  fv_pl_eq_cm_dag_ti_fv(xi_c_, SU3_1, spinor1);


		  // Positive t-direction.

		  phi_c_ = phi_c + gsi(get_index(it+1, ix, iy, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 0, phi_c_);

		  fv_mi_eq_fv(spinor1, phi_c_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 0);

		  cm_eq_cm_ti_co(SU3_1, U_, &B);
		  fv_pl_eq_cm_ti_fv(xi_c_, SU3_1, spinor1);


		  // Negative x-direction.

		  phi_c_ = phi_c + gsi(get_index(it, ix-1, iy, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 1, phi_c_);
		  fv_mi(spinor1);

		  fv_mi_eq_fv(spinor1, phi_c_);

		  U_ = gauge_field + ggi(get_index(it, ix-1, iy, iz, T, L), 1);

		  fv_pl_eq_cm_dag_ti_fv(xi_c_, U_, spinor1);


		  // Positive x-direction.

		  phi_c_ = phi_c + gsi(get_index(it, ix+1, iy, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 1, phi_c_);

		  fv_mi_eq_fv(spinor1, phi_c_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 1);

		  fv_pl_eq_cm_ti_fv(xi_c_, U_, spinor1);


		  // Negative y-direction.

		  phi_c_ = phi_c + gsi(get_index(it, ix, iy-1, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 2, phi_c_);
		  fv_mi(spinor1);

		  fv_mi_eq_fv(spinor1, phi_c_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy-1, iz, T, L), 2);

		  fv_pl_eq_cm_dag_ti_fv(xi_c_, U_, spinor1);


		  // Positive y-direction.

		  phi_c_ = phi_c + gsi(get_index(it, ix, iy+1, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 2, phi_c_);

		  fv_mi_eq_fv(spinor1, phi_c_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 2);

		  fv_pl_eq_cm_ti_fv(xi_c_, U_, spinor1);


		  // Negative z-direction.

		  phi_c_ = phi_c + gsi(get_index(it, ix, iy, iz-1, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 3, phi_c_);
		  fv_mi(spinor1);

		  fv_mi_eq_fv(spinor1, phi_c_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy, iz-1, T, L), 3);

		  fv_pl_eq_cm_dag_ti_fv(xi_c_, U_, spinor1);


		  // Positive z-direction.

		  phi_c_ = phi_c + gsi(get_index(it, ix, iy, iz+1, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 3, phi_c_);

		  fv_mi_eq_fv(spinor1, phi_c_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 3);

		  fv_pl_eq_cm_ti_fv(xi_c_, U_, spinor1);


		  if(diag) {
		    // Multiplication with 1/2.
		    fv_ti_eq_re(xi_c_, 0.5);
		    
		    // Diagonal elements.
		    phi_c_ = phi_c + index_s;
		    phi_s_ = phi_s + index_s;
		    
		    fv_eq_fv_ti_re(spinor1, phi_c_, _1_2_kappa);
		    fv_pl_eq_fv(xi_c_, spinor1);
		    
		    fv_eq_gamma_ti_fv(spinor1, 5, phi_s_);
		    fv_eq_fv_ti_im(spinor2, spinor1, mu_sigma);
		    fv_pl_eq_fv(xi_c_, spinor2);
		    
		    fv_eq_fv_ti_re(spinor1, phi_c_, mu_delta);
		    fv_pl_eq_fv(xi_c_, spinor1);
		  }
		  else {
		    fv_eq_fv_ti_re(xi_c_, xi_c_, kappa);
		  }
		}
	    }
	}
    }
}


void c_B_h_phi(double *xi_c, double *phi_c, double *phi_s, const int T, const int L, 
	       const double kappa, const double mu_sigma, const double mu_delta, 
	       bool dagger) { 

  int it, ix, iy, iz;
  double spinor1[24];
  double nrm = 1./(1+4*kappa*kappa*(mu_sigma*mu_sigma-mu_delta*mu_delta));

  for(it = 0; it < T; it++) {
    fprintf(stderr, "void c_B_h_phi(...   -->   it = %2d ...\n", it);

    for(ix = 0; ix < L; ix++) {
      for(iy = 0; iy < L; iy++) {
	for(iz = 0; iz < L; iz++) {
	  int index_s = gsi(get_index(it, ix, iy, iz, T, L));

	  double *xi_c_ = xi_c + index_s;
	  double *phi_c_, *phi_s_;
	  
	  phi_s_ = phi_s + index_s;
	  phi_c_ = phi_c + index_s;
	  
	  fv_eq_gamma_ti_fv(spinor1, 5, phi_s_);
	  if(dagger) fv_eq_fv_ti_im(xi_c_, spinor1, 2*kappa*mu_sigma*nrm);
	  else fv_eq_fv_ti_im(xi_c_, spinor1, -2*kappa*mu_sigma*nrm);
	  
	  fv_eq_fv_ti_re(spinor1, phi_c_, (1-2*kappa*mu_delta)*nrm);
	  fv_pl_eq_fv(xi_c_, spinor1);
	}
      }
    }
  }
  return;
}



// ********************



// Computes (...,xi_s) = Q_h (phi_c,phi_s), where Q_h is the heavy tm Dirac
// operator.

void s_Q_h_phi(double *xi_s, double *phi_c, double *phi_s, double *gauge_field, int T, int L, double kappa, double mu_sigma, double mu_delta, bool diag)
{
  int it, ix, iy, iz;
  double SU3_1[18];
  double spinor1[24], spinor2[24];


  double _1_2_kappa = 0.5 / kappa;

  for(it = 0; it < T; it++)
    {
      fprintf(stderr, "void Q_h_phi(...   -->   it = %2d ...\n", it);

      double angle_B = M_PI / (double)T;

      complex B;
      B.re = cos(angle_B);
      B.im = sin(angle_B);

      for(ix = 0; ix < L; ix++)
	{
	  for(iy = 0; iy < L; iy++)
	    {
	      for(iz = 0; iz < L; iz++)
		{
		  int index_s = gsi(get_index(it, ix, iy, iz, T, L));

		  double *xi_s_ = xi_s + index_s;
		  double *phi_c_, *phi_s_;
		  double *U_;


		  fv_eq_zero(xi_s_);


		  // Negative t-direction.

		  phi_s_ = phi_s + gsi(get_index(it-1, ix, iy, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 0, phi_s_);
		  fv_mi(spinor1);

		  fv_mi_eq_fv(spinor1, phi_s_);

		  U_ = gauge_field + ggi(get_index(it-1, ix, iy, iz, T, L), 0);

		  cm_eq_cm_ti_co(SU3_1, U_, &B);
		  fv_pl_eq_cm_dag_ti_fv(xi_s_, SU3_1, spinor1);


		  // Positive t-direction.

		  phi_s_ = phi_s + gsi(get_index(it+1, ix, iy, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 0, phi_s_);

		  fv_mi_eq_fv(spinor1, phi_s_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 0);

		  cm_eq_cm_ti_co(SU3_1, U_, &B);
		  fv_pl_eq_cm_ti_fv(xi_s_, SU3_1, spinor1);


		  // Negative x-direction.

		  phi_s_ = phi_s + gsi(get_index(it, ix-1, iy, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 1, phi_s_);
		  fv_mi(spinor1);

		  fv_mi_eq_fv(spinor1, phi_s_);

		  U_ = gauge_field + ggi(get_index(it, ix-1, iy, iz, T, L), 1);

		  fv_pl_eq_cm_dag_ti_fv(xi_s_, U_, spinor1);


		  // Positive x-direction.

		  phi_s_ = phi_s + gsi(get_index(it, ix+1, iy, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 1, phi_s_);

		  fv_mi_eq_fv(spinor1, phi_s_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 1);

		  fv_pl_eq_cm_ti_fv(xi_s_, U_, spinor1);


		  // Negative y-direction.

		  phi_s_ = phi_s + gsi(get_index(it, ix, iy-1, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 2, phi_s_);
		  fv_mi(spinor1);

		  fv_mi_eq_fv(spinor1, phi_s_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy-1, iz, T, L), 2);

		  fv_pl_eq_cm_dag_ti_fv(xi_s_, U_, spinor1);


		  // Positive y-direction.

		  phi_s_ = phi_s + gsi(get_index(it, ix, iy+1, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 2, phi_s_);

		  fv_mi_eq_fv(spinor1, phi_s_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 2);

		  fv_pl_eq_cm_ti_fv(xi_s_, U_, spinor1);


		  // Negative z-direction.

		  phi_s_ = phi_s + gsi(get_index(it, ix, iy, iz-1, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 3, phi_s_);
		  fv_mi(spinor1);

		  fv_mi_eq_fv(spinor1, phi_s_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy, iz-1, T, L), 3);

		  fv_pl_eq_cm_dag_ti_fv(xi_s_, U_, spinor1);


		  // Positive z-direction.

		  phi_s_ = phi_s + gsi(get_index(it, ix, iy, iz+1, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 3, phi_s_);

		  fv_mi_eq_fv(spinor1, phi_s_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 3);

		  fv_pl_eq_cm_ti_fv(xi_s_, U_, spinor1);


		  if(diag) {
		    // Multiplication with 1/2.
		    fv_ti_eq_re(xi_s_, 0.5);

		    // Diagonal elements.
		    phi_s_ = phi_s + index_s;
		    phi_c_ = phi_c + index_s;
		    
		    fv_eq_fv_ti_re(spinor1, phi_s_, _1_2_kappa);
		    fv_pl_eq_fv(xi_s_, spinor1);
		    
		    fv_eq_gamma_ti_fv(spinor1, 5, phi_c_);
		    fv_eq_fv_ti_im(spinor2, spinor1, mu_sigma);
		    fv_pl_eq_fv(xi_s_, spinor2);
		    
		    fv_eq_fv_ti_re(spinor1, phi_s_, mu_delta);
		    fv_mi_eq_fv(xi_s_, spinor1);
		  }
		  else {
		    fv_eq_fv_ti_re(xi_s_, xi_s_, kappa);
		  }
		}
	    }
	}
    }
}



void s_B_h_phi(double *xi_s, double *phi_c, double *phi_s, const int T, const int L, 
	       const double kappa, const double mu_sigma, const double mu_delta, 
	       bool dagger) { 

  int it, ix, iy, iz;
  double spinor1[24];
  double nrm = 1./(1.+4.*kappa*kappa*(mu_sigma*mu_sigma-mu_delta*mu_delta));

  for(it = 0; it < T; it++) {
    fprintf(stderr, "void s_B_h_phi(...   -->   it = %2d ...\n", it);

    for(ix = 0; ix < L; ix++) {
      for(iy = 0; iy < L; iy++) {
	for(iz = 0; iz < L; iz++) {
	  int index_s = gsi(get_index(it, ix, iy, iz, T, L));

	  double *xi_s_ = xi_s + index_s;
	  double *phi_c_, *phi_s_;
	  
	  phi_s_ = phi_s + index_s;
	  phi_c_ = phi_c + index_s;
	  
	  fv_eq_gamma_ti_fv(spinor1, 5, phi_c_);
	  if(dagger) fv_eq_fv_ti_im(xi_s_, spinor1, 2*kappa*mu_sigma*nrm);
	  else fv_eq_fv_ti_im(xi_s_, spinor1, -2*kappa*mu_sigma*nrm);
	  
	  fv_eq_fv_ti_re(spinor1, phi_s_, (1+2*kappa*mu_delta)*nrm);
	  fv_pl_eq_fv(xi_s_, spinor1);
	}
      }
    }
  }
  return;
}


void B_h_phi(double *xi_c, double *xi_s, double *phi_c, double *phi_s, const int T, const int L, const double kappa, const double mu_sigma, const double mu_delta, bool dagger)
{
  c_B_h_phi(xi_c, phi_c, phi_s, T, L, kappa, mu_sigma, mu_delta, dagger);
  s_B_h_phi(xi_s, phi_c, phi_s, T, L, kappa, mu_sigma, mu_delta, dagger);
}

// ********************



// as above, but for periodic boundary conditions in time direction

void Q_h_phi_periodic(double *xi_c, double *xi_s, double *phi_c, double *phi_s, double *gauge_field, int T, int L, double kappa, double mu_sigma, double mu_delta, bool diag)
{
  c_Q_h_phi_periodic(xi_c, phi_c, phi_s, gauge_field, T, L, kappa, mu_sigma, mu_delta, diag);
  s_Q_h_phi_periodic(xi_s, phi_c, phi_s, gauge_field, T, L, kappa, mu_sigma, mu_delta, diag);
}

void c_Q_h_phi_periodic(double *xi_c, double *phi_c, double *phi_s, double *gauge_field, int T, int L, double kappa, double mu_sigma, double mu_delta, bool diag)
{
  int it, ix, iy, iz;
  double SU3_1[18];
  double spinor1[24], spinor2[24];


  double _1_2_kappa = 0.5 / kappa;

  for(it = 0; it < T; it++)
    {
      fprintf(stderr, "void Q_h_phi(...   -->   it = %2d ...\n", it);

      double angle_B = 0.0;

      complex B;
      B.re = cos(angle_B);
      B.im = sin(angle_B);

      for(ix = 0; ix < L; ix++)
	{
	  for(iy = 0; iy < L; iy++)
	    {
	      for(iz = 0; iz < L; iz++)
		{
		  int index_s = gsi(get_index(it, ix, iy, iz, T, L));

		  double *xi_c_ = xi_c + index_s;
		  double *phi_c_, *phi_s_;
		  double *U_;


		  fv_eq_zero(xi_c_);


		  // Negative t-direction.

		  phi_c_ = phi_c + gsi(get_index(it-1, ix, iy, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 0, phi_c_);
		  fv_mi(spinor1);

		  fv_mi_eq_fv(spinor1, phi_c_);

		  U_ = gauge_field + ggi(get_index(it-1, ix, iy, iz, T, L), 0);

		  cm_eq_cm_ti_co(SU3_1, U_, &B);
		  fv_pl_eq_cm_dag_ti_fv(xi_c_, SU3_1, spinor1);


		  // Positive t-direction.

		  phi_c_ = phi_c + gsi(get_index(it+1, ix, iy, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 0, phi_c_);

		  fv_mi_eq_fv(spinor1, phi_c_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 0);

		  cm_eq_cm_ti_co(SU3_1, U_, &B);
		  fv_pl_eq_cm_ti_fv(xi_c_, SU3_1, spinor1);


		  // Negative x-direction.

		  phi_c_ = phi_c + gsi(get_index(it, ix-1, iy, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 1, phi_c_);
		  fv_mi(spinor1);

		  fv_mi_eq_fv(spinor1, phi_c_);

		  U_ = gauge_field + ggi(get_index(it, ix-1, iy, iz, T, L), 1);

		  fv_pl_eq_cm_dag_ti_fv(xi_c_, U_, spinor1);


		  // Positive x-direction.

		  phi_c_ = phi_c + gsi(get_index(it, ix+1, iy, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 1, phi_c_);

		  fv_mi_eq_fv(spinor1, phi_c_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 1);

		  fv_pl_eq_cm_ti_fv(xi_c_, U_, spinor1);


		  // Negative y-direction.

		  phi_c_ = phi_c + gsi(get_index(it, ix, iy-1, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 2, phi_c_);
		  fv_mi(spinor1);

		  fv_mi_eq_fv(spinor1, phi_c_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy-1, iz, T, L), 2);

		  fv_pl_eq_cm_dag_ti_fv(xi_c_, U_, spinor1);


		  // Positive y-direction.

		  phi_c_ = phi_c + gsi(get_index(it, ix, iy+1, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 2, phi_c_);

		  fv_mi_eq_fv(spinor1, phi_c_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 2);

		  fv_pl_eq_cm_ti_fv(xi_c_, U_, spinor1);


		  // Negative z-direction.

		  phi_c_ = phi_c + gsi(get_index(it, ix, iy, iz-1, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 3, phi_c_);
		  fv_mi(spinor1);

		  fv_mi_eq_fv(spinor1, phi_c_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy, iz-1, T, L), 3);

		  fv_pl_eq_cm_dag_ti_fv(xi_c_, U_, spinor1);


		  // Positive z-direction.

		  phi_c_ = phi_c + gsi(get_index(it, ix, iy, iz+1, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 3, phi_c_);

		  fv_mi_eq_fv(spinor1, phi_c_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 3);

		  fv_pl_eq_cm_ti_fv(xi_c_, U_, spinor1);


		  if(diag) {
		    // Multiplication with 1/2.
		    fv_ti_eq_re(xi_c_, 0.5);
		    
		    // Diagonal elements.
		    phi_c_ = phi_c + index_s;
		    phi_s_ = phi_s + index_s;
		    
		    fv_eq_fv_ti_re(spinor1, phi_c_, _1_2_kappa);
		    fv_pl_eq_fv(xi_c_, spinor1);
		    
		    fv_eq_gamma_ti_fv(spinor1, 5, phi_s_);
		    fv_eq_fv_ti_im(spinor2, spinor1, mu_sigma);
		    fv_pl_eq_fv(xi_c_, spinor2);
		    
		    fv_eq_fv_ti_re(spinor1, phi_c_, mu_delta);
		    fv_pl_eq_fv(xi_c_, spinor1);
		  }
		  else {
		    fv_eq_fv_ti_re(xi_c_, xi_c_, kappa);
		  }
		}
	    }
	}
    }
}

void s_Q_h_phi_periodic(double *xi_s, double *phi_c, double *phi_s, double *gauge_field, int T, int L, double kappa, double mu_sigma, double mu_delta, bool diag)
{
  int it, ix, iy, iz;
  double SU3_1[18];
  double spinor1[24], spinor2[24];


  double _1_2_kappa = 0.5 / kappa;

  for(it = 0; it < T; it++)
    {
      fprintf(stderr, "void Q_h_phi(...   -->   it = %2d ...\n", it);

      double angle_B = 0.0;

      complex B;
      B.re = cos(angle_B);
      B.im = sin(angle_B);

      for(ix = 0; ix < L; ix++)
	{
	  for(iy = 0; iy < L; iy++)
	    {
	      for(iz = 0; iz < L; iz++)
		{
		  int index_s = gsi(get_index(it, ix, iy, iz, T, L));

		  double *xi_s_ = xi_s + index_s;
		  double *phi_c_, *phi_s_;
		  double *U_;


		  fv_eq_zero(xi_s_);


		  // Negative t-direction.

		  phi_s_ = phi_s + gsi(get_index(it-1, ix, iy, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 0, phi_s_);
		  fv_mi(spinor1);

		  fv_mi_eq_fv(spinor1, phi_s_);

		  U_ = gauge_field + ggi(get_index(it-1, ix, iy, iz, T, L), 0);

		  cm_eq_cm_ti_co(SU3_1, U_, &B);
		  fv_pl_eq_cm_dag_ti_fv(xi_s_, SU3_1, spinor1);


		  // Positive t-direction.

		  phi_s_ = phi_s + gsi(get_index(it+1, ix, iy, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 0, phi_s_);

		  fv_mi_eq_fv(spinor1, phi_s_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 0);

		  cm_eq_cm_ti_co(SU3_1, U_, &B);
		  fv_pl_eq_cm_ti_fv(xi_s_, SU3_1, spinor1);


		  // Negative x-direction.

		  phi_s_ = phi_s + gsi(get_index(it, ix-1, iy, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 1, phi_s_);
		  fv_mi(spinor1);

		  fv_mi_eq_fv(spinor1, phi_s_);

		  U_ = gauge_field + ggi(get_index(it, ix-1, iy, iz, T, L), 1);

		  fv_pl_eq_cm_dag_ti_fv(xi_s_, U_, spinor1);


		  // Positive x-direction.

		  phi_s_ = phi_s + gsi(get_index(it, ix+1, iy, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 1, phi_s_);

		  fv_mi_eq_fv(spinor1, phi_s_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 1);

		  fv_pl_eq_cm_ti_fv(xi_s_, U_, spinor1);


		  // Negative y-direction.

		  phi_s_ = phi_s + gsi(get_index(it, ix, iy-1, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 2, phi_s_);
		  fv_mi(spinor1);

		  fv_mi_eq_fv(spinor1, phi_s_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy-1, iz, T, L), 2);

		  fv_pl_eq_cm_dag_ti_fv(xi_s_, U_, spinor1);


		  // Positive y-direction.

		  phi_s_ = phi_s + gsi(get_index(it, ix, iy+1, iz, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 2, phi_s_);

		  fv_mi_eq_fv(spinor1, phi_s_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 2);

		  fv_pl_eq_cm_ti_fv(xi_s_, U_, spinor1);


		  // Negative z-direction.

		  phi_s_ = phi_s + gsi(get_index(it, ix, iy, iz-1, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 3, phi_s_);
		  fv_mi(spinor1);

		  fv_mi_eq_fv(spinor1, phi_s_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy, iz-1, T, L), 3);

		  fv_pl_eq_cm_dag_ti_fv(xi_s_, U_, spinor1);


		  // Positive z-direction.

		  phi_s_ = phi_s + gsi(get_index(it, ix, iy, iz+1, T, L));

		  fv_eq_gamma_ti_fv(spinor1, 3, phi_s_);

		  fv_mi_eq_fv(spinor1, phi_s_);

		  U_ = gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 3);

		  fv_pl_eq_cm_ti_fv(xi_s_, U_, spinor1);


		  if(diag) {
		    // Multiplication with 1/2.
		    fv_ti_eq_re(xi_s_, 0.5);

		    // Diagonal elements.
		    phi_s_ = phi_s + index_s;
		    phi_c_ = phi_c + index_s;
		    
		    fv_eq_fv_ti_re(spinor1, phi_s_, _1_2_kappa);
		    fv_pl_eq_fv(xi_s_, spinor1);
		    
		    fv_eq_gamma_ti_fv(spinor1, 5, phi_c_);
		    fv_eq_fv_ti_im(spinor2, spinor1, mu_sigma);
		    fv_pl_eq_fv(xi_s_, spinor2);
		    
		    fv_eq_fv_ti_re(spinor1, phi_s_, mu_delta);
		    fv_mi_eq_fv(xi_s_, spinor1);
		  }
		  else {
		    fv_eq_fv_ti_re(xi_s_, xi_s_, kappa);
		  }
		}
	    }
	}
    }
}



// ********************



