// ********************



// Q_phi.hh

// Author: Marc Wagner
// Date: September 2007



// ********************



#ifndef __Q_PHI_HH__

#define __Q_PHI_HH__



// ********************



// Computes xi = Q phi, where Q is the light tm Dirac operator.
//
// if diag=false only the hopping matrix is applied
// the result is then in the hopping parameter normalisation!!

void Q_phi(double *xi, double *phi, double *gauge_field, int T, int L, double kappa, double mu, bool diag = true);

void Q_phi_timeslice(double *xi, double *phi, double *gauge_field, int T, int L, double kappa, double mu, int timeslice, int T_B);


// computes xi = B phi
// B = (1+i mu gamma_5)/(1+mu^2)
//
// result is in the hopping parameter normalisation!!

void B_phi(double * xi, double * phi, const int T, const int L, 
	   const double kappa, const double mu, bool dagger);

// Computes (xi_c,xi_s) = Q_h (phi_c,phi_s), where Q_h is the heavy tm Dirac
// operator.
// if diag=false only the hopping matrix is applied
// the result is then in the hopping parameter normalisation!!

void Q_h_phi(double *xi_c, double *xi_s, double *phi_c, double *phi_s, double *gauge_field, int T, int L, double kappa, double mu_sigma, double mu_delta, bool diag = true);

// Computes (xi_c,...) = Q_h (phi_c,phi_s), where Q_h is the heavy tm Dirac
// operator.

void c_Q_h_phi(double *xi_c, double *phi_c, double *phi_s, double *gauge_field, int T, int L, double kappa, double mu_sigma, double mu_delta, bool diag = true);

// Computes (...,xi_s) = Q_h (phi_c,phi_s), where Q_h is the heavy tm Dirac
// operator.

void s_Q_h_phi(double *xi_s, double *phi_c, double *phi_s, double *gauge_field, int T, int L, double kappa, double mu_sigma, double mu_delta, bool diag = true);

// computes (xi_c,xi_s) = B_h*(phi_c, phi_s)
// B_h = (1-i\g5\tau^1\musigma-\tau^3\mudelta)/c
// where
// c = 1 + \musigma^2 + \mudelta^2
//
// result is in the hopping parameter normalisation!!

void B_h_phi(double *xi_c, double *xi_s, double *phi_c, double *phi_s, 
	     const int T, const int L, const double kappa, 
	     const double mu_sigma, const double mu_delta, bool dagger);



// *****



// as above, but for periodic boundary conditions in time direction

void Q_h_phi_periodic(double *xi_c, double *xi_s, double *phi_c, double *phi_s, double *gauge_field, int T, int L, double kappa, double mu_sigma, double mu_delta, bool diag = true);

void c_Q_h_phi_periodic(double *xi_c, double *phi_c, double *phi_s, double *gauge_field, int T, int L, double kappa, double mu_sigma, double mu_delta, bool diag = true);

void s_Q_h_phi_periodic(double *xi_s, double *phi_c, double *phi_s, double *gauge_field, int T, int L, double kappa, double mu_sigma, double mu_delta, bool diag = true);



// ********************



#endif



// ********************
