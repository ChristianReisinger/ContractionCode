// ********************

// smearing_techniques.hh

// Author: Marc Wagner
// Date: September 2007

// ********************

#ifndef __SMEARING_TECHNIQUES_HH__

#define __SMEARING_TECHNIQUES_HH__

// ********************

#include "linear_algebra.hh"

// ********************

// Computes fat links in time direction.

void Fat_Time_Links(double *gauge_field, double *smeared_gauge_field, int T, int L, double time_link_epsilon);

// Computes HYP links in time direction.

void HYP_Time_Links(double *gauge_field, int T, int L, double time_link_alpha1, double time_link_alpha2, double time_link_alpha3);

// Performs an APE smearing step.

void APE_Smearing_Step(double *smeared_gauge_field, int T, int L, double APE_smearing_alpha);

// Performs an APE smearing step on a given timeslice.

void APE_Smearing_Step_Timeslice(double *smeared_gauge_field, int T, int L, double APE_smearing_alpha, int timeslice);

// Performs a number of Jacobi smearing steps on a given timeslice.

// psi = quark spinor
// N, kappa = Jacobi smearing parameters
// timeslice = the timeslice, on which the smearing is performed

void Jacobi_Smearing_Steps(double *smeared_gauge_field, double *psi, int T, int L, int N, double kappa, int timeslice,
		bool use_slice_0_for_psi = false);  // for __OPTIMIZE_MEM__

// ********************

#endif

// ********************
