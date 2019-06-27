#include <cstdlib>

#include "geometry.hh"
#include "smearing_techniques.hh"
#include "fields.hh"
#include "fuzz.hh"

// ********************

// written by Chris Michael

// Combines Nlong APE-smeared links to one straight fuzzed link for 1 timeslice

// as written assumes Nlong > 1  
// if Nlong==1 -- needs to copy timeslice of smeared to fuzzed


void fuzzed_links_Timeslice(double *fuzzed_gauge_field, double *smeared_gauge_field, 
			    const unsigned int T, const unsigned int L, 
			    const unsigned int Nlong, const unsigned int timeslice) {
  unsigned long int index, index_;
  unsigned long int index_px , index_py , index_pz ;

  double * fuzzed_gauge_field_old;
  
  Timeslice_Gauge_Field_Alloc(&fuzzed_gauge_field_old, L);
  Timeslice_Gauge_Field_Copy(fuzzed_gauge_field_old, smeared_gauge_field, T, L,
			     timeslice);
  
  for(unsigned int ir = 1; ir < Nlong; ir++) {
    
    for(unsigned int ix = 0; ix < L; ix++) {
      for(unsigned int iy = 0; iy < L; iy++) {
	for(unsigned int iz = 0; iz < L; iz++) {
	  
	  
	  index_ = ggi(get_index(0, ix, iy, iz, T, L), 1);
	  index  = ggi(get_index(timeslice, ix, iy, iz, T, L), 1);
	  index_px  = ggi(get_index(timeslice, ix+ir, iy, iz, T, L), 1);
	  
	  cm_eq_cm_ti_cm(fuzzed_gauge_field+index, fuzzed_gauge_field_old + index_,
			 smeared_gauge_field + index_px );
	  
	  index_ = ggi(get_index(0, ix, iy, iz, T, L), 2);
	  index  = ggi(get_index(timeslice, ix, iy, iz, T, L), 2);
	  index_py  = ggi(get_index(timeslice, ix, iy+ir, iz, T, L), 2);
	  
	  cm_eq_cm_ti_cm(fuzzed_gauge_field+index, fuzzed_gauge_field_old + index_,
			 smeared_gauge_field + index_py );
	  
	  index_ = ggi(get_index(0, ix, iy, iz, T, L), 3);
	  index  = ggi(get_index(timeslice, ix, iy, iz, T, L), 3);
	  index_pz  = ggi(get_index(timeslice, ix, iy, iz+ir, T, L), 3);
	  
	  cm_eq_cm_ti_cm(fuzzed_gauge_field+index, fuzzed_gauge_field_old + index_,
			 smeared_gauge_field + index_pz );
	  
	}
      }
    }
    if (ir < (Nlong-1)) {
      Timeslice_Gauge_Field_Copy(fuzzed_gauge_field_old, fuzzed_gauge_field, T, L,
				 timeslice);
    }
  }
  
  Gauge_Field_Free(&fuzzed_gauge_field_old);
}



// ********************

// written by Chris Michael - based on Jacobi smearing code

// Creates a fuzzed spinor - overwrites input spinor

//  fuzzed gauge field must have length Nlong

// psi = quark spinor
// Nlong = Fuzzing  parameter
// timeslice = the timeslice, on which the fuzzing is performed

void Fuzz_prop(double *fuzzed_gauge_field, double *psi, 
	       const unsigned int T, const unsigned int L, const unsigned int Nlong, 
	       const double kappa, const unsigned int timeslice) {

  double *psi_old;
  Timeslice_Spinor_Field_Alloc(&psi_old, L);
  
  // Copy the timeslice of interest to psi_old.
  Timeslice_Spinor_Field_Copy(psi_old, psi, T, L, timeslice);
  
  for(unsigned int ix = 0; ix < L; ix++) {
    for(unsigned int iy = 0; iy < L; iy++) {
      for(unsigned int iz = 0; iz < L; iz++) {
	
	// Get indices.

	int index_s = gsi(get_index_timeslice(ix, iy, iz, T, L));
	int index_s_mx = gsi(get_index_timeslice(ix-Nlong, iy, iz, T, L));
	int index_s_px = gsi(get_index_timeslice(ix+Nlong, iy, iz, T, L));
	int index_s_my = gsi(get_index_timeslice(ix, iy-Nlong, iz, T, L));
	int index_s_py = gsi(get_index_timeslice(ix, iy+Nlong, iz, T, L));
	int index_s_mz = gsi(get_index_timeslice(ix, iy, iz-Nlong, T, L));
	int index_s_pz = gsi(get_index_timeslice(ix, iy, iz+Nlong, T, L));
	
	int index_g_mx = ggi(get_index(timeslice, ix-Nlong, iy, iz, T, L), 1);
	int index_g_px = ggi(get_index(timeslice, ix, iy, iz, T, L), 1);
	int index_g_my = ggi(get_index(timeslice, ix, iy-Nlong, iz, T, L), 2);
	int index_g_py = ggi(get_index(timeslice, ix, iy, iz, T, L), 2);
	int index_g_mz = ggi(get_index(timeslice, ix, iy, iz-Nlong, T, L), 3);
	int index_g_pz = ggi(get_index(timeslice, ix, iy, iz, T, L), 3);
	

	double *s = psi + gsi(get_index(timeslice, ix, iy, iz, T, L));
	
	fv_eq_zero(s);


	// negative x-direction
	
	fv_pl_eq_cm_dag_ti_fv(s, fuzzed_gauge_field + index_g_mx,
			      psi_old + index_s_mx);
	
	// positive x-direction
	
	fv_pl_eq_cm_ti_fv(s, fuzzed_gauge_field + index_g_px,
			  psi_old + index_s_px);
	
	// negative y-direction
	
	fv_pl_eq_cm_dag_ti_fv(s, fuzzed_gauge_field + index_g_my,
			      psi_old + index_s_my);
	
	// positive y-direction
	
	fv_pl_eq_cm_ti_fv(s, fuzzed_gauge_field + index_g_py,
			  psi_old + index_s_py);
	
	// negative z-direction
	
	fv_pl_eq_cm_dag_ti_fv(s, fuzzed_gauge_field + index_g_mz,
			      psi_old + index_s_mz);
	
	// positive z-direction
	
	fv_pl_eq_cm_ti_fv(s, fuzzed_gauge_field + index_g_pz,
			  psi_old + index_s_pz);
	
	
	
	// conventionally CM did not normalise here though /6.0 makes numbers nicer 
	// double norm = 1.0 / 6.0;
	// fv_ti_eq_re(s, norm);
	
      }
    }
  }      
  Spinor_Field_Free(&psi_old);
}



// ********************

