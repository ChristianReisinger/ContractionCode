#ifndef _FUZZ_HH
#define _FUZZ_HH


void fuzzed_links_Timeslice(double *fuzzed_gauge_field, double *smeared_gauge_field, 
			    const unsigned int T, const unsigned int L, 
			    const unsigned int Nlong, const unsigned int timeslice);

void Fuzz_prop(double *fuzzed_gauge_field, double *psi, 
	       const unsigned int T, const unsigned int L, const unsigned int Nlong, 
	       const double kappa, const unsigned int timeslice);


#endif
