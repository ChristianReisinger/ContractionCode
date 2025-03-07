// ********************

// smearing_techniques.cc

// Author: Marc Wagner
// Date: September 2007

// ********************

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "fields.hh"
#include "geometry.hh"

#include "smearing_techniques.hh"

// ********************

// Computes fat links in time direction.

void Fat_Time_Links(double *gauge_field, double *smeared_gauge_field, int T, int L, double time_link_epsilon) {
	int it, ix, iy, iz;
	double M1[18], M2[18];

	for (it = 0; it < T; it++) {
		for (ix = 0; ix < L; ix++) {
			for (iy = 0; iy < L; iy++) {
				for (iz = 0; iz < L; iz++) {
					// *****
					// *****
					// *****
					// *****
					// *****

					int index = ggi(get_index(it, ix, iy, iz, T, L), 0);

					int index_mx_1 = ggi(get_index(it, ix - 1, iy, iz, T, L), 1);
					int index_mx_2 = ggi(get_index(it, ix - 1, iy, iz, T, L), 0);
					int index_mx_3 = ggi(get_index(it + 1, ix - 1, iy, iz, T, L), 1);

					int index_px_1 = ggi(get_index(it, ix, iy, iz, T, L), 1);
					int index_px_2 = ggi(get_index(it, ix + 1, iy, iz, T, L), 0);
					int index_px_3 = ggi(get_index(it + 1, ix, iy, iz, T, L), 1);

					int index_my_1 = ggi(get_index(it, ix, iy - 1, iz, T, L), 2);
					int index_my_2 = ggi(get_index(it, ix, iy - 1, iz, T, L), 0);
					int index_my_3 = ggi(get_index(it + 1, ix, iy - 1, iz, T, L), 2);

					int index_py_1 = ggi(get_index(it, ix, iy, iz, T, L), 2);
					int index_py_2 = ggi(get_index(it, ix, iy + 1, iz, T, L), 0);
					int index_py_3 = ggi(get_index(it + 1, ix, iy, iz, T, L), 2);

					int index_mz_1 = ggi(get_index(it, ix, iy, iz - 1, T, L), 3);
					int index_mz_2 = ggi(get_index(it, ix, iy, iz - 1, T, L), 0);
					int index_mz_3 = ggi(get_index(it + 1, ix, iy, iz - 1, T, L), 3);

					int index_pz_1 = ggi(get_index(it, ix, iy, iz, T, L), 3);
					int index_pz_2 = ggi(get_index(it, ix, iy, iz + 1, T, L), 0);
					int index_pz_3 = ggi(get_index(it + 1, ix, iy, iz, T, L), 3);

					double *U = smeared_gauge_field + index;

					// center

					if (time_link_epsilon != 0.0) {
						cm_eq_cm_ti_re(U, gauge_field + index, time_link_epsilon);
					} else
						cm_eq_zero(U);

					// negative x-direction

					cm_eq_cm_ti_cm(M1, gauge_field + index_mx_2, gauge_field + index_mx_3);
					cm_eq_cm_dag_ti_cm(M2, gauge_field + index_mx_1, M1);
					cm_pl_eq_cm(U, M2);

					// positive x-direction

					cm_eq_cm_ti_cm_dag(M1, gauge_field + index_px_2, gauge_field + index_px_3);
					cm_eq_cm_ti_cm(M2, gauge_field + index_px_1, M1);
					cm_pl_eq_cm(U, M2);

					// negative y-direction

					cm_eq_cm_ti_cm(M1, gauge_field + index_my_2, gauge_field + index_my_3);
					cm_eq_cm_dag_ti_cm(M2, gauge_field + index_my_1, M1);
					cm_pl_eq_cm(U, M2);

					// positive y-direction

					cm_eq_cm_ti_cm_dag(M1, gauge_field + index_py_2, gauge_field + index_py_3);
					cm_eq_cm_ti_cm(M2, gauge_field + index_py_1, M1);
					cm_pl_eq_cm(U, M2);

					// negative z-direction

					cm_eq_cm_ti_cm(M1, gauge_field + index_mz_2, gauge_field + index_mz_3);
					cm_eq_cm_dag_ti_cm(M2, gauge_field + index_mz_1, M1);
					cm_pl_eq_cm(U, M2);

					// positive z-direction

					cm_eq_cm_ti_cm_dag(M1, gauge_field + index_pz_2, gauge_field + index_pz_3);
					cm_eq_cm_ti_cm(M2, gauge_field + index_pz_1, M1);
					cm_pl_eq_cm(U, M2);

					// Projection to SU(3).
					cm_proj(U);

					// *****
					// *****
					// *****
					// *****
					// *****
				}
			}
		}
	}
}

// ********************

// Computes HYP links in time direction.

// Decorated links (temporary).

// Va_b_cd   -->   a = level (2 or 3), b = direction, c,d = staples.

double *V2_0_12, *V2_0_13, *V2_0_23, *V2_1_23, *V2_2_13, *V2_3_12;
double *V3_0_1, *V3_0_2, *V3_0_3;
double *V3_1_2, *V3_1_3;
double *V3_2_1, *V3_2_3;
double *V3_3_1, *V3_3_2;

inline void HYP_Helper_3(double *gauge_field, int T, int L, double time_link_alpha3, double *U, int it, int ix, int iy, int iz,
		int mu, int nu, int rho) {
	int eta;
	int index_1, index_2, index_3;
	double M1[18], M2[18];

	cm_eq_zero(U);

	for (eta = 0; eta < 4; eta++)
	// For all staples.
			{
		if (eta == mu || eta == nu || eta == rho)
			continue;

		// negative eta-direction

		if (mu == 0) {
			if (eta == 1) {
				index_1 = ggi(get_index(it, ix - 1, iy, iz, T, L), eta);
				index_2 = ggi(get_index(it, ix - 1, iy, iz, T, L), mu);
				index_3 = ggi(get_index(it + 1, ix - 1, iy, iz, T, L), eta);
			} else if (eta == 2) {
				index_1 = ggi(get_index(it, ix, iy - 1, iz, T, L), eta);
				index_2 = ggi(get_index(it, ix, iy - 1, iz, T, L), mu);
				index_3 = ggi(get_index(it + 1, ix, iy - 1, iz, T, L), eta);
			} else if (eta == 3) {
				index_1 = ggi(get_index(it, ix, iy, iz - 1, T, L), eta);
				index_2 = ggi(get_index(it, ix, iy, iz - 1, T, L), mu);
				index_3 = ggi(get_index(it + 1, ix, iy, iz - 1, T, L), eta);
			} else {
				fprintf(stderr, "Error: void HYP_Helper_3(...\n");
				exit(EXIT_FAILURE);
			}
		} else if (mu == 1) {
			if (eta == 0) {
				fprintf(stderr, "Error: void HYP_Helper_3(...\n");
				exit(EXIT_FAILURE);

				/*
				 index_1 = ggi(get_index(it-1, ix  , iy  , iz  , T, L), eta);
				 index_2 = ggi(get_index(it-1, ix  , iy  , iz  , T, L), mu);
				 index_3 = ggi(get_index(it-1, ix+1, iy  , iz  , T, L), eta);
				 */
			} else if (eta == 2) {
				index_1 = ggi(get_index(it, ix, iy - 1, iz, T, L), eta);
				index_2 = ggi(get_index(it, ix, iy - 1, iz, T, L), mu);
				index_3 = ggi(get_index(it, ix + 1, iy - 1, iz, T, L), eta);
			} else if (eta == 3) {
				index_1 = ggi(get_index(it, ix, iy, iz - 1, T, L), eta);
				index_2 = ggi(get_index(it, ix, iy, iz - 1, T, L), mu);
				index_3 = ggi(get_index(it, ix + 1, iy, iz - 1, T, L), eta);
			} else {
				fprintf(stderr, "Error: void HYP_Helper_3(...\n");
				exit(EXIT_FAILURE);
			}
		} else if (mu == 2) {
			if (eta == 0) {
				fprintf(stderr, "Error: void HYP_Helper_3(...\n");
				exit(EXIT_FAILURE);

				/*
				 index_1 = ggi(get_index(it-1, ix  , iy  ,iz  , T, L), eta);
				 index_2 = ggi(get_index(it-1, ix  , iy  ,iz  , T, L), mu);
				 index_3 = ggi(get_index(it-1, ix  , iy+1,iz  , T, L), eta);
				 */
			} else if (eta == 1) {
				index_1 = ggi(get_index(it, ix - 1, iy, iz, T, L), eta);
				index_2 = ggi(get_index(it, ix - 1, iy, iz, T, L), mu);
				index_3 = ggi(get_index(it, ix - 1, iy + 1, iz, T, L), eta);
			} else if (eta == 3) {
				index_1 = ggi(get_index(it, ix, iy, iz - 1, T, L), eta);
				index_2 = ggi(get_index(it, ix, iy, iz - 1, T, L), mu);
				index_3 = ggi(get_index(it, ix, iy + 1, iz - 1, T, L), eta);
			} else {
				fprintf(stderr, "Error: void HYP_Helper_3(...\n");
				exit(EXIT_FAILURE);
			}
		} else if (mu == 3) {
			if (eta == 0) {
				fprintf(stderr, "Error: void HYP_Helper_3(...\n");
				exit(EXIT_FAILURE);

				/*
				 index_1 = ggi(get_index(it-1, ix  , iy  , iz  , T, L), eta);
				 index_2 = ggi(get_index(it-1, ix  , iy  , iz  , T, L), mu);
				 index_3 = ggi(get_index(it-1, ix  , iy  , iz+1, T, L), eta);
				 */
			} else if (eta == 1) {
				index_1 = ggi(get_index(it, ix - 1, iy, iz, T, L), eta);
				index_2 = ggi(get_index(it, ix - 1, iy, iz, T, L), mu);
				index_3 = ggi(get_index(it, ix - 1, iy, iz + 1, T, L), eta);
			} else if (eta == 2) {
				index_1 = ggi(get_index(it, ix, iy - 1, iz, T, L), eta);
				index_2 = ggi(get_index(it, ix, iy - 1, iz, T, L), mu);
				index_3 = ggi(get_index(it, ix, iy - 1, iz + 1, T, L), eta);
			} else {
				fprintf(stderr, "Error: void HYP_Helper_3(...\n");
				exit(EXIT_FAILURE);
			}
		} else {
			fprintf(stderr, "Error: void HYP_Helper_3(...\n");
			exit(EXIT_FAILURE);
		}

		cm_eq_cm_ti_cm(M1, gauge_field + index_2, gauge_field + index_3);
		cm_eq_cm_dag_ti_cm(M2, gauge_field + index_1, M1);

		cm_pl_eq_cm(U, M2);

		// positive eta-direction

		if (mu == 0) {
			if (eta == 1) {
				index_1 = ggi(get_index(it, ix, iy, iz, T, L), eta);
				index_2 = ggi(get_index(it, ix + 1, iy, iz, T, L), mu);
				index_3 = ggi(get_index(it + 1, ix, iy, iz, T, L), eta);
			} else if (eta == 2) {
				index_1 = ggi(get_index(it, ix, iy, iz, T, L), eta);
				index_2 = ggi(get_index(it, ix, iy + 1, iz, T, L), mu);
				index_3 = ggi(get_index(it + 1, ix, iy, iz, T, L), eta);
			} else if (eta == 3) {
				index_1 = ggi(get_index(it, ix, iy, iz, T, L), eta);
				index_2 = ggi(get_index(it, ix, iy, iz + 1, T, L), mu);
				index_3 = ggi(get_index(it + 1, ix, iy, iz, T, L), eta);
			} else {
				fprintf(stderr, "Error: void HYP_Helper_3(...\n");
				exit(EXIT_FAILURE);
			}
		} else if (mu == 1) {
			if (eta == 0) {
				fprintf(stderr, "Error: void HYP_Helper_3(...\n");
				exit(EXIT_FAILURE);

				/*
				 index_1 = ggi(get_index(it  , ix  , iy  , iz  , T, L), eta);
				 index_2 = ggi(get_index(it+1, ix  , iy  , iz  , T, L), mu);
				 index_3 = ggi(get_index(it  , ix+1, iy  , iz  , T, L), eta);
				 */
			} else if (eta == 2) {
				index_1 = ggi(get_index(it, ix, iy, iz, T, L), eta);
				index_2 = ggi(get_index(it, ix, iy + 1, iz, T, L), mu);
				index_3 = ggi(get_index(it, ix + 1, iy, iz, T, L), eta);
			} else if (eta == 3) {
				index_1 = ggi(get_index(it, ix, iy, iz, T, L), eta);
				index_2 = ggi(get_index(it, ix, iy, iz + 1, T, L), mu);
				index_3 = ggi(get_index(it, ix + 1, iy, iz, T, L), eta);
			} else {
				fprintf(stderr, "Error: void HYP_Helper_3(...\n");
				exit(EXIT_FAILURE);
			}
		} else if (mu == 2) {
			if (eta == 0) {
				fprintf(stderr, "Error: void HYP_Helper_3(...\n");
				exit(EXIT_FAILURE);

				/*
				 index_1 = ggi(get_index(it  , ix  , iy  , iz  , T, L), eta);
				 index_2 = ggi(get_index(it+1, ix  , iy  , iz  , T, L), mu);
				 index_3 = ggi(get_index(it  , ix  , iy+1, iz  , T, L), eta);
				 */
			} else if (eta == 1) {
				index_1 = ggi(get_index(it, ix, iy, iz, T, L), eta);
				index_2 = ggi(get_index(it, ix + 1, iy, iz, T, L), mu);
				index_3 = ggi(get_index(it, ix, iy + 1, iz, T, L), eta);
			} else if (eta == 3) {
				index_1 = ggi(get_index(it, ix, iy, iz, T, L), eta);
				index_2 = ggi(get_index(it, ix, iy, iz + 1, T, L), mu);
				index_3 = ggi(get_index(it, ix, iy + 1, iz, T, L), eta);
			} else {
				fprintf(stderr, "Error: void HYP_Helper_3(...\n");
				exit(EXIT_FAILURE);
			}
		} else if (mu == 3) {
			if (eta == 0) {
				fprintf(stderr, "Error: void HYP_Helper_3(...\n");
				exit(EXIT_FAILURE);

				/*
				 index_1 = ggi(get_index(it  , ix  , iy  , iz  , T, L), eta);
				 index_2 = ggi(get_index(it+1, ix  , iy  , iz  , T, L), mu);
				 index_3 = ggi(get_index(it  , ix  , iy  , iz+1, T, L), eta);
				 */
			} else if (eta == 1) {
				index_1 = ggi(get_index(it, ix, iy, iz, T, L), eta);
				index_2 = ggi(get_index(it, ix + 1, iy, iz, T, L), mu);
				index_3 = ggi(get_index(it, ix, iy, iz + 1, T, L), eta);
			} else if (eta == 2) {
				index_1 = ggi(get_index(it, ix, iy, iz, T, L), eta);
				index_2 = ggi(get_index(it, ix, iy + 1, iz, T, L), mu);
				index_3 = ggi(get_index(it, ix, iy, iz + 1, T, L), eta);
			} else {
				fprintf(stderr, "Error: void HYP_Helper_3(...\n");
				exit(EXIT_FAILURE);
			}
		} else {
			fprintf(stderr, "Error: void HYP_Helper_3(...\n");
			exit(EXIT_FAILURE);
		}

		cm_eq_cm_ti_cm_dag(M1, gauge_field + index_2, gauge_field + index_3);
		cm_eq_cm_ti_cm(M2, gauge_field + index_1, M1);

		cm_pl_eq_cm(U, M2);
	}

	cm_ti_eq_re(U, time_link_alpha3 / 2.0);

	// center

	int index = ggi(get_index(it, ix, iy, iz, T, L), mu);

	cm_eq_cm_ti_re(M2, gauge_field + index, 1.0 - time_link_alpha3);

	cm_pl_eq_cm(U, M2);

	// Projection to SU(3).
	cm_proj(U);
}

// idt = 0   -->   current timeslice
// idt = 1   -->   current timeslice+1

inline void HYP_Ref_3(int T, int L, double **U, int idt, int ix, int iy, int iz, int mu, int nu, int rho) {
	int index = get_index(idt, ix, iy, iz, T, L) * 18;

	if (mu == 0) {
		if ((nu == 2 && rho == 3) || (nu == 3 && rho == 2))
			*U = V3_0_1 + index;
		else if ((nu == 1 && rho == 3) || (nu == 3 && rho == 1))
			*U = V3_0_2 + index;
		else if ((nu == 1 && rho == 2) || (nu == 2 && rho == 1))
			*U = V3_0_3 + index;
		else {
			fprintf(stderr, "Error: inline void HYP_Ref_3(...\n");
			exit(EXIT_FAILURE);
		}
	} else if (mu == 1) {
		if ((nu == 0 && rho == 3) || (nu == 3 && rho == 0))
			*U = V3_1_2 + index;
		else if ((nu == 0 && rho == 2) || (nu == 2 && rho == 0))
			*U = V3_1_3 + index;
		else {
			fprintf(stderr, "Error: inline void HYP_Ref_3(...\n");
			exit(EXIT_FAILURE);
		}
	} else if (mu == 2) {
		if ((nu == 0 && rho == 3) || (nu == 3 && rho == 0))
			*U = V3_2_1 + index;
		else if ((nu == 0 && rho == 1) || (nu == 1 && rho == 0))
			*U = V3_2_3 + index;
		else {
			fprintf(stderr, "Error: inline void HYP_Ref_3(...\n");
			exit(EXIT_FAILURE);
		}
	} else if (mu == 3) {
		if ((nu == 0 && rho == 2) || (nu == 2 && rho == 0))
			*U = V3_3_1 + index;
		else if ((nu == 0 && rho == 1) || (nu == 1 && rho == 0))
			*U = V3_3_2 + index;
		else {
			fprintf(stderr, "Error: inline void HYP_Ref_3(...\n");
			exit(EXIT_FAILURE);
		}
	} else {
		fprintf(stderr, "Error: inline void HYP_Ref_3(...\n");
		exit(EXIT_FAILURE);
	}
}

inline void HYP_Helper_2(double *gauge_field, int T, int L, double time_link_alpha2, double *U, int it, int idt, int ix, int iy,
		int iz, int mu, int nu) {
	int rho;
	double M1[18], M2[18], *SU3_1, *SU3_2, *SU3_3;

	cm_eq_zero(U);

	if (mu == 0 && idt == 1)
		return;

	for (rho = 0; rho < 4; rho++)
	// For all staples.
			{
		if (rho == mu || rho == nu)
			continue;

		// negative rho-direction

		if (mu == 0) {
			if (rho == 1) {
				HYP_Ref_3(T, L, &SU3_1, idt, ix - 1, iy, iz, rho, mu, nu);
				HYP_Ref_3(T, L, &SU3_2, idt, ix - 1, iy, iz, mu, rho, nu);
				HYP_Ref_3(T, L, &SU3_3, idt + 1, ix - 1, iy, iz, rho, mu, nu);
			} else if (rho == 2) {
				HYP_Ref_3(T, L, &SU3_1, idt, ix, iy - 1, iz, rho, mu, nu);
				HYP_Ref_3(T, L, &SU3_2, idt, ix, iy - 1, iz, mu, rho, nu);
				HYP_Ref_3(T, L, &SU3_3, idt + 1, ix, iy - 1, iz, rho, mu, nu);
			} else if (rho == 3) {
				HYP_Ref_3(T, L, &SU3_1, idt, ix, iy, iz - 1, rho, mu, nu);
				HYP_Ref_3(T, L, &SU3_2, idt, ix, iy, iz - 1, mu, rho, nu);
				HYP_Ref_3(T, L, &SU3_3, idt + 1, ix, iy, iz - 1, rho, mu, nu);
			} else {
				fprintf(stderr, "Error: void HYP_Helper_2(...\n");
				exit(EXIT_FAILURE);
			}
		} else if (mu == 1) {
			if (rho == 0) {
				fprintf(stderr, "Error: inline void HYP_Helper_2(\n");
				exit(EXIT_FAILURE);

				/*
				 HYP_Ref_3(T, L, &SU3_1, idt-1, ix  , iy  , iz  , rho, mu, nu);
				 HYP_Ref_3(T, L, &SU3_2, idt-1, ix  , iy  , iz  , mu, rho, nu);
				 HYP_Ref_3(T, L, &SU3_3, idt-1, ix+1, iy  , iz  , rho, mu, nu);
				 */
			} else if (rho == 2) {
				HYP_Ref_3(T, L, &SU3_1, idt, ix, iy - 1, iz, rho, mu, nu);
				HYP_Ref_3(T, L, &SU3_2, idt, ix, iy - 1, iz, mu, rho, nu);
				HYP_Ref_3(T, L, &SU3_3, idt, ix + 1, iy - 1, iz, rho, mu, nu);
			} else if (rho == 3) {
				HYP_Ref_3(T, L, &SU3_1, idt, ix, iy, iz - 1, rho, mu, nu);
				HYP_Ref_3(T, L, &SU3_2, idt, ix, iy, iz - 1, mu, rho, nu);
				HYP_Ref_3(T, L, &SU3_3, idt, ix + 1, iy, iz - 1, rho, mu, nu);
			} else {
				fprintf(stderr, "Error: void HYP_Helper_2(...\n");
				exit(EXIT_FAILURE);
			}
		} else if (mu == 2) {
			if (rho == 0) {
				fprintf(stderr, "Error: inline void HYP_Helper_2(\n");
				exit(EXIT_FAILURE);

				/*
				 HYP_Ref_3(T, L, &SU3_1, idt-1, ix  , iy  , iz  , rho, mu, nu);
				 HYP_Ref_3(T, L, &SU3_2, idt-1, ix  , iy  , iz  , mu, rho, nu);
				 HYP_Ref_3(T, L, &SU3_3, idt-1, ix  , iy+1, iz  , rho, mu, nu);
				 */
			} else if (rho == 1) {
				HYP_Ref_3(T, L, &SU3_1, idt, ix - 1, iy, iz, rho, mu, nu);
				HYP_Ref_3(T, L, &SU3_2, idt, ix - 1, iy, iz, mu, rho, nu);
				HYP_Ref_3(T, L, &SU3_3, idt, ix - 1, iy + 1, iz, rho, mu, nu);
			} else if (rho == 3) {
				HYP_Ref_3(T, L, &SU3_1, idt, ix, iy, iz - 1, rho, mu, nu);
				HYP_Ref_3(T, L, &SU3_2, idt, ix, iy, iz - 1, mu, rho, nu);
				HYP_Ref_3(T, L, &SU3_3, idt, ix, iy + 1, iz - 1, rho, mu, nu);
			} else {
				fprintf(stderr, "Error: void HYP_Helper_2(...\n");
				exit(EXIT_FAILURE);
			}
		} else if (mu == 3) {
			if (rho == 0) {
				fprintf(stderr, "Error: inline void HYP_Helper_2(\n");
				exit(EXIT_FAILURE);

				/*
				 HYP_Ref_3(T, L, &SU3_1, idt-1, ix  , iy  , iz  , rho, mu, nu);
				 HYP_Ref_3(T, L, &SU3_2, idt-1, ix  , iy  , iz  , mu, rho, nu);
				 HYP_Ref_3(T, L, &SU3_3, idt-1, ix  , iy  , iz+1, rho, mu, nu);
				 */
			} else if (rho == 1) {
				HYP_Ref_3(T, L, &SU3_1, idt, ix - 1, iy, iz, rho, mu, nu);
				HYP_Ref_3(T, L, &SU3_2, idt, ix - 1, iy, iz, mu, rho, nu);
				HYP_Ref_3(T, L, &SU3_3, idt, ix - 1, iy, iz + 1, rho, mu, nu);
			} else if (rho == 2) {
				HYP_Ref_3(T, L, &SU3_1, idt, ix, iy - 1, iz, rho, mu, nu);
				HYP_Ref_3(T, L, &SU3_2, idt, ix, iy - 1, iz, mu, rho, nu);
				HYP_Ref_3(T, L, &SU3_3, idt, ix, iy - 1, iz + 1, rho, mu, nu);
			} else {
				fprintf(stderr, "Error: void HYP_Helper_2(...\n");
				exit(EXIT_FAILURE);
			}
		} else {
			fprintf(stderr, "Error: void HYP_Helper_2(...\n");
			exit(EXIT_FAILURE);
		}

		cm_eq_cm_ti_cm(M1, SU3_2, SU3_3);
		cm_eq_cm_dag_ti_cm(M2, SU3_1, M1);

		cm_pl_eq_cm(U, M2);

		// positive rho-direction

		if (mu == 0) {
			if (rho == 1) {
				HYP_Ref_3(T, L, &SU3_1, idt, ix, iy, iz, rho, mu, nu);
				HYP_Ref_3(T, L, &SU3_2, idt, ix + 1, iy, iz, mu, rho, nu);
				HYP_Ref_3(T, L, &SU3_3, idt + 1, ix, iy, iz, rho, mu, nu);
			} else if (rho == 2) {
				HYP_Ref_3(T, L, &SU3_1, idt, ix, iy, iz, rho, mu, nu);
				HYP_Ref_3(T, L, &SU3_2, idt, ix, iy + 1, iz, mu, rho, nu);
				HYP_Ref_3(T, L, &SU3_3, idt + 1, ix, iy, iz, rho, mu, nu);
			} else if (rho == 3) {
				HYP_Ref_3(T, L, &SU3_1, idt, ix, iy, iz, rho, mu, nu);
				HYP_Ref_3(T, L, &SU3_2, idt, ix, iy, iz + 1, mu, rho, nu);
				HYP_Ref_3(T, L, &SU3_3, idt + 1, ix, iy, iz, rho, mu, nu);
			} else {
				fprintf(stderr, "Error: void HYP_Helper_2(...\n");
				exit(EXIT_FAILURE);
			}
		} else if (mu == 1) {
			if (rho == 0) {
				fprintf(stderr, "Error: inline void HYP_Helper_2(\n");
				exit(EXIT_FAILURE);

				/*
				 HYP_Ref_3(T, L, &SU3_1, idt  , ix  , iy  , iz  , rho, mu, nu);
				 HYP_Ref_3(T, L, &SU3_2, idt+1, ix  , iy  , iz  , mu, rho, nu);
				 HYP_Ref_3(T, L, &SU3_3, idt  , ix+1, iy  , iz  , rho, mu, nu);
				 */
			} else if (rho == 2) {
				HYP_Ref_3(T, L, &SU3_1, idt, ix, iy, iz, rho, mu, nu);
				HYP_Ref_3(T, L, &SU3_2, idt, ix, iy + 1, iz, mu, rho, nu);
				HYP_Ref_3(T, L, &SU3_3, idt, ix + 1, iy, iz, rho, mu, nu);
			} else if (rho == 3) {
				HYP_Ref_3(T, L, &SU3_1, idt, ix, iy, iz, rho, mu, nu);
				HYP_Ref_3(T, L, &SU3_2, idt, ix, iy, iz + 1, mu, rho, nu);
				HYP_Ref_3(T, L, &SU3_3, idt, ix + 1, iy, iz, rho, mu, nu);
			} else {
				fprintf(stderr, "Error: void HYP_Helper_2(...\n");
				exit(EXIT_FAILURE);
			}
		} else if (mu == 2) {
			if (rho == 0) {
				fprintf(stderr, "Error: inline void HYP_Helper_2(\n");
				exit(EXIT_FAILURE);

				/*
				 HYP_Ref_3(T, L, &SU3_1, idt  , ix  , iy  , iz  , rho, mu, nu);
				 HYP_Ref_3(T, L, &SU3_2, idt+1, ix  , iy  , iz  , mu, rho, nu);
				 HYP_Ref_3(T, L, &SU3_3, idt  , ix  , iy+1, iz  , rho, mu, nu);
				 */
			} else if (rho == 1) {
				HYP_Ref_3(T, L, &SU3_1, idt, ix, iy, iz, rho, mu, nu);
				HYP_Ref_3(T, L, &SU3_2, idt, ix + 1, iy, iz, mu, rho, nu);
				HYP_Ref_3(T, L, &SU3_3, idt, ix, iy + 1, iz, rho, mu, nu);
			} else if (rho == 3) {
				HYP_Ref_3(T, L, &SU3_1, idt, ix, iy, iz, rho, mu, nu);
				HYP_Ref_3(T, L, &SU3_2, idt, ix, iy, iz + 1, mu, rho, nu);
				HYP_Ref_3(T, L, &SU3_3, idt, ix, iy + 1, iz, rho, mu, nu);
			} else {
				fprintf(stderr, "Error: void HYP_Helper_2(...\n");
				exit(EXIT_FAILURE);
			}
		} else if (mu == 3) {
			if (rho == 0) {
				fprintf(stderr, "Error: inline void HYP_Helper_2(\n");
				exit(EXIT_FAILURE);

				/*
				 HYP_Ref_3(T, L, &SU3_1, idt  , ix  , iy  , iz  , rho, mu, nu);
				 HYP_Ref_3(T, L, &SU3_2, idt+1, ix  , iy  , iz  , mu, rho, nu);
				 HYP_Ref_3(T, L, &SU3_3, idt  , ix  , iy  , iz+1, rho, mu, nu);
				 */
			} else if (rho == 1) {
				HYP_Ref_3(T, L, &SU3_1, idt, ix, iy, iz, rho, mu, nu);
				HYP_Ref_3(T, L, &SU3_2, idt, ix + 1, iy, iz, mu, rho, nu);
				HYP_Ref_3(T, L, &SU3_3, idt, ix, iy, iz + 1, rho, mu, nu);
			} else if (rho == 2) {
				HYP_Ref_3(T, L, &SU3_1, idt, ix, iy, iz, rho, mu, nu);
				HYP_Ref_3(T, L, &SU3_2, idt, ix, iy + 1, iz, mu, rho, nu);
				HYP_Ref_3(T, L, &SU3_3, idt, ix, iy, iz + 1, rho, mu, nu);
			} else {
				fprintf(stderr, "Error: void HYP_Helper_2(...\n");
				exit(EXIT_FAILURE);
			}
		} else {
			fprintf(stderr, "Error: void HYP_Helper_2(...\n");
			exit(EXIT_FAILURE);
		}

		cm_eq_cm_ti_cm_dag(M1, SU3_2, SU3_3);
		cm_eq_cm_ti_cm(M2, SU3_1, M1);

		cm_pl_eq_cm(U, M2);
	}

	cm_ti_eq_re(U, time_link_alpha2 / 4.0);

	// center

	int index = ggi(get_index(it + idt, ix, iy, iz, T, L), mu);

	cm_eq_cm_ti_re(M2, gauge_field + index, 1.0 - time_link_alpha2);

	cm_pl_eq_cm(U, M2);

	// Projection to SU(3).
	cm_proj(U);
}

// idt = 0   -->   current timeslice
// idt = 1   -->   current timeslice+1

inline void HYP_Ref_2(int T, int L, double **U, int idt, int ix, int iy, int iz, int mu, int nu) {
	int index = get_index(idt, ix, iy, iz, T, L) * 18;

	if (mu == 0) {
		if (nu == 3)
			*U = V2_0_12 + index;
		else if (nu == 2)
			*U = V2_0_13 + index;
		else if (nu == 1)
			*U = V2_0_23 + index;
		else {
			fprintf(stderr, "Error: inline void HYP_Ref_2(...\n");
			exit(EXIT_FAILURE);
		}
	} else if (mu == 1) {
		if (nu == 0)
			*U = V2_1_23 + index;
		else {
			fprintf(stderr, "Error: inline void HYP_Ref_2(...\n");
			exit(EXIT_FAILURE);
		}
	} else if (mu == 2) {
		if (nu == 0)
			*U = V2_2_13 + index;
		else {
			fprintf(stderr, "Error: inline void HYP_Ref_2(...\n");
			exit(EXIT_FAILURE);
		}
	} else if (mu == 3) {
		if (nu == 0)
			*U = V2_3_12 + index;
		else {
			fprintf(stderr, "Error: inline void HYP_Ref_2(...\n");
			exit(EXIT_FAILURE);
		}
	} else {
		fprintf(stderr, "Error: inline void HYP_Ref_2(...\n");
		exit(EXIT_FAILURE);
	}
}

inline void HYP_Helper_1(double *gauge_field, int T, int L, double time_link_alpha1, int it, int ix, int iy, int iz) {
	int index = ggi(get_index(it, ix, iy, iz, T, L), 0);
	double M1[18], M2[18], *SU3_1, *SU3_2, *SU3_3;

	double U[18];

	cm_eq_zero(U);

	// negative x-direction

	HYP_Ref_2(T, L, &SU3_1, 0, ix - 1, iy, iz, 1, 0);
	HYP_Ref_2(T, L, &SU3_2, 0, ix - 1, iy, iz, 0, 1);
	HYP_Ref_2(T, L, &SU3_3, 0 + 1, ix - 1, iy, iz, 1, 0);

	cm_eq_cm_ti_cm(M1, SU3_2, SU3_3);
	cm_eq_cm_dag_ti_cm(M2, SU3_1, M1);

	cm_pl_eq_cm(U, M2);

	// positive x-direction

	HYP_Ref_2(T, L, &SU3_1, 0, ix, iy, iz, 1, 0);
	HYP_Ref_2(T, L, &SU3_2, 0, ix + 1, iy, iz, 0, 1);
	HYP_Ref_2(T, L, &SU3_3, 0 + 1, ix, iy, iz, 1, 0);

	cm_eq_cm_ti_cm_dag(M1, SU3_2, SU3_3);
	cm_eq_cm_ti_cm(M2, SU3_1, M1);

	cm_pl_eq_cm(U, M2);

	// negative y-direction

	HYP_Ref_2(T, L, &SU3_1, 0, ix, iy - 1, iz, 2, 0);
	HYP_Ref_2(T, L, &SU3_2, 0, ix, iy - 1, iz, 0, 2);
	HYP_Ref_2(T, L, &SU3_3, 0 + 1, ix, iy - 1, iz, 2, 0);

	cm_eq_cm_ti_cm(M1, SU3_2, SU3_3);
	cm_eq_cm_dag_ti_cm(M2, SU3_1, M1);

	cm_pl_eq_cm(U, M2);

	// positive y-direction

	HYP_Ref_2(T, L, &SU3_1, 0, ix, iy, iz, 2, 0);
	HYP_Ref_2(T, L, &SU3_2, 0, ix, iy + 1, iz, 0, 2);
	HYP_Ref_2(T, L, &SU3_3, 0 + 1, ix, iy, iz, 2, 0);

	cm_eq_cm_ti_cm_dag(M1, SU3_2, SU3_3);
	cm_eq_cm_ti_cm(M2, SU3_1, M1);

	cm_pl_eq_cm(U, M2);

	// negative z-direction

	HYP_Ref_2(T, L, &SU3_1, 0, ix, iy, iz - 1, 3, 0);
	HYP_Ref_2(T, L, &SU3_2, 0, ix, iy, iz - 1, 0, 3);
	HYP_Ref_2(T, L, &SU3_3, 0 + 1, ix, iy, iz - 1, 3, 0);

	cm_eq_cm_ti_cm(M1, SU3_2, SU3_3);
	cm_eq_cm_dag_ti_cm(M2, SU3_1, M1);

	cm_pl_eq_cm(U, M2);

	// positive z-direction

	HYP_Ref_2(T, L, &SU3_1, 0, ix, iy, iz, 3, 0);
	HYP_Ref_2(T, L, &SU3_2, 0, ix, iy, iz + 1, 0, 3);
	HYP_Ref_2(T, L, &SU3_3, 0 + 1, ix, iy, iz, 3, 0);

	cm_eq_cm_ti_cm_dag(M1, SU3_2, SU3_3);
	cm_eq_cm_ti_cm(M2, SU3_1, M1);

	cm_pl_eq_cm(U, M2);

	cm_ti_eq_re(U, time_link_alpha1 / 6.0);

	// center

	cm_eq_cm_ti_re(M2, gauge_field + index, 1.0 - time_link_alpha1);

	cm_pl_eq_cm(U, M2);

	// Projection to SU(3).
	cm_proj(U);

	cm_eq_cm(gauge_field + index, U);
}

double *HYP_Time_Links_alloc(int L) {
	// two timeslices are needed
	int volume = 2 * L * L * L;

	fprintf(stderr, "  double *HYP_Time_Links_alloc(...   -->   Trying to allocate %d M ...",
			volume * 18 * sizeof(double) / 1000000);

	double *p_c;

	if ((p_c = (double *) malloc(volume * 18 * sizeof(double))) == NULL) {
		fprintf(stderr, "Error: double *HYP_Time_Links_alloc(...\n");
		exit(EXIT_FAILURE);
	}

	fprintf(stderr, " o.k.\n");

	return p_c;
}

void HYP_Time_Links(double *gauge_field, int T, int L, double time_link_alpha1, double time_link_alpha2,
		double time_link_alpha3) {
	int it, idt, ix, iy, iz;

	V2_0_12 = HYP_Time_Links_alloc(L);
	V2_0_13 = HYP_Time_Links_alloc(L);
	V2_0_23 = HYP_Time_Links_alloc(L);
	V2_1_23 = HYP_Time_Links_alloc(L);
	V2_2_13 = HYP_Time_Links_alloc(L);
	V2_3_12 = HYP_Time_Links_alloc(L);

	V3_0_1 = HYP_Time_Links_alloc(L);
	V3_0_2 = HYP_Time_Links_alloc(L);
	V3_0_3 = HYP_Time_Links_alloc(L);

	V3_1_2 = HYP_Time_Links_alloc(L);
	V3_1_3 = HYP_Time_Links_alloc(L);

	V3_2_1 = HYP_Time_Links_alloc(L);
	V3_2_3 = HYP_Time_Links_alloc(L);

	V3_3_1 = HYP_Time_Links_alloc(L);
	V3_3_2 = HYP_Time_Links_alloc(L);

	for (it = 0; it < T; it++) {
		fprintf(stderr, "  it = %2d, it = %2d, it = %2d, it = %2d ...\n", it, 0, 0, 0);

		// *****
		// *****
		// *****

		for (idt = 0; idt < 2; idt++) {
			// Decorated links with a single staple.

			for (ix = 0; ix < L; ix++) {
				for (iy = 0; iy < L; iy++) {
					for (iz = 0; iz < L; iz++) {
						int index = get_index(idt, ix, iy, iz, T, L) * 18;

						HYP_Helper_3(gauge_field, T, L, time_link_alpha3, V3_0_1 + index, it + idt, ix, iy, iz, 0, 2, 3);
						HYP_Helper_3(gauge_field, T, L, time_link_alpha3, V3_0_2 + index, it + idt, ix, iy, iz, 0, 1, 3);
						HYP_Helper_3(gauge_field, T, L, time_link_alpha3, V3_0_3 + index, it + idt, ix, iy, iz, 0, 1, 2);

						HYP_Helper_3(gauge_field, T, L, time_link_alpha3, V3_1_2 + index, it + idt, ix, iy, iz, 1, 0, 3);
						HYP_Helper_3(gauge_field, T, L, time_link_alpha3, V3_1_3 + index, it + idt, ix, iy, iz, 1, 0, 2);

						HYP_Helper_3(gauge_field, T, L, time_link_alpha3, V3_2_1 + index, it + idt, ix, iy, iz, 2, 0, 3);
						HYP_Helper_3(gauge_field, T, L, time_link_alpha3, V3_2_3 + index, it + idt, ix, iy, iz, 2, 0, 1);

						HYP_Helper_3(gauge_field, T, L, time_link_alpha3, V3_3_1 + index, it + idt, ix, iy, iz, 3, 0, 2);
						HYP_Helper_3(gauge_field, T, L, time_link_alpha3, V3_3_2 + index, it + idt, ix, iy, iz, 3, 0, 1);
					}
				}
			}
		}

		for (idt = 0; idt < 2; idt++) {
			// Decorated links with two staples.

			for (ix = 0; ix < L; ix++) {
				for (iy = 0; iy < L; iy++) {
					for (iz = 0; iz < L; iz++) {
						int index = get_index(idt, ix, iy, iz, T, L) * 18;

						HYP_Helper_2(gauge_field, T, L, time_link_alpha2, V2_0_12 + index, it, idt, ix, iy, iz, 0, 3);
						HYP_Helper_2(gauge_field, T, L, time_link_alpha2, V2_0_13 + index, it, idt, ix, iy, iz, 0, 2);
						HYP_Helper_2(gauge_field, T, L, time_link_alpha2, V2_0_23 + index, it, idt, ix, iy, iz, 0, 1);
						HYP_Helper_2(gauge_field, T, L, time_link_alpha2, V2_1_23 + index, it, idt, ix, iy, iz, 1, 0);
						HYP_Helper_2(gauge_field, T, L, time_link_alpha2, V2_2_13 + index, it, idt, ix, iy, iz, 2, 0);
						HYP_Helper_2(gauge_field, T, L, time_link_alpha2, V2_3_12 + index, it, idt, ix, iy, iz, 3, 0);
					}
				}
			}
		}

		// HYP links.

		for (ix = 0; ix < L; ix++) {
			for (iy = 0; iy < L; iy++) {
				for (iz = 0; iz < L; iz++) {
					HYP_Helper_1(gauge_field, T, L, time_link_alpha1, it, ix, iy, iz);
				}
			}
		}

		// *****
		// *****
		// *****

	}

	free(V3_0_1);
	free(V3_0_2);
	free(V3_0_3);

	free(V3_1_2);
	free(V3_1_3);

	free(V3_2_1);
	free(V3_2_3);

	free(V3_3_1);
	free(V3_3_2);

	free(V2_0_12);
	free(V2_0_13);
	free(V2_0_23);
	free(V2_1_23);
	free(V2_2_13);
	free(V2_3_12);
}

// ********************

// Performs an APE smearing step.

void APE_Smearing_Step(double *smeared_gauge_field, int T, int L, double APE_smearing_alpha) {

	double *smeared_gauge_field_old;
	Gauge_Field_Alloc(&smeared_gauge_field_old, T, L);
	Gauge_Field_Copy(smeared_gauge_field_old, smeared_gauge_field, T, L);

#pragma omp parallel for collapse(4)
	for (int it = 0; it < T; it++) {
		for (int ix = 0; ix < L; ix++) {
			for (int iy = 0; iy < L; iy++) {
				for (int iz = 0; iz < L; iz++) {
					int index;

					int index_mx_1, index_mx_2, index_mx_3;
					int index_px_1, index_px_2, index_px_3;
					int index_my_1, index_my_2, index_my_3;
					int index_py_1, index_py_2, index_py_3;
					int index_mz_1, index_mz_2, index_mz_3;
					int index_pz_1, index_pz_2, index_pz_3;

					double *U;
					double M1[18], M2[18];

					// **********
					// Links in x-direction.
					// **********

					index = ggi(get_index(it, ix, iy, iz, T, L), 1);

					index_my_1 = ggi(get_index(it, ix, iy - 1, iz, T, L), 2);
					index_my_2 = ggi(get_index(it, ix, iy - 1, iz, T, L), 1);
					index_my_3 = ggi(get_index(it, ix + 1, iy - 1, iz, T, L), 2);

					index_py_1 = ggi(get_index(it, ix, iy, iz, T, L), 2);
					index_py_2 = ggi(get_index(it, ix, iy + 1, iz, T, L), 1);
					index_py_3 = ggi(get_index(it, ix + 1, iy, iz, T, L), 2);

					index_mz_1 = ggi(get_index(it, ix, iy, iz - 1, T, L), 3);
					index_mz_2 = ggi(get_index(it, ix, iy, iz - 1, T, L), 1);
					index_mz_3 = ggi(get_index(it, ix + 1, iy, iz - 1, T, L), 3);

					index_pz_1 = ggi(get_index(it, ix, iy, iz, T, L), 3);
					index_pz_2 = ggi(get_index(it, ix, iy, iz + 1, T, L), 1);
					index_pz_3 = ggi(get_index(it, ix + 1, iy, iz, T, L), 3);

					U = smeared_gauge_field + index;
					cm_eq_zero(U);

					// negative y-direction

					cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_my_2, smeared_gauge_field_old + index_my_3);

//   cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_my_1, M1);
//   cm_pl_eq_cm(U, M2);
					cm_pl_eq_cm_dag_ti_cm(U, smeared_gauge_field_old + index_my_1, M1);

					// positive y-direction

					cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_py_2, smeared_gauge_field_old + index_py_3);

//   cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_py_1, M1);
//   cm_pl_eq_cm(U, M2);
					cm_pl_eq_cm_ti_cm(U, smeared_gauge_field_old + index_py_1, M1);

					// negative z-direction

					cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_mz_2, smeared_gauge_field_old + index_mz_3);

//   cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_mz_1, M1);
//   cm_pl_eq_cm(U, M2);
					cm_pl_eq_cm_dag_ti_cm(U, smeared_gauge_field_old + index_mz_1, M1);

					// positive z-direction

					cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_pz_2, smeared_gauge_field_old + index_pz_3);

//   cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_pz_1, M1);
//   cm_pl_eq_cm(U, M2);
					cm_pl_eq_cm_ti_cm(U, smeared_gauge_field_old + index_pz_1, M1);

					cm_ti_eq_re(U, APE_smearing_alpha);

					// center

					cm_pl_eq_cm(U, smeared_gauge_field_old + index);

					// Projection to SU(3).
					cm_proj(U);

					// **********
					// Links in y-direction.
					// **********

					index = ggi(get_index(it, ix, iy, iz, T, L), 2);

					index_mx_1 = ggi(get_index(it, ix - 1, iy, iz, T, L), 1);
					index_mx_2 = ggi(get_index(it, ix - 1, iy, iz, T, L), 2);
					index_mx_3 = ggi(get_index(it, ix - 1, iy + 1, iz, T, L), 1);

					index_px_1 = ggi(get_index(it, ix, iy, iz, T, L), 1);
					index_px_2 = ggi(get_index(it, ix + 1, iy, iz, T, L), 2);
					index_px_3 = ggi(get_index(it, ix, iy + 1, iz, T, L), 1);

					index_mz_1 = ggi(get_index(it, ix, iy, iz - 1, T, L), 3);
					index_mz_2 = ggi(get_index(it, ix, iy, iz - 1, T, L), 2);
					index_mz_3 = ggi(get_index(it, ix, iy + 1, iz - 1, T, L), 3);

					index_pz_1 = ggi(get_index(it, ix, iy, iz, T, L), 3);
					index_pz_2 = ggi(get_index(it, ix, iy, iz + 1, T, L), 2);
					index_pz_3 = ggi(get_index(it, ix, iy + 1, iz, T, L), 3);

					U = smeared_gauge_field + index;
					cm_eq_zero(U);

					// negative x-direction

					cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_mx_2, smeared_gauge_field_old + index_mx_3);

					cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_mx_1, M1);

					cm_pl_eq_cm(U, M2);

					// positive x-direction

					cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_px_2, smeared_gauge_field_old + index_px_3);

					cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_px_1, M1);

					cm_pl_eq_cm(U, M2);

					// negative z-direction

					cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_mz_2, smeared_gauge_field_old + index_mz_3);

					cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_mz_1, M1);

					cm_pl_eq_cm(U, M2);

					// positive z-direction

					cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_pz_2, smeared_gauge_field_old + index_pz_3);

					cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_pz_1, M1);

					cm_pl_eq_cm(U, M2);

					cm_ti_eq_re(U, APE_smearing_alpha);

					// center

					cm_pl_eq_cm(U, smeared_gauge_field_old + index);

					// Projection to SU(3).
					cm_proj(U);

					// **********
					// Links in z-direction.
					// **********

					index = ggi(get_index(it, ix, iy, iz, T, L), 3);

					index_mx_1 = ggi(get_index(it, ix - 1, iy, iz, T, L), 1);
					index_mx_2 = ggi(get_index(it, ix - 1, iy, iz, T, L), 3);
					index_mx_3 = ggi(get_index(it, ix - 1, iy, iz + 1, T, L), 1);

					index_px_1 = ggi(get_index(it, ix, iy, iz, T, L), 1);
					index_px_2 = ggi(get_index(it, ix + 1, iy, iz, T, L), 3);
					index_px_3 = ggi(get_index(it, ix, iy, iz + 1, T, L), 1);

					index_my_1 = ggi(get_index(it, ix, iy - 1, iz, T, L), 2);
					index_my_2 = ggi(get_index(it, ix, iy - 1, iz, T, L), 3);
					index_my_3 = ggi(get_index(it, ix, iy - 1, iz + 1, T, L), 2);

					index_py_1 = ggi(get_index(it, ix, iy, iz, T, L), 2);
					index_py_2 = ggi(get_index(it, ix, iy + 1, iz, T, L), 3);
					index_py_3 = ggi(get_index(it, ix, iy, iz + 1, T, L), 2);

					U = smeared_gauge_field + index;
					cm_eq_zero(U);

					// negative x-direction

					cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_mx_2, smeared_gauge_field_old + index_mx_3);

					cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_mx_1, M1);

					cm_pl_eq_cm(U, M2);

					// positive x-direction

					cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_px_2, smeared_gauge_field_old + index_px_3);

					cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_px_1, M1);

					cm_pl_eq_cm(U, M2);

					// negative y-direction

					cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_my_2, smeared_gauge_field_old + index_my_3);

					cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_my_1, M1);

					cm_pl_eq_cm(U, M2);

					// positive y-direction

					cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_py_2, smeared_gauge_field_old + index_py_3);

					cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_py_1, M1);

					cm_pl_eq_cm(U, M2);

					cm_ti_eq_re(U, APE_smearing_alpha);

					// center

					cm_pl_eq_cm(U, smeared_gauge_field_old + index);

					// Projection to SU(3).
					cm_proj(U);
				}
			}
		}
	}

	Gauge_Field_Free(&smeared_gauge_field_old);
}

// ********************

// Performs an APE smearing step on a given timeslice.

void APE_Smearing_Step_Timeslice(double *smeared_gauge_field, int T, int L, double APE_smearing_alpha, int timeslice) {
	int ix, iy, iz;
	double M1[18], M2[18];

	double *smeared_gauge_field_old;
	Timeslice_Gauge_Field_Alloc(&smeared_gauge_field_old, L);
	Timeslice_Gauge_Field_Copy(smeared_gauge_field_old, smeared_gauge_field, T, L, timeslice);

	for (ix = 0; ix < L; ix++) {
		for (iy = 0; iy < L; iy++) {
			for (iz = 0; iz < L; iz++) {
				int index, index_;

				int index_mx_1, index_mx_2, index_mx_3;
				int index_px_1, index_px_2, index_px_3;
				int index_my_1, index_my_2, index_my_3;
				int index_py_1, index_py_2, index_py_3;
				int index_mz_1, index_mz_2, index_mz_3;
				int index_pz_1, index_pz_2, index_pz_3;

				double *U;

				// **********
				// Links in x-direction.
				// **********

				index = ggi(get_index(timeslice, ix, iy, iz, T, L), 1);

				index_ = ggi(get_index(0, ix, iy, iz, T, L), 1);

				index_my_1 = ggi(get_index(0, ix, iy - 1, iz, T, L), 2);
				index_my_2 = ggi(get_index(0, ix, iy - 1, iz, T, L), 1);
				index_my_3 = ggi(get_index(0, ix + 1, iy - 1, iz, T, L), 2);

				index_py_1 = ggi(get_index(0, ix, iy, iz, T, L), 2);
				index_py_2 = ggi(get_index(0, ix, iy + 1, iz, T, L), 1);
				index_py_3 = ggi(get_index(0, ix + 1, iy, iz, T, L), 2);

				index_mz_1 = ggi(get_index(0, ix, iy, iz - 1, T, L), 3);
				index_mz_2 = ggi(get_index(0, ix, iy, iz - 1, T, L), 1);
				index_mz_3 = ggi(get_index(0, ix + 1, iy, iz - 1, T, L), 3);

				index_pz_1 = ggi(get_index(0, ix, iy, iz, T, L), 3);
				index_pz_2 = ggi(get_index(0, ix, iy, iz + 1, T, L), 1);
				index_pz_3 = ggi(get_index(0, ix + 1, iy, iz, T, L), 3);

				U = smeared_gauge_field + index;
				cm_eq_zero(U);

				// negative y-direction

				cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_my_2, smeared_gauge_field_old + index_my_3);

				cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_my_1, M1);

				cm_pl_eq_cm(U, M2);

				// positive y-direction

				cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_py_2, smeared_gauge_field_old + index_py_3);

				cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_py_1, M1);

				cm_pl_eq_cm(U, M2);

				// negative z-direction

				cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_mz_2, smeared_gauge_field_old + index_mz_3);

				cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_mz_1, M1);

				cm_pl_eq_cm(U, M2);

				// positive z-direction

				cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_pz_2, smeared_gauge_field_old + index_pz_3);

				cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_pz_1, M1);

				cm_pl_eq_cm(U, M2);

				cm_ti_eq_re(U, APE_smearing_alpha);

				// center

				cm_pl_eq_cm(U, smeared_gauge_field_old + index_);

				// Projection to SU(3).
				cm_proj(U);

				// **********
				// Links in y-direction.
				// **********

				index = ggi(get_index(timeslice, ix, iy, iz, T, L), 2);

				index_ = ggi(get_index(0, ix, iy, iz, T, L), 2);

				index_mx_1 = ggi(get_index(0, ix - 1, iy, iz, T, L), 1);
				index_mx_2 = ggi(get_index(0, ix - 1, iy, iz, T, L), 2);
				index_mx_3 = ggi(get_index(0, ix - 1, iy + 1, iz, T, L), 1);

				index_px_1 = ggi(get_index(0, ix, iy, iz, T, L), 1);
				index_px_2 = ggi(get_index(0, ix + 1, iy, iz, T, L), 2);
				index_px_3 = ggi(get_index(0, ix, iy + 1, iz, T, L), 1);

				index_mz_1 = ggi(get_index(0, ix, iy, iz - 1, T, L), 3);
				index_mz_2 = ggi(get_index(0, ix, iy, iz - 1, T, L), 2);
				index_mz_3 = ggi(get_index(0, ix, iy + 1, iz - 1, T, L), 3);

				index_pz_1 = ggi(get_index(0, ix, iy, iz, T, L), 3);
				index_pz_2 = ggi(get_index(0, ix, iy, iz + 1, T, L), 2);
				index_pz_3 = ggi(get_index(0, ix, iy + 1, iz, T, L), 3);

				U = smeared_gauge_field + index;
				cm_eq_zero(U);

				// negative x-direction

				cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_mx_2, smeared_gauge_field_old + index_mx_3);

				cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_mx_1, M1);

				cm_pl_eq_cm(U, M2);

				// positive x-direction

				cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_px_2, smeared_gauge_field_old + index_px_3);

				cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_px_1, M1);

				cm_pl_eq_cm(U, M2);

				// negative z-direction

				cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_mz_2, smeared_gauge_field_old + index_mz_3);

				cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_mz_1, M1);

				cm_pl_eq_cm(U, M2);

				// positive z-direction

				cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_pz_2, smeared_gauge_field_old + index_pz_3);

				cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_pz_1, M1);

				cm_pl_eq_cm(U, M2);

				cm_ti_eq_re(U, APE_smearing_alpha);

				// center

				cm_pl_eq_cm(U, smeared_gauge_field_old + index_);

				// Projection to SU(3).
				cm_proj(U);

				// **********
				// Links in z-direction.
				// **********

				index = ggi(get_index(timeslice, ix, iy, iz, T, L), 3);

				index_ = ggi(get_index(0, ix, iy, iz, T, L), 3);

				index_mx_1 = ggi(get_index(0, ix - 1, iy, iz, T, L), 1);
				index_mx_2 = ggi(get_index(0, ix - 1, iy, iz, T, L), 3);
				index_mx_3 = ggi(get_index(0, ix - 1, iy, iz + 1, T, L), 1);

				index_px_1 = ggi(get_index(0, ix, iy, iz, T, L), 1);
				index_px_2 = ggi(get_index(0, ix + 1, iy, iz, T, L), 3);
				index_px_3 = ggi(get_index(0, ix, iy, iz + 1, T, L), 1);

				index_my_1 = ggi(get_index(0, ix, iy - 1, iz, T, L), 2);
				index_my_2 = ggi(get_index(0, ix, iy - 1, iz, T, L), 3);
				index_my_3 = ggi(get_index(0, ix, iy - 1, iz + 1, T, L), 2);

				index_py_1 = ggi(get_index(0, ix, iy, iz, T, L), 2);
				index_py_2 = ggi(get_index(0, ix, iy + 1, iz, T, L), 3);
				index_py_3 = ggi(get_index(0, ix, iy, iz + 1, T, L), 2);

				U = smeared_gauge_field + index;
				cm_eq_zero(U);

				// negative x-direction

				cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_mx_2, smeared_gauge_field_old + index_mx_3);

				cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_mx_1, M1);

				cm_pl_eq_cm(U, M2);

				// positive x-direction

				cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_px_2, smeared_gauge_field_old + index_px_3);

				cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_px_1, M1);

				cm_pl_eq_cm(U, M2);

				// negative y-direction

				cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_my_2, smeared_gauge_field_old + index_my_3);

				cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_my_1, M1);

				cm_pl_eq_cm(U, M2);

				// positive y-direction

				cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_py_2, smeared_gauge_field_old + index_py_3);

				cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_py_1, M1);

				cm_pl_eq_cm(U, M2);

				cm_ti_eq_re(U, APE_smearing_alpha);

				// center

				cm_pl_eq_cm(U, smeared_gauge_field_old + index_);

				// Projection to SU(3).
				cm_proj(U);
			}
		}
	}

	Gauge_Field_Free(&smeared_gauge_field_old);
}

// ********************

// Performs a number of Jacobi smearing steps on a given timeslice.

// psi = quark spinor
// N, kappa = Jacobi smearing parameters
// timeslice = the timeslice, on which the smearing is performed

void Jacobi_Smearing_Steps(double *smeared_gauge_field, double *psi, int T, int L, int N, double kappa, int timeslice,
		bool use_slice_0_for_psi)  // for __OPTIMIZE_MEM__
		{
	int i1;
	int ix, iy, iz;

	int timeslice_psi = timeslice;

	if (use_slice_0_for_psi == true)
		timeslice_psi = 0;

	double *psi_old;
	Timeslice_Spinor_Field_Alloc(&psi_old, L);

	for (i1 = 0; i1 < N; i1++) {
		// Copy the timeslice of interest to psi_old.
		Timeslice_Spinor_Field_Copy(psi_old, psi, T, L, timeslice_psi);

		for (ix = 0; ix < L; ix++) {
			for (iy = 0; iy < L; iy++) {
				for (iz = 0; iz < L; iz++) {
					// *****
					// *****
					// *****
					// *****
					// *****

// Get indices.

					int index_s = gsi(get_index_timeslice(ix, iy, iz, T, L));

					int index_s_mx = gsi(get_index_timeslice(ix - 1, iy, iz, T, L));
					int index_s_px = gsi(get_index_timeslice(ix + 1, iy, iz, T, L));
					int index_s_my = gsi(get_index_timeslice(ix, iy - 1, iz, T, L));
					int index_s_py = gsi(get_index_timeslice(ix, iy + 1, iz, T, L));
					int index_s_mz = gsi(get_index_timeslice(ix, iy, iz - 1, T, L));
					int index_s_pz = gsi(get_index_timeslice(ix, iy, iz + 1, T, L));

					int index_g_mx = ggi(get_index(timeslice, ix - 1, iy, iz, T, L), 1);
					int index_g_px = ggi(get_index(timeslice, ix, iy, iz, T, L), 1);
					int index_g_my = ggi(get_index(timeslice, ix, iy - 1, iz, T, L), 2);
					int index_g_py = ggi(get_index(timeslice, ix, iy, iz, T, L), 2);
					int index_g_mz = ggi(get_index(timeslice, ix, iy, iz - 1, T, L), 3);
					int index_g_pz = ggi(get_index(timeslice, ix, iy, iz, T, L), 3);

					double *s = psi + gsi(get_index(timeslice_psi, ix, iy, iz, T, L));

					fv_eq_zero(s);

// negative x-direction

					fv_pl_eq_cm_dag_ti_fv(s, smeared_gauge_field + index_g_mx, psi_old + index_s_mx);

// positive x-direction

					fv_pl_eq_cm_ti_fv(s, smeared_gauge_field + index_g_px, psi_old + index_s_px);

// negative y-direction

					fv_pl_eq_cm_dag_ti_fv(s, smeared_gauge_field + index_g_my, psi_old + index_s_my);

// positive y-direction

					fv_pl_eq_cm_ti_fv(s, smeared_gauge_field + index_g_py, psi_old + index_s_py);

// negative z-direction

					fv_pl_eq_cm_dag_ti_fv(s, smeared_gauge_field + index_g_mz, psi_old + index_s_mz);

// positive z-direction

					fv_pl_eq_cm_ti_fv(s, smeared_gauge_field + index_g_pz, psi_old + index_s_pz);

// Put everything together; normalization.

					fv_ti_eq_re(s, kappa);

					fv_pl_eq_fv(s, psi_old + index_s);

					double norm = 1.0 / (1.0 + 6.0 * kappa);

					fv_ti_eq_re(s, norm);

					// *****
					// *****
					// *****
					// *****
					// *****
				}
			}
		}
	}

	Spinor_Field_Free(&psi_old);
}

// ********************
