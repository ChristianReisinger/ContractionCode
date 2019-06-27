// ********************



// sl_spectrum.cc

// Author: Marc Wagner
// Date: September 2007 (v1)
// Date: January 2009 (v2)



// ********************



/*

rm sl_spectrum.e*
rm sl_spectrum.o*
rm qsub_input_*

*/

#define __Q_SUB__

#define __OPTIMIZE_MEM__



// ********************



#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "fields.hh"
#include "geometry.hh"
#include "io.hh"
#include "linear_algebra.hh"
#include "propagator_io.hh"
#include "Q_phi.hh"
#include "smearing_techniques.hh"
#include "Wilson_loops.hh"



// ********************



// s = 0.5 * t.

inline void pl_id(double *s, const double *t)
{
  double spinor1[24];

  fv_eq_fv(s, t);

  fv_ti_eq_re(s, 0.5);
}

// s = 0.5 * (+1 - gamma_0) * t.

inline void pl_id_mi_gamma0(double *s, const double *t)
{
  double spinor1[24];

  fv_eq_fv(s, t);

  fv_eq_gamma_ti_fv(spinor1, 0, t);

  fv_mi_eq_fv(s, spinor1);
  fv_ti_eq_re(s, 0.5);
}

// s = 0.5 * (+1 + gamma_0) * t.

inline void pl_id_pl_gamma0(double *s, const double *t)
{
  double spinor1[24];

  fv_eq_fv(s, t);

  fv_eq_gamma_ti_fv(spinor1, 0, t);

  fv_pl_eq_fv(s, spinor1);
  fv_ti_eq_re(s, 0.5);
}


// Six paths meson creation operator.

void ExtendedOperator_6(

  double *PSI__S__gamma_5,
  double *PSI__P_mi__1,
  double *PSI__S__gamma_5_gamma_r,
  double *PSI__P_mi__gamma_r,

  double *PSI__P_pl__gamma_x_x_mi_gamma_y_y,
  double *PSI__P_pl__gamma_y_y_mi_gamma_z_z,
  double *PSI__P_pl__gamma_z_z_mi_gamma_x_x,
  double *PSI__D_mi__gamma_5_gamma_x_x_mi_gamma_y_y,
  double *PSI__D_mi__gamma_5_gamma_y_y_mi_gamma_z_z,
  double *PSI__D_mi__gamma_5_gamma_z_z_mi_gamma_x_x,

  double *psi, int it, int ix, int iy, int iz, bool ast, int operator_width,

  int index___, int delta_t___,  // only if __PRECOMPUTE_WIDTH_6__ is defined

  bool use_slice_0_for_psi = false);  // for __OPTIMIZE_MEM__


// Eight paths meson creation operator.

void ExtendedOperator_8(

  double *PSI__S__gamma_5,
  double *PSI__P_mi__1,
  double *PSI__S__gamma_5_gamma_r,
  double *PSI__P_mi__gamma_r,

  double *PSI__P_pl__gamma_x_x_mi_gamma_y_y,
  double *PSI__P_pl__gamma_y_y_mi_gamma_z_z,
  double *PSI__P_pl__gamma_z_z_mi_gamma_x_x,
  double *PSI__D_mi__gamma_5_gamma_x_x_mi_gamma_y_y,
  double *PSI__D_mi__gamma_5_gamma_y_y_mi_gamma_z_z,
  double *PSI__D_mi__gamma_5_gamma_z_z_mi_gamma_x_x,

  double *PSI__D_pl__gamma_x_y_z_pl_gamma_y_z_x_pl_gamma_z_x_y,
  double *PSI__F_mi__gamma_5_gamma_x_y_z_pl_gamma_y_z_x_pl_gamma_z_x_y,

  double *psi, int it, int ix, int iy, int iz, bool ast, int operator_width,

  int index___, int delta_t___,  // only if __PRECOMPUTE_WIDTH_8__ is defined

  bool use_slice_0_for_psi = false);  // for __OPTIMIZE_MEM__



// ********************



// Lattice parameters.


// Temporal extension of the lattice.
// int T = 4;
int T = 96;
// int T = 48;

// Spatial extension of the lattice.
// int L = 4;
int L = 48;
// int L = 24;

// kappa.
// double kappa = 0.160856;
double kappa = 0.154073;

// mu.
// double mu = 0.006400;
double mu = 0.0020;



// *****



// Path.
const char path_gauges[] = "XXXXXXXXXX";
const char path_propagators[] = "XXXXXXXXXX";
const char path_results[] = "XXXXXXXXXX";
// const char path_gauges[] = "../conf_L4T4_random_tslice2/";
// const char path_propagators[] = "../conf_L4T4_random_tslice2/";
// const char path_results[] = "./tmp/";

/*
const char path_gauges[] =
  "/work/mcwagner/mu0.0064/";
// const char path_propagators_arch[] =
//   "/arch/hch02/hch025/b3.9/L24T48/k0.160856_mu0.0064/prop/";
const char path_propagators[] =
  "/work/mcwagner/mu0.0064/propagators/";
*/


// Perform a check with another gauge field configuration obtained by a random
// gauge transformation.

// #define __RANDOM_GAUGE_TRANSFORMATION__
// #define __RANDOM_GAUGE_TRANSFORMATION_2__


// Gauge field configuration id.
int gauge_field_id = 0;

// Number of stochastic sources.
// const int phi_num = 3;
const int phi_num = 1;

// Timeslice of the stochastic source.
int timeslice = 2;

// Maximum temporal separation of the meson correlation functions.
// const int cor_meson_T_max = 1;
const int cor_meson_T_max = 18;

// Index of the propagator.
int index_prop = 0;



// *****



// Time link smearing parameters.

// #define __FAT_TIME_LINKS__
double time_link_epsilon = 0.0;

#define __HYP_TIME_LINKS__

/*
// HYP 1.
double time_link_alpha1 = 0.75;
double time_link_alpha2 = 0.60;
double time_link_alpha3 = 0.30;
*/

// /*
// HYP 2.
double time_link_alpha1 = 1.0;
double time_link_alpha2 = 1.0;
double time_link_alpha3 = 0.5;
// */

#ifdef __FAT_TIME_LINKS__
#ifdef __HYP_TIME_LINKS__
--> Error!!!
#endif
#endif


// APE smearing parameters.

// int APE_smearing_N = 2;
int APE_smearing_N = 10;

double APE_smearing_alpha = 0.5;
// double APE_smearing_alpha = 0.5;


// Number of different Jacobi smearing levels / operator widths.
// const int num_states = 1;
const int num_states = 5;


// Jacobi smearing parameters.

// int Jacobi_smearing_N[num_states] = {3};
int Jacobi_smearing_N[num_states] = {0, 30, 60, 60, 90};

// double Jacobi_smearing_kappa = 0.75;
double Jacobi_smearing_kappa = 0.5;


// Width of the meson creation operators.

// const int operator_width_6[num_states] = {2};
const int operator_width_6[num_states] = {0, 3, 3, 6, 3};

// #define __PRECOMPUTE_WIDTH_6__

#ifdef __PRECOMPUTE_WIDTH_6__
// index of a state in the following "precomputed paths tables"
int meson_paths_6_index[num_states];
// number of precomputed "path types"
int meson_paths_6_num;
// widths of the precomputed "path types" (i.g. not all entries are filled)
int meson_paths_6_width[num_states];
// precomputed paths
double *meson_paths_6;

const int meson_paths_6_gi_helper = (2 * cor_meson_T_max + 1) * L*L*L * 3 * 18;
inline int meson_paths_6_gi(const int index, const int delta_t, const int x, const int y, const int z, const int dir, const int T, const int L)
{
  return index * meson_paths_6_gi_helper +
    (3 * get_index(cor_meson_T_max + delta_t, x, y, z, T, L) + dir) * 18;
}
#endif

// const int operator_width_8[num_states] = {1};
const int operator_width_8[num_states] = {0, 2, 2, 4, 2};

// #define __PRECOMPUTE_WIDTH_8__

#ifdef __PRECOMPUTE_WIDTH_8__
// index of a state in the following "precomputed paths tables"
int meson_paths_8_index[num_states];
// number of precomputed "path types"
int meson_paths_8_num;
// widths of the precomputed "path types" (i.g. not all entries are filled)
int meson_paths_8_width[num_states];
// precomputed paths
double *meson_paths_8;

const int meson_paths_8_gi_helper = (2 * cor_meson_T_max + 1) * L*L*L * 4 * 18;
inline int meson_paths_8_gi(const int index, const int delta_t, const int x, const int y, const int z, const int dir, const int T, const int L)
{
  return index * meson_paths_8_gi_helper +
    (4 * get_index(cor_meson_T_max + delta_t, x, y, z, T, L) + dir) * 18;
}
#endif



// ********************



// Volume of the lattice.
int volume = T*L*L*L;


// Gauge field.

double *gauge_field;



// ********************



unsigned int time_cur, time_old = 0, time_start;

void Update_Time(const char *message)
{
  int i1;

  if(time_old == 0)
    {
      time_old = (unsigned int)time(NULL);
      time_start = time_old;
      return;
    }

  time_cur = (unsigned int)time(NULL);
  fprintf(stderr, "********** TIME **********\n");
  fprintf(stderr, "********** TIME **********\n");
  fprintf(stderr, "********** TIME **********\n");
  int time_delta = (int)(time_cur - time_old);
  int time_hh = time_delta / 3600;
  i1 = (time_delta - time_hh * 3600);
  int time_mm = i1 / 60;
  int time_ss = i1 - time_mm * 60;
  fprintf(stderr, "*** %02d:%02d:%02d hh:mm:ss  *** (%s)\n", time_hh, time_mm, time_ss, message);
  time_delta = (int)(time_cur - time_start);
  time_hh = time_delta / 3600;
  i1 = (time_delta - time_hh * 3600);
  time_mm = i1 / 60;
  time_ss = i1 - time_mm * 60;
  fprintf(stderr, "*** %02d:%02d:%02d total     ***\n", time_hh, time_mm, time_ss);
  fprintf(stderr, "********** TIME **********\n");
  fprintf(stderr, "********** TIME **********\n");
  fprintf(stderr, "********** TIME **********\n");
  fprintf(stderr, "\n");
  time_old = time_cur;
}

int main(int argc, char **argv)
{
  complex c1;
  int i1, i2, i3;
  int it, ix, iy, iz;
  int is, ic;
  char string1[1000], string2[1000];
  double SU3_1[18];
  double spinor1[24], spinor2[24];


  // **********
  Update_Time("");
  // **********


#ifdef __Q_SUB__

  // Format: gauge_field_id timeslice

  // Read the qsub index.

  if(fscanf(stdin, "%d %d %d", &gauge_field_id, &timeslice, &index_prop) != 3)
    {
      fprintf(stderr, "Error: int main(...\n");
      exit(EXIT_FAILURE);
    }

  fprintf(stderr, "--> gauge_field_id = %4d.\n", gauge_field_id);
  fprintf(stderr, "--> timeslice = %2d.\n", timeslice);
  fprintf(stderr, "--> index_prop = %2d.\n", index_prop);

#endif


  // **********


  // Gamma matrices.

  fprintf(stderr, "gamma_0 = \n");
  gamma_fprintf(0, stderr);
  fprintf(stderr, "\n");

  fprintf(stderr, "gamma_1 = \n");
  gamma_fprintf(1, stderr);
  fprintf(stderr, "\n");

  fprintf(stderr, "gamma_2 = \n");
  gamma_fprintf(2, stderr);
  fprintf(stderr, "\n");

  fprintf(stderr, "gamma_3 = \n");
  gamma_fprintf(3, stderr);
  fprintf(stderr, "\n");

  fprintf(stderr, "id = \n");
  gamma_fprintf(4, stderr);
  fprintf(stderr, "\n");

  fprintf(stderr, "gamma_5 = \n");
  gamma_fprintf(5, stderr);
  fprintf(stderr, "\n");


  // **********


  /*
  // Copy the propagators (to avoid problems with the Juelich file system).

  for(i1 = 0; i1 < 4; i1++)
    {
      for(i2 = 0; i2 < 10; i2++)
	{
	  char filename_src[1000], filename_dst[1000];

	  sprintf(filename_src, "%sconf.%04d.%02d.%d.inverted",
		  path_propagators_arch, gauge_field_id, timeslice, i1);

	  sprintf(filename_dst, "%sconf.%04d.%02d.%d.inverted",
		  path_propagators, gauge_field_id, timeslice, i1);

	  struct stat buf;
	  buf.st_size = 0;
	  stat(filename_dst, &buf);

	  fprintf(stderr, "stat: %s ...\n", filename_dst);

	  if(buf.st_size == 63700992)
	    {
	      fprintf(stderr, "  --> %d bytes (o.k.).\n", (int)(buf.st_size));
	      break;
	    }
	  else
	    {
	      fprintf(stderr, "  --> %d bytes (not o.k., 63700992 required).\n",
		      (int)(buf.st_size));
	    }

	  if(i2 == 10)
	    {
	      fprintf(stderr, "!!!!! ERROR !!!!!\n");
	      exit(EXIT_FAILURE);
	    }

	  sprintf(string1, "cp %s %s", filename_src, filename_dst);
	  fprintf(stderr, "%s ...\n", string1);
	  system(string1);
	}
    }

  fprintf(stderr, "\n");
  */


  // **********


  // Read the gauge field configuration.

#ifndef __OPTIMIZE_MEM__
  Gauge_Field_Alloc(&gauge_field, T, L);
#else
  Gauge_Field_Alloc(&gauge_field, 2*cor_meson_T_max+1, L);
#endif

  fprintf(stderr, "\n");

#ifndef __RANDOM_GAUGE_TRANSFORMATION_2__
#ifndef __RANDOM_GAUGE_TRANSFORMATION__
  sprintf(string1, "%sconf.%04d", path_gauges, gauge_field_id);
#else
  sprintf(string1, "%sconf.%04d.rnd", path_gauges, gauge_field_id);
#endif
#else
  sprintf(string1, "%sconf.%04d.rnd.rnd", path_gauges, gauge_field_id);
#endif

  fprintf(stderr, "Reading gauge field configuration '%s' ...\n", string1);

#ifndef __OPTIMIZE_MEM__

  read_lime_gauge_field_doubleprec(gauge_field, string1, T, L, L, L);

#else

  size_t offset = 0;

  if(timeslice - cor_meson_T_max < 0)
    {
      read_lime_gauge_field_doubleprec_timeslices(gauge_field + offset, string1, T, L, L, L, T + (timeslice - cor_meson_T_max), T - 1);

      offset += (cor_meson_T_max - timeslice)*L*L*L * 4 * 18;
    }

  int slice_i, slice_f;

  if(timeslice - cor_meson_T_max < 0)
    slice_i = 0;
  else
    slice_i = timeslice - cor_meson_T_max;

  if(timeslice + cor_meson_T_max >= T)
    slice_f = T - 1;
  else
    slice_f = timeslice + cor_meson_T_max;

  read_lime_gauge_field_doubleprec_timeslices(gauge_field + offset, string1, T, L, L, L, slice_i, slice_f);

  offset += (slice_f - slice_i + 1)*L*L*L * 4 * 18;

  if(timeslice + cor_meson_T_max >= T)
    {
      read_lime_gauge_field_doubleprec_timeslices(gauge_field + offset, string1, T, L, L, L, 0, (timeslice + cor_meson_T_max) - T);

      offset += (timeslice + cor_meson_T_max - T + 1)*L*L*L * 4 * 18;
    }

  if(offset != (2*cor_meson_T_max+1)*L*L*L * 4 * 18)
    {
      fprintf(stderr, "Error: int main(...\n");
      exit(EXIT_FAILURE);
    }

#endif

  fprintf(stderr, "  ... o.k.\n\n");


  // **********
  Update_Time("essentially reading the gauge field configuration");
  // **********


  // **********
  // **********
  // **********
  // Compute the static-light pseudoscalar correlation function.
  // **********
  // **********
  // **********


  fprintf(stderr, "Computing the static-light pseudoscalar correlation function ...\n\n");

  FILE *cor_meson_fd;

  sprintf(string1, "%scor_meson.%4d.%2d.%d", path_results, gauge_field_id, timeslice, index_prop);

  for(i1 = 0; i1 < strlen(string1); i1++)
    {
      if(string1[i1] == ' ')
	string1[i1] = '0';
    }

  if((cor_meson_fd = fopen(string1, "w")) == NULL)
    {
      fprintf(stderr, "Error: int main(...\n");
      exit(EXIT_FAILURE);
    }

  fprintf(cor_meson_fd, "%d\n\n", cor_meson_T_max);


  // six paths operators:

  //  0  -->   S__gamma_5   (A1)
  //  1  -->   P_mi__1   (A1)
  //  2  -->   S__gamma_5_gamma_r   (A1)
  //  3  -->   P_mi__gamma_r   (A1)

  //  4  -->   P_pl__gamma_x_x_mi_gamma_y_y   (E)
  //  5  -->   P_pl__gamma_y_y_mi_gamma_z_z   (E)
  //  6  -->   P_pl__gamma_z_z_mi_gamma_x_x   (E)
  //  7  -->   D_mi__gamma_5_gamma_x_x_mi_gamma_y_y   (E)
  //  8  -->   D_mi__gamma_5_gamma_y_y_mi_gamma_z_z   (E)
  //  9  -->   D_mi__gamma_5_gamma_z_z_mi_gamma_x_x   (E)

  // eight paths operators:

  // 10  -->   S__gamma_5   (A1)
  // 11  -->   P_mi__1   (A1)
  // 12  -->   S__gamma_5_gamma_r   (A1)
  // 13  -->   P_mi__gamma_r   (A1)

  // 14  -->   P_pl__gamma_x_x_mi_gamma_y_y   (E)
  // 15  -->   P_pl__gamma_y_y_mi_gamma_z_z   (E)
  // 16  -->   P_pl__gamma_z_z_mi_gamma_x_x   (E)
  // 17  -->   D_mi__gamma_5_gamma_x_x_mi_gamma_y_y   (E)
  // 18  -->   D_mi__gamma_5_gamma_y_y_mi_gamma_z_z   (E)
  // 19  -->   D_mi__gamma_5_gamma_z_z_mi_gamma_x_x   (E)

  // 20  -->   D_pl__gamma_x_y_z_pl_gamma_y_z_x_pl_gamma_z_x_y   (A2)
  // 21  -->   F_mi__gamma_5_gamma_x_y_z_pl_gamma_y_z_x_pl_gamma_z_x_y   (A2)

  const int num_operators = 22;

  for(i1 = 0; i1 < num_states; i1++)
    // For all smearing levels ...
    {
      for(i2 = 0; i2 < num_operators; i2++)
	// For all "different operators" ...
	{
	  int operator_index = i1 * num_operators + i2;


	  char operator_name[100];

	  if(i2 == 0)
	    sprintf(operator_name, "6__S__gamma_5");
	  if(i2 == 1)
	    sprintf(operator_name, "6__P_mi__1");
	  if(i2 == 2)
	    sprintf(operator_name, "6__S__gamma_5_gamma_r");
	  if(i2 == 3)
	    sprintf(operator_name, "6__P_mi__gamma_r");

	  if(i2 == 4)
	    sprintf(operator_name, "6__P_pl__gamma_x_x_mi_gamma_y_y");
	  if(i2 == 5)
	    sprintf(operator_name, "6__P_pl__gamma_y_y_mi_gamma_z_z");
	  if(i2 == 6)
	    sprintf(operator_name, "6__P_pl__gamma_z_z_mi_gamma_x_x");
	  if(i2 == 7)
	    sprintf(operator_name, "6__D_mi__gamma_5_gamma_x_x_mi_gamma_y_y");
	  if(i2 == 8)
	    sprintf(operator_name, "6__D_mi__gamma_5_gamma_y_y_mi_gamma_z_z");
	  if(i2 == 9)
	    sprintf(operator_name, "6__D_mi__gamma_5_gamma_z_z_mi_gamma_x_x");

	  if(i2 == 10)
	    sprintf(operator_name, "8__S__gamma_5");
	  if(i2 == 11)
	    sprintf(operator_name, "8__P_mi__1");
	  if(i2 == 12)
	    sprintf(operator_name, "8__S__gamma_5_gamma_r");
	  if(i2 == 13)
	    sprintf(operator_name, "8__P_mi__gamma_r");

	  if(i2 == 14)
	    sprintf(operator_name, "8__P_pl__gamma_x_x_mi_gamma_y_y");
	  if(i2 == 15)
	    sprintf(operator_name, "8__P_pl__gamma_y_y_mi_gamma_z_z");
	  if(i2 == 16)
	    sprintf(operator_name, "8__P_pl__gamma_z_z_mi_gamma_x_x");
	  if(i2 == 17)
	    sprintf(operator_name, "8__D_mi__gamma_5_gamma_x_x_mi_gamma_y_y");
	  if(i2 == 18)
	    sprintf(operator_name, "8__D_mi__gamma_5_gamma_y_y_mi_gamma_z_z");
	  if(i2 == 19)
	    sprintf(operator_name, "8__D_mi__gamma_5_gamma_z_z_mi_gamma_x_x");

	  if(i2 == 20)
	    sprintf(operator_name,
		    "8__D_pl__gamma_x_y_z_pl_gamma_y_z_x_pl_gamma_z_x_y");
	  if(i2 == 21)
	    sprintf(operator_name,
		    "8__F_mi__gamma_5_gamma_x_y_z_pl_gamma_y_z_x_pl_gamma_z_x_y");


	  if(i2 < 10)
	    {
	      fprintf(cor_meson_fd,
		      "%3d   %.6lf %.6lf %.6lf   %3d %.6lf   %3d %.6lf   %d   %s\n",
		      operator_index,
		      time_link_alpha1, time_link_alpha2, time_link_alpha3,
		      APE_smearing_N, APE_smearing_alpha,
		      Jacobi_smearing_N[i1], Jacobi_smearing_kappa,
		      operator_width_6[i1],
		      operator_name);
	    }
	  else
	    {
	      fprintf(cor_meson_fd,
		      "%3d   %.6lf %.6lf %.6lf   %3d %.6lf   %3d %.6lf   %d   %s\n",
		      operator_index,
		      time_link_alpha1, time_link_alpha2, time_link_alpha3,
		      APE_smearing_N, APE_smearing_alpha,
		      Jacobi_smearing_N[i1], Jacobi_smearing_kappa,
		      operator_width_8[i1],
		      operator_name);
	    }
	}
    }


  // The Wilson lines in time direction.

  double *Wilson_lines_pos, *Wilson_lines_neg;

  if((Wilson_lines_pos = (double *)malloc(L*L*L * 18 * sizeof(double))) == NULL)
    {
      fprintf(stderr, "Error: int main(...\n");
      exit(EXIT_FAILURE);
    }

  if((Wilson_lines_neg = (double *)malloc(L*L*L * 18 * sizeof(double))) == NULL)
    {
      fprintf(stderr, "Error: int main(...\n");
      exit(EXIT_FAILURE);
    }


  // Read the "propagators phi".

  double *phi[phi_num];

#ifndef __OPTIMIZE_MEM__

  for(i1 = 0; i1 < phi_num; i1++)
    {
      Spinor_Field_Alloc(phi + i1, T, L);

      // /*
#ifndef __RANDOM_GAUGE_TRANSFORMATION_2__
#ifndef __RANDOM_GAUGE_TRANSFORMATION__
      if(i1 == 0)
	sprintf(string1, "%spoint_source.c0.inverted", path_propagators);
      else if(i1 == 1)
	sprintf(string1, "%spoint_source.c1.inverted", path_propagators);
      else
	sprintf(string1, "%spoint_source.c2.inverted", path_propagators);
#else
      if(i1 == 0)
	sprintf(string1, "%spoint_source.rnd.c0.inverted", path_propagators);
      else if(i1 == 1)
	sprintf(string1, "%spoint_source.rnd.c1.inverted", path_propagators);
      else
	sprintf(string1, "%spoint_source.rnd.c2.inverted", path_propagators);
#endif
#else
      if(i1 == 0)
	sprintf(string1, "%spoint_source.rnd.rnd.c0.inverted", path_propagators);
      else if(i1 == 1)
	sprintf(string1, "%spoint_source.rnd.rnd.c1.inverted", path_propagators);
      else
	sprintf(string1, "%spoint_source.rnd.rnd.c2.inverted", path_propagators);
#endif
      // */

      /*
      sprintf(string1, "%sconf.%4d.%2d.%d.inverted",
	      path_propagators, gauge_field_id, timeslice, i1);

      for(i2 = 0; i2 < strlen(string1); i2++)
	{
	  if(string1[i2] == ' ')
	    string1[i2] = '0';
	}
      */

      fprintf(stderr, "Reading timeslice propagator '%s' ...\n", string1);

      // read_cmi(phi[i1], L, T, string1, kappa);

      // sprintf(string2, "%s", string1);
      sprintf(string2, "%s.lime", string1);
      read_lime_spinor(phi[i1], string2, 0, -1, T, L, L, L);

      /*
      // convert to a different format
      sprintf(string2, "%s.lime", string1);
      write_lime_spinor(phi[i1], string2, 0, 32, T, L, L, L);
      */

      fprintf(stderr, "  ... o.k.\n\n");
    }

#else

  for(i1 = 0; i1 < phi_num; i1++)
    {
      Spinor_Field_Alloc(phi + i1, 2*cor_meson_T_max+1, L);

      /*
#ifndef __RANDOM_GAUGE_TRANSFORMATION_2__
#ifndef __RANDOM_GAUGE_TRANSFORMATION__
      if(i1 == 0)
	sprintf(string1, "%spoint_source.c0.inverted", path_propagators);
      else if(i1 == 1)
	sprintf(string1, "%spoint_source.c1.inverted", path_propagators);
      else
	sprintf(string1, "%spoint_source.c2.inverted", path_propagators);
#else
      if(i1 == 0)
	sprintf(string1, "%spoint_source.rnd.c0.inverted", path_propagators);
      else if(i1 == 1)
	sprintf(string1, "%spoint_source.rnd.c1.inverted", path_propagators);
      else
	sprintf(string1, "%spoint_source.rnd.c2.inverted", path_propagators);
#endif
#else
      if(i1 == 0)
	sprintf(string1, "%spoint_source.rnd.rnd.c0.inverted", path_propagators);
      else if(i1 == 1)
	sprintf(string1, "%spoint_source.rnd.rnd.c1.inverted", path_propagators);
      else
	sprintf(string1, "%spoint_source.rnd.rnd.c2.inverted", path_propagators);
#endif
      */

      // /*
      // sprintf(string1, "%sconf.%4d.%2d.%d.inverted",
      //   path_propagators, gauge_field_id, timeslice, i1);
      sprintf(string1, "%sconf.%4d.%2d.%d.inverted",
	      path_propagators, gauge_field_id, timeslice, index_prop);

      for(i2 = 0; i2 < strlen(string1); i2++)
	{
	  if(string1[i2] == ' ')
	    string1[i2] = '0';
	}
      // */

      fprintf(stderr, "Reading timeslice propagator '%s' ...\n", string1);

      // sprintf(string2, "%s.lime", string1);
      sprintf(string2, "%s", string1);

      size_t offset = 0;

      for(it = 0; it <= 2*cor_meson_T_max; it++)
	{
	  read_lime_spinor_timeslice(phi[i1] + offset, string2, 0, (timeslice - cor_meson_T_max + it + T) % T, T, L, L, L);
	  offset += L*L*L * 24;
	}

      fprintf(stderr, "  ... o.k.\n\n");
    }

#endif


#ifndef __OPTIMIZE_MEM__

  // Compute the stochastic sources xi.

  double *xi[phi_num];

  fprintf(stderr, "Computing stochastic sources ...\n");

  for(i1 = 0; i1 < phi_num; i1++)
    {
      Spinor_Field_Alloc(xi + i1, T, L);

      Q_phi(xi[i1], phi[i1], gauge_field, T, L, kappa, mu);
    }

  fprintf(stderr, "  ... o.k.\n\n");


  /*
  // write to disk (Lime format)

  for(i1 = 0; i1 < phi_num; i1++)
    {
#ifndef __RANDOM_GAUGE_TRANSFORMATION_2__
#ifndef __RANDOM_GAUGE_TRANSFORMATION__
      if(i1 == 0)
	sprintf(string1, "%spoint_source.c0", path_propagators);
      else if(i1 == 1)
	sprintf(string1, "%spoint_source.c1", path_propagators);
      else
	sprintf(string1, "%spoint_source.c2", path_propagators);
#else
      if(i1 == 0)
	sprintf(string1, "%spoint_source.rnd.c0", path_propagators);
      else if(i1 == 1)
	sprintf(string1, "%spoint_source.rnd.c1", path_propagators);
      else
	sprintf(string1, "%spoint_source.rnd.c2", path_propagators);
#endif
#else
      if(i1 == 0)
	sprintf(string1, "%spoint_source.rnd.rnd.c0", path_propagators);
      else if(i1 == 1)
	sprintf(string1, "%spoint_source.rnd.rnd.c1", path_propagators);
      else
	sprintf(string1, "%spoint_source.rnd.rnd.c2", path_propagators);
#endif

      sprintf(string2, "%s.lime", string1);
      write_lime_spinor(xi[i1], string2, 0, 32, T, L, L, L);
    }
  */


  // Check, whether the stochastic sources are located on a single timeslice
  // and determine that timeslice.

  fprintf(stderr, "Checking stochastic sources ...\n");

  for(i1 = 0; i1 < phi_num; i1++)
    {
      for(it = 0; it < T; it++)
	{
	  for(ix = 0; ix < L; ix++)
	    {
	      for(iy = 0; iy < L; iy++)
		{
		  for(iz = 0; iz < L; iz++)
		    {
		      double *spinor =
			xi[i1] + gsi(get_index(it, ix, iy, iz, T, L));

		      // double epsilon = 0.000015;
		      double epsilon = 0.00050;

		      for(i2 = 0; i2 < 24; i2++)
			{
			  if(fabs(spinor[i2]) >= epsilon)
			    {
			      if(timeslice != it && timeslice != -1)
				{
				  fprintf(stderr, "it = %2d, ix = %2d, iy = %2d, iz = %2d, index = %2d   -->   xi = %+.11lf.\n", it, ix, iy, iz, i2, spinor[i2]);
				  fprintf(stderr, "  ... the stochastic source (index %d) is not a timeslice source!\n", i1);
				  exit(EXIT_FAILURE);
				}

			      timeslice = it;
			    }
			}
		    }
		}
	    }
	}
    }

  fprintf(stderr,
	  "  ... the stochastic sources are located on timeslice %2d.\n\n",
	  timeslice);

#else

  // Compute the stochastic sources xi.

  double *xi[phi_num];

  fprintf(stderr, "Computing stochastic sources ...\n");

  for(i1 = 0; i1 < phi_num; i1++)
    {
      Spinor_Field_Alloc(xi + i1, 1, L);

      Q_phi_timeslice(xi[i1], phi[i1], gauge_field, 2 * cor_meson_T_max + 1, L, kappa, mu, cor_meson_T_max, T);
    }

  fprintf(stderr, "  ... o.k.\n\n");


  // Check, whether the stochastic sources are located on a single timeslice
  // and determine that timeslice.

  fprintf(stderr, "Checking stochastic sources ...\n");

  for(i1 = 0; i1 < phi_num; i1++)
    {
      for(ix = 0; ix < L; ix++)
	{
	  for(iy = 0; iy < L; iy++)
	    {
	      for(iz = 0; iz < L; iz++)
		{
		  double *spinor =
		    xi[i1] + gsi(get_index(0, ix, iy, iz, T, L));

		  // double epsilon = 0.000015;
		  double epsilon = 0.001;

		  for(i2 = 0; i2 < 24; i2++)
		    {
		      if(fabs(fabs(spinor[i2]) - 1.0 / sqrt(2.0)) >= epsilon &&
			 fabs(spinor[i2]) >= epsilon)
			{
			  fprintf(stderr, "ix = %2d, iy = %2d, iz = %2d, index = %2d   -->   xi = %+.11lf (%+.11lf).\n", ix, iy, iz, i2, spinor[i2], fabs(spinor[i2]) - 1.0 / sqrt(2.0));
		      fprintf(stderr, "  ... the stochastic source (index %d) is not a timeslice source!\n", i1);
		      exit(EXIT_FAILURE);
			}
		    }
		}
	    }
	}
    }






  /*
  // Read the stochastic sources xi.

  double *xi[phi_num];

  fprintf(stderr, "Reading stochastic sources ...\n");

  for(i1 = 0; i1 < phi_num; i1++)
    {
      Spinor_Field_Alloc(xi + i1, 1, L);

#ifndef __RANDOM_GAUGE_TRANSFORMATION_2__
#ifndef __RANDOM_GAUGE_TRANSFORMATION__
      if(i1 == 0)
	sprintf(string1, "%spoint_source.c0", path_propagators);
      else if(i1 == 1)
	sprintf(string1, "%spoint_source.c1", path_propagators);
      else
	sprintf(string1, "%spoint_source.c2", path_propagators);
#else
      if(i1 == 0)
	sprintf(string1, "%spoint_source.rnd.c0", path_propagators);
      else if(i1 == 1)
	sprintf(string1, "%spoint_source.rnd.c1", path_propagators);
      else
	sprintf(string1, "%spoint_source.rnd.c2", path_propagators);
#endif
#else
      if(i1 == 0)
	sprintf(string1, "%spoint_source.rnd.rnd.c0", path_propagators);
      else if(i1 == 1)
	sprintf(string1, "%spoint_source.rnd.rnd.c1", path_propagators);
      else
	sprintf(string1, "%spoint_source.rnd.rnd.c2", path_propagators);
#endif

      sprintf(string2, "%s.lime", string1);
      read_lime_spinor_timeslice(xi[i1], string2, 0, timeslice, T, L, L, L);
    }

  fprintf(stderr, "  ... o.k.\n\n");
  */

#endif


  // Store the unsmeared stochastic sources.

  double *xi_org[phi_num];

  for(i1 = 0; i1 < phi_num; i1++)
    {
      Timeslice_Spinor_Field_Alloc(xi_org + i1, L);

#ifndef __OPTIMIZE_MEM__
      Timeslice_Spinor_Field_Copy(xi_org[i1], xi[i1], T, L, timeslice);
#else
      Timeslice_Spinor_Field_Copy(xi_org[i1], xi[i1], T, L, 0);
#endif
    }

  fprintf(stderr, "\n");


  // **********
  Update_Time("essentially reading propagators and computing sources");
  // **********


  // Apply link smearing techniques.


#ifdef __FAT_TIME_LINKS__

  // Fat links in time direction.

  fprintf(stderr, "Fat links are not optimized yet ...\n");
  fprintf(stderr, "Error: int main(...\n");
  exit(EXIT_FAILURE);

  /*
  fprintf(stderr,
          "Computing fat links in time direction (epsilon = %.3lf) ...\n",
	  time_link_epsilon);

  Fat_Time_Links(gauge_field, smeared_gauge_field, T, L, time_link_epsilon);

  fprintf(stderr, "  ... o.k.\n\n");
  */

#endif


#ifdef __HYP_TIME_LINKS__

  // HYP links in time direction.

  fprintf(stderr, "Computing HYP links in time direction (alpha_1 = %.3lf, alpha_2 = %.3lf, alpha_3 = %.3lf) ...\n", time_link_alpha1, time_link_alpha2, time_link_alpha3);

#ifndef __OPTIMIZE_MEM__

  HYP_Time_Links(gauge_field, T, L,
		 time_link_alpha1, time_link_alpha2, time_link_alpha3);

#else

  HYP_Time_Links(gauge_field, 2*cor_meson_T_max+1, L,
		 time_link_alpha1, time_link_alpha2, time_link_alpha3);

#endif

  fprintf(stderr, "  ... o.k.\n\n");

#endif


  // APE smearing of spatial links.

  for(i1 = 0; i1 < APE_smearing_N; i1++)
    {
      fprintf(stderr, "APE smearing (step %3d [%3d], alpha = %.3lf) ...\n",
	      i1+1, APE_smearing_N, APE_smearing_alpha);

#ifndef __OPTIMIZE_MEM__
      APE_Smearing_Step(gauge_field, T, L, APE_smearing_alpha);
#else
      APE_Smearing_Step(gauge_field, 2*cor_meson_T_max+1, L, APE_smearing_alpha);
#endif
    }

  if(APE_smearing_N > 0)
    fprintf(stderr, "  ... o.k.\n\n");


  // **********
  Update_Time("essentially link smearing");
  // **********


#ifdef __PRECOMPUTE_WIDTH_6__

#ifdef __OPTIMIZE_MEM__
--> Error!!!!!
#endif

  // Initialize.

  meson_paths_6_num = 0;

  for(i1 = 0; i1 < num_states; i1++)
    {
      if(operator_width_6[i1] == 0)
	{
	  meson_paths_6_index[i1] = -1;
	  continue;
	}

      for(i2 = 0; i2 < i1; i2++)
	{
	  if(operator_width_6[i1] == operator_width_6[i2])
	    break;
	}

      if(i2 < i1)
	{
	  meson_paths_6_index[i1] = meson_paths_6_index[i2];
	  continue;
	}

      meson_paths_6_index[i1] = meson_paths_6_num;
      meson_paths_6_width[meson_paths_6_num] = operator_width_6[i1];
      meson_paths_6_num++;
    }

  /*
  fprintf(stderr, "meson_paths_6_index = ");
  for(i1 = 0; i1 < num_states; i1++)
    fprintf(stderr, "%d ", meson_paths_6_index[i1]);
  fprintf(stderr, "\n");
  fprintf(stderr, "meson_paths_6_num = %d\n", meson_paths_6_num);
  fprintf(stderr, "meson_paths_6_width = ");
  for(i1 = 0; i1 < meson_paths_6_num; i1++)
    fprintf(stderr, "%d ", meson_paths_6_width[i1]);
  fprintf(stderr, "\n");
  */

  // Allocate memory.

  i1 = meson_paths_6_num * (2 * cor_meson_T_max + 1) * L*L*L * 3 * 18;

  fprintf(stderr,
	  "precomputing meson paths   -->   Trying to allocate %d M ...",
	  i1 * sizeof(double) / 1000000);

  if((meson_paths_6 = (double *)malloc(i1 * sizeof(double))) == NULL)
    {
      fprintf(stderr, "\nError (precomputing meson paths): int main(...\n");
      exit(EXIT_FAILURE);
    }

  fprintf(stderr, " o.k.\n\n");

  // Compute meson paths.

  for(i1 = 0; i1 < meson_paths_6_num; i1++)
    {
      for(it = -cor_meson_T_max; it <= +cor_meson_T_max; it++)
	{
	  for(ix = 0; ix < L; ix++)
	    {
	      for(iy = 0; iy < L; iy++)
		{
		  for(iz = 0; iz < L; iz++)
		    {
		      Path_1(meson_paths_6 + meson_paths_6_gi(i1, it, ix, iy, iz, 0, T, L), gauge_field, T, L, timeslice+it, ix, iy, iz, 1, 0, 0, meson_paths_6_width[i1]);
		      Path_1(meson_paths_6 + meson_paths_6_gi(i1, it, ix, iy, iz, 1, T, L), gauge_field, T, L, timeslice+it, ix, iy, iz, 0, 1, 0, meson_paths_6_width[i1]);
		      Path_1(meson_paths_6 + meson_paths_6_gi(i1, it, ix, iy, iz, 2, T, L), gauge_field, T, L, timeslice+it, ix, iy, iz, 0, 0, 1, meson_paths_6_width[i1]);
		    }
		}
	    }
	}
    }

#endif

#ifdef __PRECOMPUTE_WIDTH_8__

#ifdef __OPTIMIZE_MEM__
--> Error!!!!!
#endif

  // Initialize.

  meson_paths_8_num = 0;

  for(i1 = 0; i1 < num_states; i1++)
    {
      if(operator_width_8[i1] == 0)
	{
	  meson_paths_8_index[i1] = -1;
	  continue;
	}

      for(i2 = 0; i2 < i1; i2++)
	{
	  if(operator_width_8[i1] == operator_width_8[i2])
	    break;
	}

      if(i2 < i1)
	{
	  meson_paths_8_index[i1] = meson_paths_8_index[i2];
	  continue;
	}

      meson_paths_8_index[i1] = meson_paths_8_num;
      meson_paths_8_width[meson_paths_8_num] = operator_width_8[i1];
      meson_paths_8_num++;
    }

  /*
  fprintf(stderr, "meson_paths_8_index = ");
  for(i1 = 0; i1 < num_states; i1++)
    fprintf(stderr, "%d ", meson_paths_8_index[i1]);
  fprintf(stderr, "\n");
  fprintf(stderr, "meson_paths_8_num = %d\n", meson_paths_8_num);
  fprintf(stderr, "meson_paths_8_width = ");
  for(i1 = 0; i1 < meson_paths_8_num; i1++)
    fprintf(stderr, "%d ", meson_paths_8_width[i1]);
  fprintf(stderr, "\n");
  */

  // Allocate memory.

  i1 = meson_paths_8_num * (2 * cor_meson_T_max + 1) * L*L*L * 4 * 18;

  fprintf(stderr,
	  "precomputing meson paths   -->   Trying to allocate %d M ...",
	  i1 * sizeof(double) / 1000000);

  if((meson_paths_8 = (double *)malloc(i1 * sizeof(double))) == NULL)
    {
      fprintf(stderr, "\nError (precomputing meson paths): int main(...\n");
      exit(EXIT_FAILURE);
    }

  fprintf(stderr, " o.k.\n\n");

  // Compute meson paths.

  for(i1 = 0; i1 < meson_paths_8_num; i1++)
    {
      for(it = -cor_meson_T_max; it <= +cor_meson_T_max; it++)
	{
	  for(ix = 0; ix < L; ix++)
	    {
	      for(iy = 0; iy < L; iy++)
		{
		  for(iz = 0; iz < L; iz++)
		    {
		      // dir == 0   -->   (+--)
		      Diagonal_Path_3(meson_paths_8 + meson_paths_8_gi(i1, it, ix, iy, iz, 0, T, L), gauge_field, T, L, timeslice+it, ix, iy, iz, +1, -1, -1, meson_paths_8_width[i1]);
		      // dir == 1   -->   (+-+)
		      Diagonal_Path_3(meson_paths_8 + meson_paths_8_gi(i1, it, ix, iy, iz, 1, T, L), gauge_field, T, L, timeslice+it, ix, iy, iz, +1, -1, +1, meson_paths_8_width[i1]);
		      // dir == 2   -->   (++-)
		      Diagonal_Path_3(meson_paths_8 + meson_paths_8_gi(i1, it, ix, iy, iz, 2, T, L), gauge_field, T, L, timeslice+it, ix, iy, iz, +1, +1, -1, meson_paths_8_width[i1]);
		      // dir == 3   -->   (+++)
		      Diagonal_Path_3(meson_paths_8 + meson_paths_8_gi(i1, it, ix, iy, iz, 3, T, L), gauge_field, T, L, timeslice+it, ix, iy, iz, +1, +1, +1, meson_paths_8_width[i1]);
		    }
		}
	    }
	}
    }

#endif


  // **********
  Update_Time("precomputing meson paths");
  // **********


  // Quark smearing loop (Jacobi smearing).

  int index_states_1, index_states_2;

  for(index_states_1 = 0; index_states_1 < num_states; index_states_1++)
    // Different smearing levels for phi ...
    {
      // **********
      Update_Time("Jacobi smearing, phi [begin]");
      // **********

      fprintf(stderr, "Jacobi smearing (phi): N = %d, kappa = %.3lf ...\n",
	      Jacobi_smearing_N[index_states_1], Jacobi_smearing_kappa);

      for(i1 = 0; i1 < phi_num; i1++)
	{
	  for(i2 = -cor_meson_T_max; i2 <= cor_meson_T_max; i2++)
	    {
	      fprintf(stderr, "  index %d, timeslice %d ...\n  ", i1, i2);

#ifndef __OPTIMIZE_MEM__

	      if(index_states_1 == 0)
		{
  Jacobi_Smearing_Steps(gauge_field, phi[i1], T, L,
    Jacobi_smearing_N[0],
    Jacobi_smearing_kappa, timeslice+i2);
		}
	      else
		{
  Jacobi_Smearing_Steps(gauge_field, phi[i1], T, L,
    Jacobi_smearing_N[index_states_1] - Jacobi_smearing_N[index_states_1 - 1],
    Jacobi_smearing_kappa, timeslice+i2);
		}

#else

	      if(index_states_1 == 0)
		{
		  /*
		  if(i2 == 1)
		    {
		      fv_fprintf(phi[i1] + gsi(get_index(cor_meson_T_max + 1, 1, 2, 3, 2 * cor_meson_T_max + 1, L)), stderr);
		    }
		  */

  Jacobi_Smearing_Steps(gauge_field, phi[i1], 2 * cor_meson_T_max + 1, L,
    Jacobi_smearing_N[0],
    Jacobi_smearing_kappa, cor_meson_T_max+i2);

                  /*
		  if(i2 == 1)
		    {
		      fv_fprintf(phi[i1] + gsi(get_index(cor_meson_T_max + 1, 1, 2, 3, 2 * cor_meson_T_max + 1, L)), stderr);
		      exit(0);
		    }
		  */
		}
	      else
		{
  Jacobi_Smearing_Steps(gauge_field, phi[i1], 2 * cor_meson_T_max + 1, L,
    Jacobi_smearing_N[index_states_1] - Jacobi_smearing_N[index_states_1 - 1],
    Jacobi_smearing_kappa, cor_meson_T_max+i2);
		}

#endif
	    }
	}

      fprintf(stderr, "  ... o.k.\n\n");

      // **********
      Update_Time("Jacobi smearing, phi [end]");
      // **********


      // Restore the unsmeared stochastic source.

#ifndef __OPTIMIZE_MEM__

      int index_ts = gsi(get_index(timeslice, 0, 0, 0, T, L));

      for(i1 = 0; i1 < phi_num; i1++)
	{
	  memcpy(xi[i1] + index_ts, xi_org[i1], L*L*L * 24 * sizeof(double));
	}

#else

      for(i1 = 0; i1 < phi_num; i1++)
	{
	  memcpy(xi[i1], xi_org[i1], L*L*L * 24 * sizeof(double));
	}

#endif


      for(index_states_2 = 0; index_states_2 < num_states; index_states_2++)
	// Different smearing levels for xi ...
	{
	  // **********
	  Update_Time("Jacobi smearing, xi [begin]");
	  // **********

	  fprintf(stderr, "Jacobi smearing (xi ): N = %d, kappa = %.3lf ...\n",
	       Jacobi_smearing_N[index_states_2], Jacobi_smearing_kappa);

	  for(i1 = 0; i1 < phi_num; i1++)
	    {
	      fprintf(stderr, "  index %d ...\n  ", i1);

#ifndef __OPTIMIZE_MEM__

	      if(index_states_2 == 0)
		{
  Jacobi_Smearing_Steps(gauge_field, xi[i1], T, L,
    Jacobi_smearing_N[0],
    Jacobi_smearing_kappa, timeslice);
		}
	      else
		{
  Jacobi_Smearing_Steps(gauge_field, xi[i1], T, L,
    Jacobi_smearing_N[index_states_2] - Jacobi_smearing_N[index_states_2 - 1],
    Jacobi_smearing_kappa, timeslice);
		}

#else

	      if(index_states_2 == 0)
		{
  Jacobi_Smearing_Steps(gauge_field, xi[i1], 2 * cor_meson_T_max + 1, L,
    Jacobi_smearing_N[0],
			Jacobi_smearing_kappa, cor_meson_T_max, true);
		}
	      else
		{
  Jacobi_Smearing_Steps(gauge_field, xi[i1], 2 * cor_meson_T_max + 1, L,
    Jacobi_smearing_N[index_states_2] - Jacobi_smearing_N[index_states_2 - 1],
			Jacobi_smearing_kappa, cor_meson_T_max, true);
		}

#endif
	    }

	  fprintf(stderr, "  ... o.k.\n\n");

	  // **********
	  Update_Time("Jacobi smearing, xi [end]");
	  // **********


  // ********************
  // ********************
  // ********************
  // ********************
  // ********************
  // contraction loop
  // ********************
  // ********************
  // ********************
  // ********************
  // ********************

  // **********
  Update_Time("contraction [begin]");
  // **********


  // Compute the meson correlation function for the given smearing level.


  // Initialize the Wilson lines in time direction.

  for(ix = 0; ix < L; ix++)
    {
      for(iy = 0; iy < L; iy++)
	{
	  for(iz = 0; iz < L; iz++)
	    {
	      int index = get_index_timeslice(ix, iy, iz, T, L) * 18;

	      cm_eq_id(Wilson_lines_pos + index);
	      cm_eq_id(Wilson_lines_neg + index);
	    }
	}
    }


  complex
    cor_meson_pos[num_operators][num_operators][cor_meson_T_max+1],
    cor_meson_neg[num_operators][num_operators][cor_meson_T_max+1];


  int cor_meson_T;

  for(cor_meson_T = 0; cor_meson_T <= cor_meson_T_max; cor_meson_T++)
    {
      if(cor_meson_T > 0)
	{
	  // Compute Wilson lines in time direction.

	  for(ix = 0; ix < L; ix++)
	    {
	      for(iy = 0; iy < L; iy++)
		{
		  for(iz = 0; iz < L; iz++)
		    {
  int index_WL = get_index_timeslice(ix, iy, iz, T, L) * 18;
  int index_gf;

#ifndef __OPTIMIZE_MEM__
  index_gf = ggi(get_index(timeslice+cor_meson_T-1, ix, iy, iz, T, L), 0);
#else
  index_gf = ggi(get_index(cor_meson_T_max+cor_meson_T-1, ix, iy, iz, 2 * cor_meson_T_max + 1, L), 0);
#endif

  cm_eq_cm_ti_cm(SU3_1, Wilson_lines_pos + index_WL,
		 gauge_field + index_gf);

  cm_eq_cm(Wilson_lines_pos + index_WL, SU3_1);

#ifndef __OPTIMIZE_MEM__
  index_gf = ggi(get_index(timeslice-cor_meson_T, ix, iy, iz, T, L), 0);
#else
  index_gf = ggi(get_index(cor_meson_T_max-cor_meson_T, ix, iy, iz, 2 * cor_meson_T_max + 1, L), 0);
#endif

  cm_eq_cm_ti_cm_dag(SU3_1, Wilson_lines_neg + index_WL,
		     gauge_field + index_gf);

  cm_eq_cm(Wilson_lines_neg + index_WL, SU3_1);
		    }
		}
	    }
	}


      for(i1 = 0; i1 < num_operators; i1++)
	{
	  for(i2 = 0; i2 < num_operators; i2++)
	    {
	      cor_meson_pos[i1][i2][cor_meson_T].re = 0.0;
	      cor_meson_pos[i1][i2][cor_meson_T].im = 0.0;

	      cor_meson_neg[i1][i2][cor_meson_T].re = 0.0;
	      cor_meson_neg[i1][i2][cor_meson_T].im = 0.0;
	    }
	}


      for(ix = 0; ix < L; ix++)
	{
	  fprintf(stderr, "  ... x = %2d, y = %2d, z = %2d ...\n", ix, 0, 0);

	  for(iy = 0; iy < L; iy++)
	    {
	      for(iz = 0; iz < L; iz++)
		{
                  // *****
                  // *****
                  // *****
                  // *****
                  // *****

int index_WL = get_index_timeslice(ix, iy, iz, T, L) * 18;

double *WL_pos_ = Wilson_lines_pos + index_WL;
double *WL_neg_ = Wilson_lines_neg + index_WL;

for(i1 = 0; i1 < phi_num; i1++)
  {
    double XI[num_operators][24];

    int index___ = 0;
    int delta_t___ = 0;

#ifdef __PRECOMPUTE_WIDTH_6__
    index___ = meson_paths_6_index[index_states_2];
    delta_t___ = 0;
#endif

#ifndef __OPTIMIZE_MEM__

    ExtendedOperator_6(XI[ 0], XI[ 1], XI[ 2], XI[ 3], XI[ 4],
		       XI[ 5], XI[ 6], XI[ 7], XI[ 8], XI[ 9],
		       xi[i1], timeslice, ix, iy, iz, true,
		       operator_width_6[index_states_2], index___, delta_t___);

#else

    ExtendedOperator_6(XI[ 0], XI[ 1], XI[ 2], XI[ 3], XI[ 4],
		       XI[ 5], XI[ 6], XI[ 7], XI[ 8], XI[ 9],
		       xi[i1], cor_meson_T_max, ix, iy, iz, true,
		       operator_width_6[index_states_2], index___, delta_t___,
		       true);

#endif

#ifdef __PRECOMPUTE_WIDTH_8__
    index___ = meson_paths_8_index[index_states_2];
    delta_t___ = 0;
#endif

#ifndef __OPTIMIZE_MEM__

    ExtendedOperator_8(XI[10], XI[11], XI[12], XI[13], XI[14],
		       XI[15], XI[16], XI[17], XI[18], XI[19],
		       XI[20], XI[21],
		       xi[i1], timeslice, ix, iy, iz, true,
		       operator_width_8[index_states_2], index___, delta_t___);

#else

    ExtendedOperator_8(XI[10], XI[11], XI[12], XI[13], XI[14],
		       XI[15], XI[16], XI[17], XI[18], XI[19],
		       XI[20], XI[21],
		       xi[i1], cor_meson_T_max, ix, iy, iz, true,
		       operator_width_8[index_states_2], index___, delta_t___,
		       true);

#endif


    // positive t direction (\bar{\psi} \phi correlation)

    double PHI_pos[num_operators][24];

#ifdef __PRECOMPUTE_WIDTH_6__
    index___ = meson_paths_6_index[index_states_1];
    delta_t___ = +cor_meson_T;
#endif

#ifndef __OPTIMIZE_MEM__

    ExtendedOperator_6(PHI_pos[ 0], PHI_pos[ 1], PHI_pos[ 2], PHI_pos[ 3],
		       PHI_pos[ 4], PHI_pos[ 5], PHI_pos[ 6], PHI_pos[ 7],
		       PHI_pos[ 8], PHI_pos[ 9],
		       phi[i1], timeslice+cor_meson_T, ix, iy, iz, false,
		       operator_width_6[index_states_1], index___, delta_t___);

#else

    ExtendedOperator_6(PHI_pos[ 0], PHI_pos[ 1], PHI_pos[ 2], PHI_pos[ 3],
		       PHI_pos[ 4], PHI_pos[ 5], PHI_pos[ 6], PHI_pos[ 7],
		       PHI_pos[ 8], PHI_pos[ 9],
		       phi[i1], cor_meson_T_max+cor_meson_T, ix, iy, iz, false,
		       operator_width_6[index_states_1], index___, delta_t___);
    /*
    if(cor_meson_T == 1)
      {
	fv_fprintf(PHI_pos[0], stderr);
	exit(0);
      }
    */

#endif

#ifdef __PRECOMPUTE_WIDTH_8__
    index___ = meson_paths_8_index[index_states_1];
    delta_t___ = +cor_meson_T;
#endif

#ifndef __OPTIMIZE_MEM__

    ExtendedOperator_8(PHI_pos[10], PHI_pos[11], PHI_pos[12], PHI_pos[13],
		       PHI_pos[14], PHI_pos[15], PHI_pos[16], PHI_pos[17],
		       PHI_pos[18], PHI_pos[19], PHI_pos[20], PHI_pos[21],
		       phi[i1], timeslice+cor_meson_T, ix, iy, iz, false,
		       operator_width_8[index_states_1], index___, delta_t___);

#else

    ExtendedOperator_8(PHI_pos[10], PHI_pos[11], PHI_pos[12], PHI_pos[13],
		       PHI_pos[14], PHI_pos[15], PHI_pos[16], PHI_pos[17],
		       PHI_pos[18], PHI_pos[19], PHI_pos[20], PHI_pos[21],
		       phi[i1], cor_meson_T_max+cor_meson_T, ix, iy, iz, false,
		       operator_width_8[index_states_1], index___, delta_t___);

#endif

    for(i2 = 0; i2 < num_operators; i2++)
      {
	fv_eq_cm_ti_fv(spinor1, WL_pos_, PHI_pos[i2]);

	if(cor_meson_T == 0)
	  pl_id(spinor2, spinor1);
	else
	  pl_id_mi_gamma0(spinor2, spinor1);

	for(i3 = 0; i3 < num_operators; i3++)
	  {
	    co_eq_fv_dag_ti_fv(&c1, XI[i3], spinor2);
	    // fprintf(stdout, "%2d   %2d %2d %2d   %d   %3d   %+9.6lf %+9.6lf I\n", cor_meson_T, ix, iy, iz, i1, i2, c1.re, c1.im);

	    cor_meson_pos[i2][i3][cor_meson_T].re -= c1.re;
	    cor_meson_pos[i2][i3][cor_meson_T].im -= c1.im;
	  }
      }


    // negative t direction (\bar{\phi} \psi correlation)

    double PHI_neg[num_operators][24];

#ifdef __PRECOMPUTE_WIDTH_6__
    index___ = meson_paths_6_index[index_states_1];
    delta_t___ = -cor_meson_T;
#endif

#ifndef __OPTIMIZE_MEM__

    ExtendedOperator_6(PHI_neg[ 0], PHI_neg[ 1], PHI_neg[ 2], PHI_neg[ 3],
		       PHI_neg[ 4], PHI_neg[ 5], PHI_neg[ 6], PHI_neg[ 7],
		       PHI_neg[ 8], PHI_neg[ 9],
		       phi[i1], timeslice-cor_meson_T, ix, iy, iz, false, 
		       operator_width_6[index_states_1], index___, delta_t___);

#else

    ExtendedOperator_6(PHI_neg[ 0], PHI_neg[ 1], PHI_neg[ 2], PHI_neg[ 3],
		       PHI_neg[ 4], PHI_neg[ 5], PHI_neg[ 6], PHI_neg[ 7],
		       PHI_neg[ 8], PHI_neg[ 9],
		       phi[i1], cor_meson_T_max-cor_meson_T, ix, iy, iz, false, 
		       operator_width_6[index_states_1], index___, delta_t___);

#endif

#ifdef __PRECOMPUTE_WIDTH_8__
    index___ = meson_paths_8_index[index_states_1];
    delta_t___ = -cor_meson_T;
#endif

#ifndef __OPTIMIZE_MEM__

    ExtendedOperator_8(PHI_neg[10], PHI_neg[11], PHI_neg[12], PHI_neg[13],
		       PHI_neg[14], PHI_neg[15], PHI_neg[16], PHI_neg[17],
		       PHI_neg[18], PHI_neg[19], PHI_neg[20], PHI_neg[21],
		       phi[i1], timeslice-cor_meson_T, ix, iy, iz, false,
		       operator_width_8[index_states_1], index___, delta_t___);

#else

    ExtendedOperator_8(PHI_neg[10], PHI_neg[11], PHI_neg[12], PHI_neg[13],
		       PHI_neg[14], PHI_neg[15], PHI_neg[16], PHI_neg[17],
		       PHI_neg[18], PHI_neg[19], PHI_neg[20], PHI_neg[21],
		       phi[i1], cor_meson_T_max-cor_meson_T, ix, iy, iz, false,
		       operator_width_8[index_states_1], index___, delta_t___);

#endif

    for(i2 = 0; i2 < num_operators; i2++)
      {
	// (-1)^(# gamma_5) for operator 1
	int sign_op1 = +1;
	if(i2 ==  0 || i2 ==  2 ||
	   i2 ==  7 || i2 ==  8 || i2 ==  9 ||
	   i2 == 10 || i2 == 12 ||
	   i2 == 17 || i2 == 18 || i2 == 19 ||
	   i2 == 21)
	  sign_op1 = -1;

	fv_eq_cm_ti_fv(spinor1, WL_neg_, PHI_neg[i2]);

	if(cor_meson_T == 0)
	  pl_id(spinor2, spinor1);
	else
	  pl_id_pl_gamma0(spinor2, spinor1);

	for(i3 = 0; i3 < num_operators; i3++)
	  {
	    // (-1)^(# gamma_5) for operator 2
	    int sign_op2 = +1;
	    if(i3 ==  0 || i3 ==  2 ||
	       i3 ==  7 || i3 ==  8 || i3 ==  9 ||
	       i3 == 10 || i3 == 12 ||
	       i3 == 17 || i3 == 18 || i3 == 19 ||
	       i3 == 21)
	      sign_op2 = -1;

	    co_eq_fv_dag_ti_fv(&c1, XI[i3], spinor2);

	    if(sign_op1 == sign_op2)
	      { 
		cor_meson_neg[i2][i3][cor_meson_T].re -= c1.re;
		cor_meson_neg[i2][i3][cor_meson_T].im += c1.im;
	      }
	    else
	      { 
		cor_meson_neg[i2][i3][cor_meson_T].re += c1.re;
		cor_meson_neg[i2][i3][cor_meson_T].im -= c1.im;
	      }
	  }
      }
  }

                  // *****
                  // *****
                  // *****
                  // *****
                  // *****
		}
	    }
	}


      // Divide by the spatial volume.

      double spatial_volume = (double)(L*L*L);

      for(i2 = 0; i2 < num_operators; i2++)
	{
	  for(i3 = 0; i3 < num_operators; i3++)
	    {
	      cor_meson_pos[i2][i3][cor_meson_T].re /= spatial_volume;
	      cor_meson_pos[i2][i3][cor_meson_T].im /= spatial_volume;

	      cor_meson_neg[i2][i3][cor_meson_T].re /= spatial_volume;
	      cor_meson_neg[i2][i3][cor_meson_T].im /= spatial_volume;
	    }
	}


      // Include antiperiodic boundary conditions.

      complex B;
      B.re = cos(M_PI * cor_meson_T / T);
      B.im = sin(M_PI * cor_meson_T / T);

      for(i2 = 0; i2 < num_operators; i2++)
	{
	  for(i3 = 0; i3 < num_operators; i3++)
	    co_ti_eq_co(cor_meson_pos[i2][i3] + cor_meson_T, &B);
	}

      // B.im = -B.im;

      for(i2 = 0; i2 < num_operators; i2++)
	{
	  for(i3 = 0; i3 < num_operators; i3++)
	    co_ti_eq_co(cor_meson_neg[i2][i3] + cor_meson_T, &B);
	}


      // screen output

      fprintf(stderr,
	      "  cor_meson_pos[0][0][%2d] = %+9.6lf %+9.6lf I (%+.11e %+.11e I).\n",
	      cor_meson_T,
	      cor_meson_pos[0][0][cor_meson_T].re,
	      cor_meson_pos[0][0][cor_meson_T].im,
	      cor_meson_pos[0][0][cor_meson_T].re,
	      cor_meson_pos[0][0][cor_meson_T].im);
      /*
      fprintf(stderr,
	      "cor_meson_pos[2][2][%2d] = %+9.6lf %+9.6lf I (%+.11e %+.11e I).\n",
	      cor_meson_T,
	      cor_meson_pos[2][2][cor_meson_T].re,
	      cor_meson_pos[2][2][cor_meson_T].im,
	      cor_meson_pos[2][2][cor_meson_T].re,
	      cor_meson_pos[2][2][cor_meson_T].im);
      */

      fprintf(stderr,
	      "  cor_meson_neg[0][0][%2d] = %+9.6lf %+9.6lf I (%+.11e %+.11e I).\n",
	      cor_meson_T,
	      cor_meson_neg[0][0][cor_meson_T].re,
	      cor_meson_neg[0][0][cor_meson_T].im,
	      cor_meson_neg[0][0][cor_meson_T].re,
	      cor_meson_neg[0][0][cor_meson_T].im);
      /*
      fprintf(stderr,
	      "cor_meson_neg[2][2][%2d] = %+9.6lf %+9.6lf I (%+.11e %+.11e I).\n",
	      cor_meson_T,
	      cor_meson_neg[2][2][cor_meson_T].re,
	      cor_meson_neg[2][2][cor_meson_T].im,
	      cor_meson_neg[2][2][cor_meson_T].re,
	      cor_meson_neg[2][2][cor_meson_T].im);
      */


      // file output

      if(cor_meson_T == 0)
	fprintf(cor_meson_fd, "\n");

      for(i2 = 0; i2 < num_operators; i2++)
	{
	  for(i3 = 0; i3 < num_operators; i3++)
	    {
	      // positive time direction

	      // index_states_1   -->   bra (phi)
	      // index_states_2   -->   ket (xi)

	      // i2   -->   phi (bra)
	      // i3   -->   xi (ket)

	      int index_bra = index_states_1 * num_operators + i2;
	      int index_ket = index_states_2 * num_operators + i3;

	      // fprintf(cor_meson_fd, "%2d pos %3d %3d %+.11e %+.11e\n",
	      fprintf(cor_meson_fd, "%2d pos %3d %3d %+.7e %+.7e\n",
		      cor_meson_T, index_bra, index_ket,
		      cor_meson_pos[i2][i3][cor_meson_T].re,
		      cor_meson_pos[i2][i3][cor_meson_T].im);


	      // negative time direction

	      // index_states_1   -->   ket (phi)
	      // index_states_2   -->   bra (xi)

	      // i2   -->   phi (ket)
	      // i3   -->   xi (bra)

	      // fprintf(cor_meson_fd, "%2d neg %3d %3d %+.11e %+.11e\n",
	      fprintf(cor_meson_fd, "%2d neg %3d %3d %+.7e %+.7e\n",
		      cor_meson_T, index_ket, index_bra,
		      cor_meson_neg[i2][i3][cor_meson_T].re,
		      cor_meson_neg[i2][i3][cor_meson_T].im);
	    }
	}
    }

  fprintf(stderr, "\n");


  // **********
  Update_Time("contraction [end]");
  // **********

  // ********************
  // ********************
  // ********************
  // ********************
  // ********************
  // contraction loop (end)
  // ********************
  // ********************
  // ********************
  // ********************
  // ********************


	}
    }


  fclose(cor_meson_fd);


  free(Wilson_lines_pos);
  free(Wilson_lines_neg);

  for(i1 = 0; i1 < phi_num; i1++)
    {
      Spinor_Field_Free(phi + i1);
      Spinor_Field_Free(xi + i1);
      Spinor_Field_Free(xi_org + i1);
    }

  fprintf(stderr, "  ... o.k.\n");


  // **********


  Gauge_Field_Free(&gauge_field);

#ifdef __PRECOMPUTE_WIDTH_6__
  free(meson_paths_6);
#endif

#ifdef __PRECOMPUTE_WIDTH_8__
  free(meson_paths_8);
#endif


  // **********
  fprintf(stderr, "\n");
  Update_Time("contractions finished");
  // **********


  return EXIT_SUCCESS;
}



// ********************



// Six paths meson creation operator.

// psi: spinor field
// (it, ix, iy, iz): spacetime point
// ast: "which side" --> cf. my personal notes.

void ExtendedOperator_6(

  double *PSI__S__gamma_5,
  double *PSI__P_mi__1,
  double *PSI__S__gamma_5_gamma_r,
  double *PSI__P_mi__gamma_r,

  double *PSI__P_pl__gamma_x_x_mi_gamma_y_y,
  double *PSI__P_pl__gamma_y_y_mi_gamma_z_z,
  double *PSI__P_pl__gamma_z_z_mi_gamma_x_x,
  double *PSI__D_mi__gamma_5_gamma_x_x_mi_gamma_y_y,
  double *PSI__D_mi__gamma_5_gamma_y_y_mi_gamma_z_z,
  double *PSI__D_mi__gamma_5_gamma_z_z_mi_gamma_x_x,

  double *psi, int it, int ix, int iy, int iz, bool ast, int operator_width,

  int index___, int delta_t___,  // only if __PRECOMPUTE_WIDTH_6__ is defined

  bool use_slice_0_for_psi)  // for __OPTIMIZE_MEM__
{
  int idx, idy, idz;

  double spinor1[24], spinor2[24];


  if(operator_width < 0 || operator_width > L)
    {
      fprintf(stderr, "Error: void ExtendedOperator_6(...\n");
      exit(EXIT_FAILURE);
    }


  fv_eq_zero(PSI__S__gamma_5);
  fv_eq_zero(PSI__P_mi__1);
  fv_eq_zero(PSI__S__gamma_5_gamma_r);
  fv_eq_zero(PSI__P_mi__gamma_r);

  fv_eq_zero(PSI__P_pl__gamma_x_x_mi_gamma_y_y);
  fv_eq_zero(PSI__P_pl__gamma_y_y_mi_gamma_z_z);
  fv_eq_zero(PSI__P_pl__gamma_z_z_mi_gamma_x_x);
  fv_eq_zero(PSI__D_mi__gamma_5_gamma_x_x_mi_gamma_y_y);
  fv_eq_zero(PSI__D_mi__gamma_5_gamma_y_y_mi_gamma_z_z);
  fv_eq_zero(PSI__D_mi__gamma_5_gamma_z_z_mi_gamma_x_x);


  double U[18], U_[18];


  int it_psi = it;
  int T_lat = T;

  if(use_slice_0_for_psi == true)
    it_psi = 0;

#ifdef __OPTIMIZE_MEM__
  T_lat = 2 * cor_meson_T_max + 1;
#endif


  // --------------------
  // --------------------
  // negative x-direction
  // --------------------
  // --------------------

#ifdef __PRECOMPUTE_WIDTH_6__
  if(index___ == -1)
    {
      cm_eq_id(U);
    }
  else
    {
      cm_eq_cm_dag(U, meson_paths_6 + meson_paths_6_gi(index___, delta_t___, ix-operator_width, iy, iz, 0, T, L));
    }
#else
  cm_eq_id(U);

  for(idx = 0; idx < operator_width; idx++)
    {
      cm_eq_cm(U_, U);

      cm_eq_cm_ti_cm_dag(U, U_,
	gauge_field + ggi(get_index(it, ix-1-idx, iy, iz, T_lat, L), 1));
    }
#endif

  // *****************************************************************************
  // !!!!! Gamma = 1 !!!!!
  // PHI (ast = false)   -->   Gamma = +1
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = +1
  // 
  // - PSI__P_mi__1
  // *****************************************************************************

  fv_pl_eq_cm_ti_fv(PSI__P_mi__1, U,
		    psi + gsi(get_index(it_psi, ix-operator_width, iy, iz, T_lat, L)));

  // *****************************************************************************
  // !!!!! Gamma = gamma_5 !!!!!
  // PHI (ast = false)   -->   Gamma = +gamma_5
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = -gamma_5
  // 
  // - PSI__S__gamma_5
  // *****************************************************************************

  // spinor1 = +/- gamma_5 * psi
  fv_eq_gamma_ti_fv(spinor1, 5,
		    psi + gsi(get_index(it_psi, ix-operator_width, iy, iz, T_lat, L)));
  if(ast == true)
    fv_mi(spinor1);

  fv_pl_eq_cm_ti_fv(PSI__S__gamma_5, U, spinor1);

  if(operator_width != 0)
    {
  // *****************************************************************************
  // !!!!! Gamma = gamma_j !!!!!
  // PHI (ast = false)   -->   Gamma = +gamma_j
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = -gamma_j
  // 
  // - PSI__P_mi__gamma_r
  // - PSI__P_pl__gamma_x_x_mi_gamma_y_y
  // - PSI__P_pl__gamma_y_y_mi_gamma_z_z
  // - PSI__P_pl__gamma_z_z_mi_gamma_x_x
  // *****************************************************************************

  // spinor1 = +/- gamma_j * psi
  fv_eq_gamma_ti_fv(spinor1, 1,
		    psi + gsi(get_index(it_psi, ix-operator_width, iy, iz, T_lat, L)));
  if(ast == true)
    fv_mi(spinor1);

  fv_mi_eq_cm_ti_fv(PSI__P_mi__gamma_r, U, spinor1);
  fv_mi_eq_cm_ti_fv(PSI__P_pl__gamma_x_x_mi_gamma_y_y, U, spinor1);
  fv_pl_eq_cm_ti_fv(PSI__P_pl__gamma_z_z_mi_gamma_x_x, U, spinor1);

  // *****************************************************************************
  // !!!!! Gamma = gamma_5 gamma_j !!!!!
  // PHI (ast = false)   -->   Gamma = +gamma_5 gamma_j
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = +gamma_5 gamma_j
  // 
  // - PSI__S__gamma_5_gamma_r
  // - PSI__D_mi__gamma_5_gamma_x_x_mi_gamma_y_y
  // - PSI__D_mi__gamma_5_gamma_y_y_mi_gamma_z_z
  // - PSI__D_mi__gamma_5_gamma_z_z_mi_gamma_x_x
  // *****************************************************************************

  // spinor1 = gamma_j * psi
  fv_eq_gamma_ti_fv(spinor1, 1,
		    psi + gsi(get_index(it_psi, ix-operator_width, iy, iz, T_lat, L)));

  // spinor2 = gamma_5 * spinor1
  fv_eq_gamma_ti_fv(spinor2, 5, spinor1);

  fv_mi_eq_cm_ti_fv(PSI__S__gamma_5_gamma_r, U, spinor2);
  fv_mi_eq_cm_ti_fv(PSI__D_mi__gamma_5_gamma_x_x_mi_gamma_y_y, U, spinor2);
  fv_pl_eq_cm_ti_fv(PSI__D_mi__gamma_5_gamma_z_z_mi_gamma_x_x, U, spinor2);
    }


  // --------------------
  // --------------------
  // positive x-direction
  // --------------------
  // --------------------

#ifdef __PRECOMPUTE_WIDTH_6__
  if(index___ == -1)
    {
      cm_eq_id(U);
    }
  else
    {
      cm_eq_cm(U, meson_paths_6 + meson_paths_6_gi(index___, delta_t___, ix, iy, iz, 0, T, L));
    }
#else
  cm_eq_id(U);

  for(idx = 0; idx < operator_width; idx++)
    {
      cm_eq_cm(U_, U);

      cm_eq_cm_ti_cm(U, U_,
	gauge_field + ggi(get_index(it, ix+idx, iy, iz, T_lat, L), 1));
    }
#endif

  // *****************************************************************************
  // !!!!! Gamma = 1 !!!!!
  // PHI (ast = false)   -->   Gamma = +1
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = +1
  // 
  // - PSI__P_mi__1
  // *****************************************************************************

  fv_pl_eq_cm_ti_fv(PSI__P_mi__1, U,
		    psi + gsi(get_index(it_psi, ix+operator_width, iy, iz, T_lat, L)));

  // *****************************************************************************
  // !!!!! Gamma = gamma_5 !!!!!
  // PHI (ast = false)   -->   Gamma = +gamma_5
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = -gamma_5
  // 
  // - PSI__S__gamma_5
  // *****************************************************************************

  // spinor1 = +/- gamma_5 * psi
  fv_eq_gamma_ti_fv(spinor1, 5,
		    psi + gsi(get_index(it_psi, ix+operator_width, iy, iz, T_lat, L)));
  if(ast == true)
    fv_mi(spinor1);

  fv_pl_eq_cm_ti_fv(PSI__S__gamma_5, U, spinor1);

  if(operator_width != 0)
    {
  // *****************************************************************************
  // !!!!! Gamma = gamma_j !!!!!
  // PHI (ast = false)   -->   Gamma = +gamma_j
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = -gamma_j
  // 
  // - PSI__P_mi__gamma_r
  // - PSI__P_pl__gamma_x_x_mi_gamma_y_y
  // - PSI__P_pl__gamma_y_y_mi_gamma_z_z
  // - PSI__P_pl__gamma_z_z_mi_gamma_x_x
  // *****************************************************************************

  // spinor1 = +/- gamma_j * psi
  fv_eq_gamma_ti_fv(spinor1, 1,
		    psi + gsi(get_index(it_psi, ix+operator_width, iy, iz, T_lat, L)));
  if(ast == true)
    fv_mi(spinor1);

  fv_pl_eq_cm_ti_fv(PSI__P_mi__gamma_r, U, spinor1);
  fv_pl_eq_cm_ti_fv(PSI__P_pl__gamma_x_x_mi_gamma_y_y, U, spinor1);
  fv_mi_eq_cm_ti_fv(PSI__P_pl__gamma_z_z_mi_gamma_x_x, U, spinor1);

  // *****************************************************************************
  // !!!!! Gamma = gamma_5 gamma_j !!!!!
  // PHI (ast = false)   -->   Gamma = +gamma_5 gamma_j
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = +gamma_5 gamma_j
  // 
  // - PSI__S__gamma_5_gamma_r
  // - PSI__D_mi__gamma_5_gamma_x_x_mi_gamma_y_y
  // - PSI__D_mi__gamma_5_gamma_y_y_mi_gamma_z_z
  // - PSI__D_mi__gamma_5_gamma_z_z_mi_gamma_x_x
  // *****************************************************************************

  // spinor1 = gamma_j * psi
  fv_eq_gamma_ti_fv(spinor1, 1,
		    psi + gsi(get_index(it_psi, ix+operator_width, iy, iz, T_lat, L)));

  // spinor2 = gamma_5 * spinor1
  fv_eq_gamma_ti_fv(spinor2, 5, spinor1);

  fv_pl_eq_cm_ti_fv(PSI__S__gamma_5_gamma_r, U, spinor2);
  fv_pl_eq_cm_ti_fv(PSI__D_mi__gamma_5_gamma_x_x_mi_gamma_y_y, U, spinor2);
  fv_mi_eq_cm_ti_fv(PSI__D_mi__gamma_5_gamma_z_z_mi_gamma_x_x, U, spinor2);
    }


  // --------------------
  // --------------------
  // negative y-direction
  // --------------------
  // --------------------

#ifdef __PRECOMPUTE_WIDTH_6__
  if(index___ == -1)
    {
      cm_eq_id(U);
    }
  else
    {
      cm_eq_cm_dag(U, meson_paths_6 + meson_paths_6_gi(index___, delta_t___, ix, iy-operator_width, iz, 1, T, L));
    }
#else
  cm_eq_id(U);

  for(idy = 0; idy < operator_width; idy++)
    {
      cm_eq_cm(U_, U);

      cm_eq_cm_ti_cm_dag(U, U_,
	gauge_field + ggi(get_index(it, ix, iy-1-idy, iz, T_lat, L), 2));
    }
#endif

  // *****************************************************************************
  // !!!!! Gamma = 1 !!!!!
  // PHI (ast = false)   -->   Gamma = +1
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = +1
  // 
  // - PSI__P_mi__1
  // *****************************************************************************

  fv_pl_eq_cm_ti_fv(PSI__P_mi__1, U,
		    psi + gsi(get_index(it_psi, ix, iy-operator_width, iz, T_lat, L)));

  // *****************************************************************************
  // !!!!! Gamma = gamma_5 !!!!!
  // PHI (ast = false)   -->   Gamma = +gamma_5
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = -gamma_5
  // 
  // - PSI__S__gamma_5
  // *****************************************************************************

  // spinor1 = +/- gamma_5 * psi
  fv_eq_gamma_ti_fv(spinor1, 5,
		    psi + gsi(get_index(it_psi, ix, iy-operator_width, iz, T_lat, L)));
  if(ast == true)
    fv_mi(spinor1);

  fv_pl_eq_cm_ti_fv(PSI__S__gamma_5, U, spinor1);

  if(operator_width != 0)
    {
  // *****************************************************************************
  // !!!!! Gamma = gamma_j !!!!!
  // PHI (ast = false)   -->   Gamma = +gamma_j
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = -gamma_j
  // 
  // - PSI__P_mi__gamma_r
  // - PSI__P_pl__gamma_x_x_mi_gamma_y_y
  // - PSI__P_pl__gamma_y_y_mi_gamma_z_z
  // - PSI__P_pl__gamma_z_z_mi_gamma_x_x
  // *****************************************************************************

  // spinor1 = +/- gamma_j * psi
  fv_eq_gamma_ti_fv(spinor1, 2,
		    psi + gsi(get_index(it_psi, ix, iy-operator_width, iz, T_lat, L)));
  if(ast == true)
    fv_mi(spinor1);

  fv_mi_eq_cm_ti_fv(PSI__P_mi__gamma_r, U, spinor1);
  fv_mi_eq_cm_ti_fv(PSI__P_pl__gamma_y_y_mi_gamma_z_z, U, spinor1);
  fv_pl_eq_cm_ti_fv(PSI__P_pl__gamma_x_x_mi_gamma_y_y, U, spinor1);

  // *****************************************************************************
  // !!!!! Gamma = gamma_5 gamma_j !!!!!
  // PHI (ast = false)   -->   Gamma = +gamma_5 gamma_j
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = +gamma_5 gamma_j
  // 
  // - PSI__S__gamma_5_gamma_r
  // - PSI__D_mi__gamma_5_gamma_x_x_mi_gamma_y_y
  // - PSI__D_mi__gamma_5_gamma_y_y_mi_gamma_z_z
  // - PSI__D_mi__gamma_5_gamma_z_z_mi_gamma_x_x
  // *****************************************************************************

  // spinor1 = gamma_j * psi
  fv_eq_gamma_ti_fv(spinor1, 2,
		    psi + gsi(get_index(it_psi, ix, iy-operator_width, iz, T_lat, L)));

  // spinor2 = gamma_5 * spinor1
  fv_eq_gamma_ti_fv(spinor2, 5, spinor1);

  fv_mi_eq_cm_ti_fv(PSI__S__gamma_5_gamma_r, U, spinor2);
  fv_mi_eq_cm_ti_fv(PSI__D_mi__gamma_5_gamma_y_y_mi_gamma_z_z, U, spinor2);
  fv_pl_eq_cm_ti_fv(PSI__D_mi__gamma_5_gamma_x_x_mi_gamma_y_y, U, spinor2);
    }


  // --------------------
  // --------------------
  // positive y-direction
  // --------------------
  // --------------------

#ifdef __PRECOMPUTE_WIDTH_6__
  if(index___ == -1)
    {
      cm_eq_id(U);
    }
  else
    {
      cm_eq_cm(U, meson_paths_6 + meson_paths_6_gi(index___, delta_t___, ix, iy, iz, 1, T, L));
    }
#else
  cm_eq_id(U);

  for(idy = 0; idy < operator_width; idy++)
    {
      cm_eq_cm(U_, U);

      cm_eq_cm_ti_cm(U, U_,
	gauge_field + ggi(get_index(it, ix, iy+idy, iz, T_lat, L), 2));
    }
#endif

  // *****************************************************************************
  // !!!!! Gamma = 1 !!!!!
  // PHI (ast = false)   -->   Gamma = +1
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = +1
  // 
  // - PSI__P_mi__1
  // *****************************************************************************

  fv_pl_eq_cm_ti_fv(PSI__P_mi__1, U,
		    psi + gsi(get_index(it_psi, ix, iy+operator_width, iz, T_lat, L)));

  // *****************************************************************************
  // !!!!! Gamma = gamma_5 !!!!!
  // PHI (ast = false)   -->   Gamma = +gamma_5
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = -gamma_5
  // 
  // - PSI__S__gamma_5
  // *****************************************************************************

  // spinor1 = +/- gamma_5 * psi
  fv_eq_gamma_ti_fv(spinor1, 5,
		    psi + gsi(get_index(it_psi, ix, iy+operator_width, iz, T_lat, L)));
  if(ast == true)
    fv_mi(spinor1);

  fv_pl_eq_cm_ti_fv(PSI__S__gamma_5, U, spinor1);

  if(operator_width != 0)
    {
  // *****************************************************************************
  // !!!!! Gamma = gamma_j !!!!!
  // PHI (ast = false)   -->   Gamma = +gamma_j
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = -gamma_j
  // 
  // - PSI__P_mi__gamma_r
  // - PSI__P_pl__gamma_x_x_mi_gamma_y_y
  // - PSI__P_pl__gamma_y_y_mi_gamma_z_z
  // - PSI__P_pl__gamma_z_z_mi_gamma_x_x
  // *****************************************************************************

  // spinor1 = +/- gamma_j * psi
  fv_eq_gamma_ti_fv(spinor1, 2,
		    psi + gsi(get_index(it_psi, ix, iy+operator_width, iz, T_lat, L)));
  if(ast == true)
    fv_mi(spinor1);

  fv_pl_eq_cm_ti_fv(PSI__P_mi__gamma_r, U, spinor1);
  fv_pl_eq_cm_ti_fv(PSI__P_pl__gamma_y_y_mi_gamma_z_z, U, spinor1);
  fv_mi_eq_cm_ti_fv(PSI__P_pl__gamma_x_x_mi_gamma_y_y, U, spinor1);

  // *****************************************************************************
  // !!!!! Gamma = gamma_5 gamma_j !!!!!
  // PHI (ast = false)   -->   Gamma = +gamma_5 gamma_j
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = +gamma_5 gamma_j
  // 
  // - PSI__S__gamma_5_gamma_r
  // - PSI__D_mi__gamma_5_gamma_x_x_mi_gamma_y_y
  // - PSI__D_mi__gamma_5_gamma_y_y_mi_gamma_z_z
  // - PSI__D_mi__gamma_5_gamma_z_z_mi_gamma_x_x
  // *****************************************************************************

  // spinor1 = gamma_j * psi
  fv_eq_gamma_ti_fv(spinor1, 2,
		    psi + gsi(get_index(it_psi, ix, iy+operator_width, iz, T_lat, L)));

  // spinor2 = gamma_5 * spinor1
  fv_eq_gamma_ti_fv(spinor2, 5, spinor1);

  fv_pl_eq_cm_ti_fv(PSI__S__gamma_5_gamma_r, U, spinor2);
  fv_pl_eq_cm_ti_fv(PSI__D_mi__gamma_5_gamma_y_y_mi_gamma_z_z, U, spinor2);
  fv_mi_eq_cm_ti_fv(PSI__D_mi__gamma_5_gamma_x_x_mi_gamma_y_y, U, spinor2);
    }


  // --------------------
  // --------------------
  // negative z-direction
  // --------------------
  // --------------------

#ifdef __PRECOMPUTE_WIDTH_6__
  if(index___ == -1)
    {
      cm_eq_id(U);
    }
  else
    {
      cm_eq_cm_dag(U, meson_paths_6 + meson_paths_6_gi(index___, delta_t___, ix, iy, iz-operator_width, 2, T, L));
    }
#else
  cm_eq_id(U);

  for(idz = 0; idz < operator_width; idz++)
    {
      cm_eq_cm(U_, U);

      cm_eq_cm_ti_cm_dag(U, U_,
	gauge_field + ggi(get_index(it, ix, iy, iz-1-idz, T_lat, L), 3));
    }
#endif

  // *****************************************************************************
  // !!!!! Gamma = 1 !!!!!
  // PHI (ast = false)   -->   Gamma = +1
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = +1
  // 
  // - PSI__P_mi__1
  // *****************************************************************************

  fv_pl_eq_cm_ti_fv(PSI__P_mi__1, U,
		    psi + gsi(get_index(it_psi, ix, iy, iz-operator_width, T_lat, L)));

  // *****************************************************************************
  // !!!!! Gamma = gamma_5 !!!!!
  // PHI (ast = false)   -->   Gamma = +gamma_5
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = -gamma_5
  // 
  // - PSI__S__gamma_5
  // *****************************************************************************

  // spinor1 = +/- gamma_5 * psi
  fv_eq_gamma_ti_fv(spinor1, 5,
		    psi + gsi(get_index(it_psi, ix, iy, iz-operator_width, T_lat, L)));
  if(ast == true)
    fv_mi(spinor1);

  fv_pl_eq_cm_ti_fv(PSI__S__gamma_5, U, spinor1);

  if(operator_width != 0)
    {
  // *****************************************************************************
  // !!!!! Gamma = gamma_j !!!!!
  // PHI (ast = false)   -->   Gamma = +gamma_j
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = -gamma_j
  // 
  // - PSI__P_mi__gamma_r
  // - PSI__P_pl__gamma_x_x_mi_gamma_y_y
  // - PSI__P_pl__gamma_y_y_mi_gamma_z_z
  // - PSI__P_pl__gamma_z_z_mi_gamma_x_x
  // *****************************************************************************

  // spinor1 = +/- gamma_j * psi
  fv_eq_gamma_ti_fv(spinor1, 3,
		    psi + gsi(get_index(it_psi, ix, iy, iz-operator_width, T_lat, L)));
  if(ast == true)
    fv_mi(spinor1);

  fv_mi_eq_cm_ti_fv(PSI__P_mi__gamma_r, U, spinor1);
  fv_mi_eq_cm_ti_fv(PSI__P_pl__gamma_z_z_mi_gamma_x_x, U, spinor1);
  fv_pl_eq_cm_ti_fv(PSI__P_pl__gamma_y_y_mi_gamma_z_z, U, spinor1);

  // *****************************************************************************
  // !!!!! Gamma = gamma_5 gamma_j !!!!!
  // PHI (ast = false)   -->   Gamma = +gamma_5 gamma_j
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = +gamma_5 gamma_j
  // 
  // - PSI__S__gamma_5_gamma_r
  // - PSI__D_mi__gamma_5_gamma_x_x_mi_gamma_y_y
  // - PSI__D_mi__gamma_5_gamma_y_y_mi_gamma_z_z
  // - PSI__D_mi__gamma_5_gamma_z_z_mi_gamma_x_x
  // *****************************************************************************

  // spinor1 = gamma_j * psi
  fv_eq_gamma_ti_fv(spinor1, 3,
		    psi + gsi(get_index(it_psi, ix, iy, iz-operator_width, T_lat, L)));

  // spinor2 = gamma_5 * spinor1
  fv_eq_gamma_ti_fv(spinor2, 5, spinor1);

  fv_mi_eq_cm_ti_fv(PSI__S__gamma_5_gamma_r, U, spinor2);
  fv_mi_eq_cm_ti_fv(PSI__D_mi__gamma_5_gamma_z_z_mi_gamma_x_x, U, spinor2);
  fv_pl_eq_cm_ti_fv(PSI__D_mi__gamma_5_gamma_y_y_mi_gamma_z_z, U, spinor2);
    }


  // --------------------
  // --------------------
  // positive z-direction
  // --------------------
  // --------------------

#ifdef __PRECOMPUTE_WIDTH_6__
  if(index___ == -1)
    {
      cm_eq_id(U);
    }
  else
    {
      cm_eq_cm(U, meson_paths_6 + meson_paths_6_gi(index___, delta_t___, ix, iy, iz, 2, T, L));
    }
#else
  cm_eq_id(U);

  for(idz = 0; idz < operator_width; idz++)
    {
      cm_eq_cm(U_, U);

      cm_eq_cm_ti_cm(U, U_,
	gauge_field + ggi(get_index(it, ix, iy, iz+idz, T_lat, L), 3));
    }
#endif

  // *****************************************************************************
  // !!!!! Gamma = 1 !!!!!
  // PHI (ast = false)   -->   Gamma = +1
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = +1
  // 
  // - PSI__P_mi__1
  // *****************************************************************************

  fv_pl_eq_cm_ti_fv(PSI__P_mi__1, U,
		    psi + gsi(get_index(it_psi, ix, iy, iz+operator_width, T_lat, L)));

  // *****************************************************************************
  // !!!!! Gamma = gamma_5 !!!!!
  // PHI (ast = false)   -->   Gamma = +gamma_5
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = -gamma_5
  // 
  // - PSI__S__gamma_5
  // *****************************************************************************

  // spinor1 = +/- gamma_5 * psi
  fv_eq_gamma_ti_fv(spinor1, 5,
		    psi + gsi(get_index(it_psi, ix, iy, iz+operator_width, T_lat, L)));
  if(ast == true)
    fv_mi(spinor1);

  fv_pl_eq_cm_ti_fv(PSI__S__gamma_5, U, spinor1);

  if(operator_width != 0)
    {
  // *****************************************************************************
  // !!!!! Gamma = gamma_j !!!!!
  // PHI (ast = false)   -->   Gamma = +gamma_j
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = -gamma_j
  // 
  // - PSI__P_mi__gamma_r
  // - PSI__P_pl__gamma_x_x_mi_gamma_y_y
  // - PSI__P_pl__gamma_y_y_mi_gamma_z_z
  // - PSI__P_pl__gamma_z_z_mi_gamma_x_x
  // *****************************************************************************

  // spinor1 = +/- gamma_j * psi
  fv_eq_gamma_ti_fv(spinor1, 3,
		    psi + gsi(get_index(it_psi, ix, iy, iz+operator_width, T_lat, L)));
  if(ast == true)
    fv_mi(spinor1);

  fv_pl_eq_cm_ti_fv(PSI__P_mi__gamma_r, U, spinor1);
  fv_pl_eq_cm_ti_fv(PSI__P_pl__gamma_z_z_mi_gamma_x_x, U, spinor1);
  fv_mi_eq_cm_ti_fv(PSI__P_pl__gamma_y_y_mi_gamma_z_z, U, spinor1);

  // *****************************************************************************
  // !!!!! Gamma = gamma_5 gamma_j !!!!!
  // PHI (ast = false)   -->   Gamma = +gamma_5 gamma_j
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = +gamma_5 gamma_j
  // 
  // - PSI__S__gamma_5_gamma_r
  // - PSI__D_mi__gamma_5_gamma_x_x_mi_gamma_y_y
  // - PSI__D_mi__gamma_5_gamma_y_y_mi_gamma_z_z
  // - PSI__D_mi__gamma_5_gamma_z_z_mi_gamma_x_x
  // *****************************************************************************

  // spinor1 = gamma_j * psi
  fv_eq_gamma_ti_fv(spinor1, 3,
		    psi + gsi(get_index(it_psi, ix, iy, iz+operator_width, T_lat, L)));

  // spinor2 = gamma_5 * spinor1
  fv_eq_gamma_ti_fv(spinor2, 5, spinor1);

  fv_pl_eq_cm_ti_fv(PSI__S__gamma_5_gamma_r, U, spinor2);
  fv_pl_eq_cm_ti_fv(PSI__D_mi__gamma_5_gamma_z_z_mi_gamma_x_x, U, spinor2);
  fv_mi_eq_cm_ti_fv(PSI__D_mi__gamma_5_gamma_y_y_mi_gamma_z_z, U, spinor2);
    }
}



// ********************



// Eight paths meson creation operator.

// psi: spinor field
// (it, ix, iy, iz): spacetime point
// ast: "which side" --> cf. my personal notes.

void ExtendedOperator_8(

  double *PSI__S__gamma_5,
  double *PSI__P_mi__1,
  double *PSI__S__gamma_5_gamma_r,
  double *PSI__P_mi__gamma_r,

  double *PSI__P_pl__gamma_x_x_mi_gamma_y_y,
  double *PSI__P_pl__gamma_y_y_mi_gamma_z_z,
  double *PSI__P_pl__gamma_z_z_mi_gamma_x_x,
  double *PSI__D_mi__gamma_5_gamma_x_x_mi_gamma_y_y,
  double *PSI__D_mi__gamma_5_gamma_y_y_mi_gamma_z_z,
  double *PSI__D_mi__gamma_5_gamma_z_z_mi_gamma_x_x,

  double *PSI__D_pl__gamma_x_y_z_pl_gamma_y_z_x_pl_gamma_z_x_y,
  double *PSI__F_mi__gamma_5_gamma_x_y_z_pl_gamma_y_z_x_pl_gamma_z_x_y,

  double *psi, int it, int ix, int iy, int iz, bool ast, int operator_width,

  int index___, int delta_t___,  // only if __PRECOMPUTE_WIDTH_8__ is defined

  bool use_slice_0_for_psi)  // for __OPTIMIZE_MEM__
{
  int idx, idy, idz;

  double spinor1[24], spinor2[24];


  if(operator_width < 0 || operator_width > L)
    {
      fprintf(stderr, "Error: void ExtendedOperator_8(...\n");
      exit(EXIT_FAILURE);
    }


  fv_eq_zero(PSI__S__gamma_5);
  fv_eq_zero(PSI__P_mi__1);
  fv_eq_zero(PSI__S__gamma_5_gamma_r);
  fv_eq_zero(PSI__P_mi__gamma_r);

  fv_eq_zero(PSI__P_pl__gamma_x_x_mi_gamma_y_y);
  fv_eq_zero(PSI__P_pl__gamma_y_y_mi_gamma_z_z);
  fv_eq_zero(PSI__P_pl__gamma_z_z_mi_gamma_x_x);
  fv_eq_zero(PSI__D_mi__gamma_5_gamma_x_x_mi_gamma_y_y);
  fv_eq_zero(PSI__D_mi__gamma_5_gamma_y_y_mi_gamma_z_z);
  fv_eq_zero(PSI__D_mi__gamma_5_gamma_z_z_mi_gamma_x_x);

  fv_eq_zero(PSI__D_pl__gamma_x_y_z_pl_gamma_y_z_x_pl_gamma_z_x_y);
  fv_eq_zero(PSI__F_mi__gamma_5_gamma_x_y_z_pl_gamma_y_z_x_pl_gamma_z_x_y);


  double U[18];


  int it_psi = it;
  int T_lat = T;

  if(use_slice_0_for_psi == true)
    it_psi = 0;

#ifdef __OPTIMIZE_MEM__
  T_lat = 2 * cor_meson_T_max + 1;
#endif


  int dirx, diry, dirz;

  for(dirx = -1; dirx <= +1; dirx+=2)
    {
      for(diry = -1; diry <= +1; diry+=2)
	{
	  for(dirz = -1; dirz <= +1; dirz+=2)
	    {

  // --------------------
  // --------------------
  // (dirx,diry,dirz)-direction
  // --------------------
  // --------------------

#ifdef __PRECOMPUTE_WIDTH_8__
  if(index___ == -1)
    {
      cm_eq_id(U);
    }
  else
    {
      if(dirx == +1)
	{
	  if(diry == +1)
	    {
	      if(dirz == +1)
		{
		  cm_eq_cm(U, meson_paths_8 + meson_paths_8_gi(index___, delta_t___, ix, iy, iz, 3, T, L));
		}
	      else
		{
		  cm_eq_cm(U, meson_paths_8 + meson_paths_8_gi(index___, delta_t___, ix, iy, iz, 2, T, L));
		}
	    }
	  else
	    {
	      if(dirz == +1)
		{
		  cm_eq_cm(U, meson_paths_8 + meson_paths_8_gi(index___, delta_t___, ix, iy, iz, 1, T, L));
		}
	      else
		{
		  cm_eq_cm(U, meson_paths_8 + meson_paths_8_gi(index___, delta_t___, ix, iy, iz, 0, T, L));
		}
	    }
	}
      else
	{
	  if(diry == +1)
	    {
	      if(dirz == +1)
		{
		  cm_eq_cm_dag(U, meson_paths_8 + meson_paths_8_gi(index___, delta_t___, ix-operator_width, iy+operator_width, iz+operator_width, 0, T, L));
		}
	      else
		{
		  cm_eq_cm_dag(U, meson_paths_8 + meson_paths_8_gi(index___, delta_t___, ix-operator_width, iy+operator_width, iz-operator_width, 1, T, L));
		}
	    }
	  else
	    {
	      if(dirz == +1)
		{
		  cm_eq_cm_dag(U, meson_paths_8 + meson_paths_8_gi(index___, delta_t___, ix-operator_width, iy-operator_width, iz+operator_width, 2, T, L));
		}
	      else
		{
		  cm_eq_cm_dag(U, meson_paths_8 + meson_paths_8_gi(index___, delta_t___, ix-operator_width, iy-operator_width, iz-operator_width, 3, T, L));
		}
	    }
	}
    }
#else
  Diagonal_Path_3(U, gauge_field, T, L, it, ix, iy, iz, dirx, diry, dirz, operator_width);
#endif

  // *****************************************************************************
  // !!!!! Gamma = 1 !!!!!
  // PHI (ast = false)   -->   Gamma = +1
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = +1
  // 
  // - PSI__P_mi__1
  // *****************************************************************************

  fv_pl_eq_cm_ti_fv(PSI__P_mi__1, U, psi + gsi(get_index(it_psi, ix + dirx*operator_width, iy + diry*operator_width, iz + dirz*operator_width, T_lat, L)));

  // *****************************************************************************
  // !!!!! Gamma = gamma_5 !!!!!
  // PHI (ast = false)   -->   Gamma = +gamma_5
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = -gamma_5
  // 
  // - PSI__S__gamma_5
  // *****************************************************************************

  // spinor1 = +/- gamma_5 * psi
  fv_eq_gamma_ti_fv(spinor1, 5, psi + gsi(get_index(it_psi, ix + dirx*operator_width, iy + diry*operator_width, iz + dirz*operator_width, T_lat, L)));
  if(ast == true)
    fv_mi(spinor1);

  fv_pl_eq_cm_ti_fv(PSI__S__gamma_5, U, spinor1);

  if(operator_width != 0)
    {
  // *****************************************************************************
  // !!!!! Gamma = gamma_j !!!!!
  // PHI (ast = false)   -->   Gamma = +gamma_j
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = -gamma_j
  // 
  // - PSI__P_mi__gamma_r
  // - PSI__P_pl__gamma_x_x_mi_gamma_y_y
  // - PSI__P_pl__gamma_y_y_mi_gamma_z_z
  // - PSI__P_pl__gamma_z_z_mi_gamma_x_x
  // - PSI__D_pl__gamma_x_y_z_pl_gamma_y_z_x_pl_gamma_z_x_y
  // *****************************************************************************

  // spinor1 = +/- gamma_1 * psi
  fv_eq_gamma_ti_fv(spinor1, 1, psi + gsi(get_index(it_psi, ix + dirx*operator_width, iy + diry*operator_width, iz + dirz*operator_width, T_lat, L)));
  if(ast == true)
    fv_mi(spinor1);

  if(dirx == -1)
    {
      fv_mi_eq_cm_ti_fv(PSI__P_mi__gamma_r, U, spinor1);
      fv_mi_eq_cm_ti_fv(PSI__P_pl__gamma_x_x_mi_gamma_y_y, U, spinor1);
      fv_pl_eq_cm_ti_fv(PSI__P_pl__gamma_z_z_mi_gamma_x_x, U, spinor1);
    }
  else
    {
      fv_pl_eq_cm_ti_fv(PSI__P_mi__gamma_r, U, spinor1);
      fv_pl_eq_cm_ti_fv(PSI__P_pl__gamma_x_x_mi_gamma_y_y, U, spinor1);
      fv_mi_eq_cm_ti_fv(PSI__P_pl__gamma_z_z_mi_gamma_x_x, U, spinor1);
    }

  if(diry*dirz == -1)
    {
      fv_mi_eq_cm_ti_fv(PSI__D_pl__gamma_x_y_z_pl_gamma_y_z_x_pl_gamma_z_x_y, U,
			spinor1);
    }
  else
    {
      fv_pl_eq_cm_ti_fv(PSI__D_pl__gamma_x_y_z_pl_gamma_y_z_x_pl_gamma_z_x_y, U,
			spinor1);
    }

  // spinor1 = +/- gamma_2 * psi
  fv_eq_gamma_ti_fv(spinor1, 2, psi + gsi(get_index(it_psi, ix + dirx*operator_width, iy + diry*operator_width, iz + dirz*operator_width, T_lat, L)));
  if(ast == true)
    fv_mi(spinor1);

  if(diry == -1)
    {
      fv_mi_eq_cm_ti_fv(PSI__P_mi__gamma_r, U, spinor1);
      fv_pl_eq_cm_ti_fv(PSI__P_pl__gamma_x_x_mi_gamma_y_y, U, spinor1);
      fv_mi_eq_cm_ti_fv(PSI__P_pl__gamma_y_y_mi_gamma_z_z, U, spinor1);
    }
  else
    {
      fv_pl_eq_cm_ti_fv(PSI__P_mi__gamma_r, U, spinor1);
      fv_mi_eq_cm_ti_fv(PSI__P_pl__gamma_x_x_mi_gamma_y_y, U, spinor1);
      fv_pl_eq_cm_ti_fv(PSI__P_pl__gamma_y_y_mi_gamma_z_z, U, spinor1);
    }

  if(dirx*dirz == -1)
    {
      fv_mi_eq_cm_ti_fv(PSI__D_pl__gamma_x_y_z_pl_gamma_y_z_x_pl_gamma_z_x_y, U,
			spinor1);
    }
  else
    {
      fv_pl_eq_cm_ti_fv(PSI__D_pl__gamma_x_y_z_pl_gamma_y_z_x_pl_gamma_z_x_y, U,
			spinor1);
    }

  // spinor1 = +/- gamma_3 * psi
  fv_eq_gamma_ti_fv(spinor1, 3, psi + gsi(get_index(it_psi, ix + dirx*operator_width, iy + diry*operator_width, iz + dirz*operator_width, T_lat, L)));
  if(ast == true)
    fv_mi(spinor1);

  if(dirz == -1)
    {
      fv_mi_eq_cm_ti_fv(PSI__P_mi__gamma_r, U, spinor1);
      fv_pl_eq_cm_ti_fv(PSI__P_pl__gamma_y_y_mi_gamma_z_z, U, spinor1);
      fv_mi_eq_cm_ti_fv(PSI__P_pl__gamma_z_z_mi_gamma_x_x, U, spinor1);
    }
  else
    {
      fv_pl_eq_cm_ti_fv(PSI__P_mi__gamma_r, U, spinor1);
      fv_mi_eq_cm_ti_fv(PSI__P_pl__gamma_y_y_mi_gamma_z_z, U, spinor1);
      fv_pl_eq_cm_ti_fv(PSI__P_pl__gamma_z_z_mi_gamma_x_x, U, spinor1);
    }

  if(dirx*diry == -1)
    {
      fv_mi_eq_cm_ti_fv(PSI__D_pl__gamma_x_y_z_pl_gamma_y_z_x_pl_gamma_z_x_y, U,
			spinor1);
    }
  else
    {
      fv_pl_eq_cm_ti_fv(PSI__D_pl__gamma_x_y_z_pl_gamma_y_z_x_pl_gamma_z_x_y, U,
			spinor1);
    }


  // *****************************************************************************
  // !!!!! Gamma = gamma_5 gamma_j !!!!!
  // PHI (ast = false)   -->   Gamma = +gamma_5 gamma_j
  // XI  (ast = true )   -->   gamma_0 Gamma gamma_0 = +gamma_5 gamma_j
  // 
  // - PSI__S__gamma_5_gamma_r
  // - PSI__D_mi__gamma_5_gamma_x_x_mi_gamma_y_y
  // - PSI__D_mi__gamma_5_gamma_y_y_mi_gamma_z_z
  // - PSI__D_mi__gamma_5_gamma_z_z_mi_gamma_x_x
  // - PSI__F_mi__gamma_5_gamma_x_y_z_pl_gamma_y_z_x_pl_gamma_z_x_y
  // *****************************************************************************

  // spinor1 = gamma_1 * psi
  fv_eq_gamma_ti_fv(spinor1, 1, psi + gsi(get_index(it_psi, ix + dirx*operator_width, iy + diry*operator_width, iz + dirz*operator_width, T_lat, L)));

  // spinor2 = gamma_5 * spinor1
  fv_eq_gamma_ti_fv(spinor2, 5, spinor1);

  if(dirx == -1)
    {
      fv_mi_eq_cm_ti_fv(PSI__S__gamma_5_gamma_r, U, spinor2);
      fv_mi_eq_cm_ti_fv(PSI__D_mi__gamma_5_gamma_x_x_mi_gamma_y_y, U, spinor2);
      fv_pl_eq_cm_ti_fv(PSI__D_mi__gamma_5_gamma_z_z_mi_gamma_x_x, U, spinor2);
    }
  else
    {
      fv_pl_eq_cm_ti_fv(PSI__S__gamma_5_gamma_r, U, spinor2);
      fv_pl_eq_cm_ti_fv(PSI__D_mi__gamma_5_gamma_x_x_mi_gamma_y_y, U, spinor2);
      fv_mi_eq_cm_ti_fv(PSI__D_mi__gamma_5_gamma_z_z_mi_gamma_x_x, U, spinor2);
    }

  if(diry*dirz == -1)
    fv_mi_eq_cm_ti_fv(PSI__F_mi__gamma_5_gamma_x_y_z_pl_gamma_y_z_x_pl_gamma_z_x_y, U, spinor2);
  else
    fv_pl_eq_cm_ti_fv(PSI__F_mi__gamma_5_gamma_x_y_z_pl_gamma_y_z_x_pl_gamma_z_x_y, U, spinor2);

  // spinor1 = gamma_2 * psi
  fv_eq_gamma_ti_fv(spinor1, 2, psi + gsi(get_index(it_psi, ix + dirx*operator_width, iy + diry*operator_width, iz + dirz*operator_width, T_lat, L)));

  // spinor2 = gamma_5 * spinor1
  fv_eq_gamma_ti_fv(spinor2, 5, spinor1);

  if(diry == -1)
    {
      fv_mi_eq_cm_ti_fv(PSI__S__gamma_5_gamma_r, U, spinor2);
      fv_pl_eq_cm_ti_fv(PSI__D_mi__gamma_5_gamma_x_x_mi_gamma_y_y, U, spinor2);
      fv_mi_eq_cm_ti_fv(PSI__D_mi__gamma_5_gamma_y_y_mi_gamma_z_z, U, spinor2);
    }
  else
    {
      fv_pl_eq_cm_ti_fv(PSI__S__gamma_5_gamma_r, U, spinor2);
      fv_mi_eq_cm_ti_fv(PSI__D_mi__gamma_5_gamma_x_x_mi_gamma_y_y, U, spinor2);
      fv_pl_eq_cm_ti_fv(PSI__D_mi__gamma_5_gamma_y_y_mi_gamma_z_z, U, spinor2);
    }

  if(dirx*dirz == -1)
    fv_mi_eq_cm_ti_fv(PSI__F_mi__gamma_5_gamma_x_y_z_pl_gamma_y_z_x_pl_gamma_z_x_y, U, spinor2);
  else
    fv_pl_eq_cm_ti_fv(PSI__F_mi__gamma_5_gamma_x_y_z_pl_gamma_y_z_x_pl_gamma_z_x_y, U, spinor2);

  // spinor1 = gamma_3 * psi
  fv_eq_gamma_ti_fv(spinor1, 3, psi + gsi(get_index(it_psi, ix + dirx*operator_width, iy + diry*operator_width, iz + dirz*operator_width, T_lat, L)));

  // spinor2 = gamma_5 * spinor1
  fv_eq_gamma_ti_fv(spinor2, 5, spinor1);

  if(dirz == -1)
    {
      fv_mi_eq_cm_ti_fv(PSI__S__gamma_5_gamma_r, U, spinor2);
      fv_pl_eq_cm_ti_fv(PSI__D_mi__gamma_5_gamma_y_y_mi_gamma_z_z, U, spinor2);
      fv_mi_eq_cm_ti_fv(PSI__D_mi__gamma_5_gamma_z_z_mi_gamma_x_x, U, spinor2);
    }
  else
    {
      fv_pl_eq_cm_ti_fv(PSI__S__gamma_5_gamma_r, U, spinor2);
      fv_mi_eq_cm_ti_fv(PSI__D_mi__gamma_5_gamma_y_y_mi_gamma_z_z, U, spinor2);
      fv_pl_eq_cm_ti_fv(PSI__D_mi__gamma_5_gamma_z_z_mi_gamma_x_x, U, spinor2);
    }

  if(dirx*diry == -1)
    fv_mi_eq_cm_ti_fv(PSI__F_mi__gamma_5_gamma_x_y_z_pl_gamma_y_z_x_pl_gamma_z_x_y, U, spinor2);
  else
    fv_pl_eq_cm_ti_fv(PSI__F_mi__gamma_5_gamma_x_y_z_pl_gamma_y_z_x_pl_gamma_z_x_y, U, spinor2);
    }

  // --------------------
  // --------------------

	    }
	}
    }
}
