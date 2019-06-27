// ********************



// generate_color_sources.cc

// Generates a set of three spacetime and color diluted Z2-sources (useful for
// checking gauge invariance).

// Author: Marc Wagner
// Date: December 2007



// ********************



#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "fields.hh"
#include "geometry.hh"
#include "linear_algebra.hh"
#include "propagator_io.hh"



// ********************



// Temporal extension of the lattice.
int T = 4;

// Spatial extension of the lattice.
int L = 4;


// Prefix of the sources.
const char filename_prefix[] = "../conf_2p1p1_L4T4_random_tslice2/point_source.c";


// Spacetime indices of the sources.

int source_it = 2;
int source_ix = 1;
int source_iy = 2;
int source_iz = 3;



// ********************



int IRand(int min, int max)
{
  return min + (int)( ((double)(max-min+1)) * ((double)(rand())) /
		      ((double)(RAND_MAX) + 1.0) );
}



double Random_Z2()
{
  if(IRand(0, 1) == 0)
    return +1.0 / sqrt(2.0);

  return -1.0 / sqrt(2.0);
}



// Generates a source, which is localized in spacetime and in color.
void Generate_Source(double *spinor_field, int source_ic);



// ********************



int main(int argc, char **argv)
{
  int it, ix, iy, iz;
  char string1[1000];


  // **********


  srand(0);
  // srand((unsigned int)time(NULL));


  // **********


  double *spinor_field;

  Spinor_Field_Alloc(&spinor_field, T, L);

  
  Generate_Source(spinor_field, 0);

  sprintf(string1, "%s0", filename_prefix);
  write_lime_spinor(spinor_field, string1, 0, 32, T, L, L, L);


  Generate_Source(spinor_field, 1);

  sprintf(string1, "%s1", filename_prefix);
  write_lime_spinor(spinor_field, string1, 0, 32, T, L, L, L);


  Generate_Source(spinor_field, 2);

  sprintf(string1, "%s2", filename_prefix);
  write_lime_spinor(spinor_field, string1, 0, 32, T, L, L, L);


  Spinor_Field_Free(&spinor_field);

  fprintf(stderr, "Spacetime and color diluted Z2-sources (t=%2d,x=%2d,y=%2d,z=%2d) have been generated (%s?).\n", source_it, source_ix, source_iy, source_iz, filename_prefix);


  // **********


  return EXIT_SUCCESS;
}



// ********************



// Generates a source, which is localized in spacetime and color.

void Generate_Source(double *spinor_field, int source_ic)
{
  int it, ix, iy, iz;

  for(it = 0; it < T; it++)
    {
      for(ix = 0; ix < L; ix++)
	{
	  for(iy = 0; iy < L; iy++)
	    {
	      for(iz = 0; iz < L; iz++)
		{
		  int index = gsi(get_index(it, ix, iy, iz, T, L));

		  fv_eq_zero(spinor_field + index);

		  if(it == source_it &&
		     ix == source_ix &&
		     iy == source_iy &&
		     iz == source_iz)
		    {
		      if(source_ic == 0)
			{
			  spinor_field[index+ 0] = Random_Z2();
			  spinor_field[index+ 1] = Random_Z2();
			  spinor_field[index+ 6] = Random_Z2();
			  spinor_field[index+ 7] = Random_Z2();
			  spinor_field[index+12] = Random_Z2();
			  spinor_field[index+13] = Random_Z2();
			  spinor_field[index+18] = Random_Z2();
			  spinor_field[index+19] = Random_Z2();
			}

		      if(source_ic == 1)
			{
			  spinor_field[index+ 2] = Random_Z2();
			  spinor_field[index+ 3] = Random_Z2();
			  spinor_field[index+ 8] = Random_Z2();
			  spinor_field[index+ 9] = Random_Z2();
			  spinor_field[index+14] = Random_Z2();
			  spinor_field[index+15] = Random_Z2();
			  spinor_field[index+20] = Random_Z2();
			  spinor_field[index+21] = Random_Z2();
			}

		      if(source_ic == 2)
			{
			  spinor_field[index+ 4] = Random_Z2();
			  spinor_field[index+ 5] = Random_Z2();
			  spinor_field[index+10] = Random_Z2();
			  spinor_field[index+11] = Random_Z2();
			  spinor_field[index+16] = Random_Z2();
			  spinor_field[index+17] = Random_Z2();
			  spinor_field[index+22] = Random_Z2();
			  spinor_field[index+23] = Random_Z2();
			}
		    }
		}
	    }
	}
    }
}



// ********************
