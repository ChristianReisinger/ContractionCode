// ********************



// Wilson_loops.hh

// Author: Marc Wagner
// Date: October 2007



// ********************



#ifndef __WILSON_LOOPS_HH__

#define __WILSON_LOOPS_HH__



// ********************



#include "linear_algebra.hh"



// ********************



// Computes a Wilson loop.

// (it,ix,iy,iz): starting point.
// (dir1,dir2): the two directions.
// (ext1, ext2): the two extensions.

complex Wilson_Loop(double *gauge_field, int T, int L, int it, int ix, int iy, int iz, int dir1, int dir2, int ext1, int ext2);

// gauge_field_1: spatial path at T
// gauge_field_2: spatial path at 0

complex Wilson_Loop_(double *gauge_field_1, double *gauge_field_2, int T, int L, int it, int ix, int iy, int iz, int dir2, int ext1, int ext2);



// Computes a 2d-diagonal (off-axis) Wilson loop.

// (it,ix,iy,iz): starting point.
// (dirx,diry,dirz): spatial direction.
// (ext_T,ext_L): extension in temporal and spatial direction.

// Slow version.
complex Wilson_Loop_Diagonal_2(double *gauge_field, int T, int L, int it, int ix, int iy, int iz, int dirx, int diry, int dirz, int ext_T, int ext_L);

// Fast version.
complex Wilson_Loop_Diagonal_2_(double *gauge_field, double *diagonal_links, int T, int L, int it, int ix, int iy, int iz, int dirx, int diry, int dirz, int ext_T, int ext_L);



// Computes a 3d-diagonal (off-axis) Wilson loop.

// (it,ix,iy,iz): starting point.
// (dirx,diry,dirz): spatial direction.
// (ext_T,ext_L): extension in temporal and spatial direction.

// Slow version.
complex Wilson_Loop_Diagonal_3(double *gauge_field, int T, int L, int it, int ix, int iy, int iz, int dirx, int diry, int dirz, int ext_T, int ext_L);

// Fast version.
complex Wilson_Loop_Diagonal_3_(double *gauge_field, double *diagonal_links, int T, int L, int it, int ix, int iy, int iz, int dirx, int diry, int dirz, int ext_T, int ext_L);



// ********************



// Computes a 1d path of spatial links.

// A: the resultung path.
// (it,ix,iy,iz): starting point.
// (dirx,diry,dirz): spatial direction.
// length: length of the path.

void Path_1(double *A, double *gauge_field, int T, int L, int it, int ix, int iy, int iz, int dirx, int diry, int dirz, int length);



// Computes a 2d diagonal path of spatial links.

// A: the resultung diagonal path.
// (it,ix,iy,iz): starting point.
// (dirx,diry,dirz): spatial direction.
// length: length of the path / sqrt(2).

void Diagonal_Path_2(double *A, double *gauge_field, int T, int L, int it, int ix, int iy, int iz, int dirx, int diry, int dirz, int length);



// Computes a 3d diagonal path of spatial links.

// A: the resultung diagonal path.
// (it,ix,iy,iz): starting point.
// (dirx,diry,dirz): spatial direction.
// length: length of the path / sqrt(3).

void Diagonal_Path_3(double *A, double *gauge_field, int T, int L, int it, int ix, int iy, int iz, int dirx, int diry, int dirz, int length);



// ********************



// Computes a 2d diagonal spatial link.

void Diagonal_Link_2(double *A, double *gauge_field, int T, int L, int it, int ix, int iy, int iz, int dirx, int diry, int dirz);

// Computes 2d diagonal spatial links.

void Diagonal_Links_2(double *diagonal_links, double *gauge_field, int T, int L, int dirx, int diry, int dirz);



// Computes a 3d diagonal spatial link.

void Diagonal_Link_3(double *A, double *gauge_field, int T, int L, int it, int ix, int iy, int iz, int dirx, int diry, int dirz);

// Computes 3d diagonal spatial links.

void Diagonal_Links_3(double *diagonal_links, double *gauge_field, int T, int L, int dirx, int diry, int dirz);



// ********************



#endif



// ********************
