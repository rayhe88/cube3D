/**
  @file mathTools.h
  @brief Set of functions to characterize critical points and to determine 
  the bond paths.
 
  These functions enters a 3x3 matrix and is obtained only the eigenvalues,
  or eigenvalues and eigenvectors to determine the bond paths.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef _MATHTOOLS_H_
 #define _MATHTOOLS_H_

 #define N 3       /**< Size of the matrices to be diagonalized in the Jacobi function.*/
 #define TOL 1E-09 /**< Tolerance in the diagonalization by the method of Jacobi.*/
 #define NMAX 10   /**< Max number of iterations in the method of Jacobi.*/

 double mayor( double matAA[][N], int *p, int *q);
 int jacobi(double matAA[][N], double valores[], double eigenvectors[][N]);

#endif
