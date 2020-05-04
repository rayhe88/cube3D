/**
  @file  mathTools.h
  @brief Set of functions to characterize critical points and to determine 
  the bond paths.
 
  These functions enters a 3x3 matrix and is obtained only the eigenvalues,
  or eigenvalues and eigenvectors to determine the bond paths.
  
  @author Raymundo Hernández-Esparza.
  @date   August 2018.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef _MATHTOOLS_H_
 #define _MATHTOOLS_H_

 #define N 3       /**< Size of the matrices to be diagonalized in the Jacobi function.*/
 #define TOL 1E-09 /**< Tolerance in the diagonalization by the method of Jacobi.*/
 #define NMAX 10   /**< Max number of iterations in the method of Jacobi.*/

 #define ND 20
 #define GRAD(X)  ( (X)*180. / M_PI)

 double determinant3 (double matA[3][3]);
 double mayor( double matAA[][N], int *p, int *q);
 int jacobi(double matAA[][N], double valores[], double eigenvectors[][N]);
 int valoresPropios3x3(double *matA, double *val);
 int valoresPropios3x3_v0(double *matA, double *val);

 int eValNci(double *val);
 int  eigenVV(double mat[9],double val[3],double vec[9]);

 double distance(double r1[3],double r2[3]);
 double distance2(double x,double y, double z,double r2[3]);

 double detMat(double mat[9]);
 double getNorm(double,double,double);
 double getNormVec(double vec[3]);
 double dotProduct(double vecA[3],double vecB[3]);

 void crossProduct(double vecA[3],double vecB[3],double vecOut[3]);
 void matVecProduct(double v[3],double m[9],double vecOut[3]);

 ////void transform(dataCube cube1,dataCube *cube2,double *coor,double *coor2,double matTinv[9]);

#endif

