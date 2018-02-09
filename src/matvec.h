#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lectura.h"
#include "struct.h"
#include "file.h"
#include "array.h"

#ifndef _MATVEC_H_
 #define _MATVEC_H_

 #define ND 20
 #define GRAD(X)  ( (X)*180. / M_PI)

// void barra(int,int);


 int  eigenVV(double mat[9],double val[3],double vec[9]);

 double detMat(double mat[9]);
 double getNormVec(double vec[3]);
 double dotProduct(double vecA[3],double vecB[3]);

 void getMatT(double *mvec, double *matT,double *matTinv);
 void crossProduct(double vecA[3],double vecB[3],double vecOut[3]);
 void matVecProduct(double v[3],double m[9],double vecOut[3]);

 void transform(dataCube cube1,dataCube *cube2,double *coor,double *coor2,double matTinv[9]);

#endif

