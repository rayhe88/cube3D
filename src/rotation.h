#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "struct.h"

#ifndef _ROTATION_H_
  #define _ROTATION_H_

  int checkRotation(double *mvec);

  void rotationCube(double *mvec, dataCube* cube);

  void getAngles( double* vec, double *theta);
  void transform( double *in, double *theta, double *out);
  void transformInv( double *in, double *theta, double *out);

  double getTheta( int i, double *vecV, double *vecU );
  void rotInX ( double *v, double theta, double *out);
  void rotInY ( double *v, double theta, double *out);
  void rotInZ ( double *v, double theta, double *out);

// 1) getAngles(vecOriginales,angles);
// 2) transform(vecIn,angles,vecOut);


#endif
