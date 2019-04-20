#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "struct.h"

#ifndef _CUBEINDEX_H_
  #define _CUBEINDEX_H_

  int getLeftRight (int pol, int *min, int *max);
  int getIndex     (int ptq, double q, double q0, double hq);
  int getIndex3D   (int *pts, double *r, double *min, double *h, int *index);

  int getPeriodicIndex( int p, int npx);

  int checkIndex(int *index, int *vecn,dataRun);

  int loadLocalField(int *, dataCube, dataRun, double *, double *,double *, double *);


#endif
