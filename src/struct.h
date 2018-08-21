#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef _STRUCT_H_
  #define _STRUCT_H_

  int *zatm;
  double *coor;
  double *field;

  typedef struct{

    int natm;
    int npt;
    int pts[3];
    double min[3];
    double hvec[3];
    double mvec[9];

  } dataCube;


#endif
