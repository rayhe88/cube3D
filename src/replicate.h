#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "struct.h"


#ifndef _REPLICATE_H_
  #define _REPLICATE_H_

  int replicate(int*,dataCube, double[], const double *, const double *, const char*);

  int replicateCoor(int, dataCube,int*,double*,double*, const double*,const double*);

  int getIndexOld(int,int,int[],int[]);

#endif
