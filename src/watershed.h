#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "struct.h"

#ifndef _WATERSHED_H_
  #define _WATERSHED_H_

  void sortField(dataCube,dataRun,int*,double*);

  int waterShed(dataCube,dataRun,const double*,int*,int*,double*);
  int asignaCenters( int, int*,int*,int[],int[],double[],double[],double[]);



#endif

