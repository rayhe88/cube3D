#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "struct.h"
#include "mathTools.h"


#ifndef _GEOM_DATA_H_
  #define _GEOM_DATA_H_

  #define NPUA1 100

  int getFieldInLine  (double,dataCube, dataRun, const double*, char*);
  int getFieldInPlane (double,dataCube, dataRun, const double*, char*);
  
  double getVecLMN ( double[3],double[3],double[3],double[3],double[3],double[3],double[3],double[3]);

  void sortCoor(double[3],double[3],double[2]);

  void reOrderAtoms(int*,int*,int*,dataCube);

  void origin(double[3],double[3],double[3],double[2],double[3]);

#endif

