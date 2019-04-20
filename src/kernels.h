#include <stdio.h>
#include <stdlib.h>

#include "file.h"
#include "array.h"
#include "utils.h"
#include "fields.h"
#include "struct.h"

#ifndef _KERNELS_H_
  #define _KERNELS_H_

  void selectExec (dataCube, dataRun, double *, char*);
  int evalRGL     (dataCube, dataRun, double *, char*);
  int evalNCI     (dataCube, dataRun, double *, char*);
  int evalVoidVol (dataCube, dataRun, char*);
  int evalRepCube (dataCube, dataRun, char*);

  double getDenInCube (int,int,int,int,int,double*);

  void getLogField( dataCube, double);

#endif
