#include <stdio.h>

#include "lagForm.h"
#include "transU.h"

#ifndef _BASINS2_H_
  #define _BASINS2_H_
  #define SET_VOID (-5)
  #define SET_NULL (-6)

  #define TOLERANCE 1.E-4

  typedef struct{
      int attr;
      int idx;
      int max;
      double r[3];
      double fun0;
      double fun1;
      double fun2;
  } dataCells2;

void getHash(dataCells2*, int*, int);

void printXYZBasins(int, dataCells2*,dataCube,int*,double*,char*);
double getVolumenCell( double m[9]);

int ascendingGLine(int, double r[3], double hv, double *coor,dataCube,dataRun,double,double*,double*,int,double[3]);
int ascendingGLinePrint(int, double r[3], double hv, double *coor,dataCube,dataRun,double,double*,double*,int,double[3]);

int ascendingGLine2(int, double r[3], double hv, double *coor,dataCube,dataRun,double,double*,double*,int,double[3]);

void evalBasins   (dataCube cube, dataRun param, double *matU,double min0, char* name);

void createArrayCells( int, dataCells2 **, const char*);

void getCells  (dataCube cube, dataRun param, double min, double *matU, char* name);

int loadFieldBasins(int,int,int,int,int,double*, double*,double*);

int assignAttr(int,int,dataCells2*, dataCube, dataRun,double, double*, double*);


int kernelBasinsLoadPer(dataCells2*,dataCube,dataRun,double*,double *matU, double *charge);

int kernelBasinsLoadNoPer(dataCells2*,dataCube,dataRun,double*,double *matU, double *charge);

#endif
