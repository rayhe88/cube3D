#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "file.h"
#include "struct.h"
#include "mathTools.h"

#ifndef _UTILS_H_
  #define _UTILS_H_
  
  #define NCHAR 60

  void printBar      (FILE *);
  void printBanner   (char *text, FILE *);


  void printBar82    (FILE *);
  void printBanner82 (char *text, FILE *);

  void printCube     (char *text, dataCube cube,FILE *);
  void printCubeRot  (int rotate, const char *nameo, dataCube cube);
  void printTime     (char *text);
  void printInfoCube (dataCube cube);
  void printHead     (char *namei, char *namef, char *nameo);
  void printRunning  (dataRun );

  void fieldMinMax   (dataCube, double*,double*);

  void cpyDataCube   (dataCube, dataCube*);
  void printSizeCube (dataCube);

  void printTapas    (dataCube);

  int  sizeUnits     (int,char*);

  void printDataNCI  (int,dataRun, double*,double*, char*);

  void SkewParameters(double *);

  void checkBoundaryCond (dataCube, int *);
#endif