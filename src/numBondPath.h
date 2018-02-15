#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifndef _NUMBONDPATH_H_
 #define _NUMBONDPATH_H_

 int myIsNanInf(double);
 int myIsNanInf_V3(double*);

 int logFile(int,int,int,int*,int*,double*,double*,int*,double*,double*,double*,double*,const char *);
 int bondPath(int,int,int,int,int*,int*,double*,double*,int*,double*,double*,double*,double*,const char *);

 int SortCoordinates(int,int*,double*);

 int getData(int,int*,int*,double,double,double,double*,double*,double*,double*);
 int getData2(int,int*,int*,double,double,double,double*,double*,double*,double*);


 void centraMess(char*,FILE*);
  
 void lineDiv(FILE*);

 void printLog1(FILE*,int*);

 void getEnergies(double rho,double lap,double *kinG,double *kinK, double *virial);

 void sortJacobi(double vec[3]);




#endif

