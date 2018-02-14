#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "struct.h"

#ifndef _CRITICAL_H_
 #define _CRITICAL_H_

 #define TOLFUN 0.5E-3
 //#define TOLGRD 1.E-11
 #define TOLGRD 1.E-13
 #define TOLGRD2 1.E-14
 //#define TOLGRD2 1.E-12
 //#define TOLNRM 1.E-7
 #define TOLNRM 1.E-9
 //#define TOLDIS 0.075
 #define TOLDIS 0.001

 #define MAXITER 30

 #define CUTRHO 0.005

 int getIndex(char,int,double,double,double);

 int deleteRepeated(int, double*,double*);

 int getCube(int*,double,double,double,double*,double*,int*,int*,int*);

 int signature(double*,int*,int*);

 int findCriticalPoints(int,dataCube,int*,double*,double*,double*,const char *);

 int findCritical(int,int,int,int*,int*,double,double,double,double,double,double,double*,double*);

 int cubeFail(int,int,int,int,int,int*,double,double,double,double,double*,double*,double*);

 int refineCritical(dataCube,int,int*,double*,double*,double*,double*,FILE*,const char*);

 int inCube(double,double,double,double*);

 double getNorm(double,double,double);

 double mayor3(double,double,double);

 double distance(double,double,double,double,double,double);

 double hessGrad(double*,double*,double*,double*,double*);




#endif
