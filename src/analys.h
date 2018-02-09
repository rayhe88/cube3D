#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "struct.h"

#ifndef _ANALYS_H_
 #define _ANALYS_H_

 #define IDX(I,J,K,NPY,NPZ) ( (I)*(NPY)*(NPZ)+(J)*(NPZ)+(K)  )
 #define IDX3(I,J,K)  ( (I)*(9)+(J)*(3)+(K)  )
 //#define IDX4(I,J,K) ( (I)*(16)+(J)*(4)+(K)  )
 //#define IDX5(I,J,K) ( (I)*(25)+(J)*(5)+(K)  )

 int limitsMacro(int,int,int,int,dataCube,double*);
 double fieldNum(double,double,double,double*);
 double fieldNum2(double,double,double,double*);
 int valoresPropios3x3(double*,double*);
 int fieldNCI(double,double,double,double*,double*,double*,double*);

 int newCube1(dataCube,int*,double*,double*,double*,char *);
 int newCube2(dataCube,int*,double*,double*,double*,char *);
 int newCube3(dataCube,int*,double*,double*,double*,double*,char *);


 int coefficients(double,double,double,double,double,double,double*,double*);

 int NumericalCrit01(double,double,double,double*,double*);
 int NumericalCrit02(double,double,double,double*,double*);

#endif
