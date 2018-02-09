#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke/lapacke_utils.h>


#ifndef _LAGRANGEPOL_H_
  #define _LAGRANGEPOL_H_

  int evalPol0(int, double,double,double,double*,double*);
  int evalPol1(int, double,double,double,double*,double*);
  int evalPol2(int, double,double,double,double*,double*);

#endif

