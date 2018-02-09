#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <lapacke/lapacke_utils.h>

#include "lagrangePol.h"

#ifndef _PRUEBA_COEF_H_
#define _PRUEBA_COEF_H_


int getCoeff(int poly,double *hvec, double x0,double y0, double z0,double *f, double *c);
int loadField(int poly,int i, int j, int k,int npx, int npy, int npz, double *field, double *fun);

#endif
