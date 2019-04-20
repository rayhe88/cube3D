/**
  @file  transU.h
  @brief 
  @author Raymundo Hernández-Esparza.
  @date   August 2018.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "struct.h"
#include "mathTools.h"

#ifndef _TRANSU_H_
 #define _TRANSU_H_

 void getMatT(dataCube cube, int *rec, double *matT);

 void getMatInv(double *matT, double *matU);

 void transformCube (dataCube inp, dataCube *out,
                    int **zatm2, double **coor2,
                    double **field2, double* matU);

 void itrans00( double *vec, const double *matU);
 void trans00 ( double *vec, const double *matU);
 void trans01 ( double *vec, const double *matU);
 void trans02 ( double *vec, const double *matU);

#endif

