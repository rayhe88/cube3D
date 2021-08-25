/**
  @file  transU.h
  @brief
  @author Raymundo Hernández-Esparza.
  @date   August 2018.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mathTools.h"
#include "struct.h"
#include "utils.h"

#ifndef _TRANSU_H_
#define _TRANSU_H_

void getMatT(dataCube cube, int *rec, double *matT);

void getMatInv(const double *matT, double *matU);

void transformCube(dataCube inp, dataCube *out, int **zatm2, double **coor2,
                   double **field2, double *matU);

void itrans00(double *vec, const double *matU);
void trans00(double *vec, const double *matU);
void trans01(double *vec, const double *matU);
void trans02(double *vec, const double *matU);

#endif
