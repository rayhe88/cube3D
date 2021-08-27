/**
 * @file   numGradLines.h
 * @brief
 * @author Raymundo Hernandez-Esparza.
 * @date   April 2020.
 */

#include "findCrit.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "struct.h"

#ifndef _NUMGRADLINES_H_
#define _NUMGRADLINES_H_

//#define STEP_EPS 0.005

int gradientLines(int, dataCritP *, dataCube, dataRun, double, const double *,
                  char *);
int gradientLines2(int, dataCritP *, dataCube, dataRun, double, const double *,
                   char *);
int gradientLines3(int, int, int, dataCritP *, dataCritP *, dataCritP *,
                   dataCube, dataRun, double, const double *, char *);

int gradientLines4(int, int, int, dataCritP *, dataCritP *, dataCritP *,
                   dataCube, dataRun, double, const double *, char *);
int sortBondCritP(int, double[3], dataCritP *);
void cpyCritP(int, dataCritP *, dataCritP *);
#endif
