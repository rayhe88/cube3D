/**
 * @file   fields.h
 * @brief
 * @author Raymundo Hernández-Esparza.
 * @date   August 2018.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "array.h"
#include "cubeIndex.h"
#include "lagForm.h"
#include "lagrange2.h"
#include "mathTools.h"
#include "struct.h"
#include "transU.h"

#ifndef _FIELDS_H_
#define _FIELDS_H_

#define CF 0.16162045967399548133
#define FOT -1.33333333333333333333

int myTernary(int, int, int, int);

double getGrd(double *val);
double getRed(double *val);
double getLap(double *val);
double getKin(double *val);
double getKEW(double *val);
double getVir(double *val);

double campo(int *, int, int, int, int, int, int, double *, double *);

int getFieldPer(dataCube cube, dataRun param, const double *matU,
                double *field2);
int getFieldNoPer(dataCube cube, dataRun param, const double *matU,
                  double *field2);

double campoPer(int *index, int *vecn, dataRun param, double *h, double *field,
                const double *matU,
                int (*f)(double, double, double, double *, double *),
                double (*g)(double *));

double campoNoPer(double x, double y, double z, int *index, int *vecn,
                  dataRun param, double *h, double *field, const double *matU,
                  double (*)(double *), double *min);

int gradientVec(int *index, int *vecn, dataRun param, double *h, double *field,
                const double *matU,
                int (*f)(double, double, double, double *, double *), double[]);

int NCIPer(int *index, int *vecn, dataRun param, double *h, double *field,
           const double *matU, double *ret);

int NCINoPer(double x, double y, double z, int *index, int *vecn, dataRun param,
             double *h, double *field, const double *matU, double *min,
             double *ret);
int getNCIPer(dataCube cube, dataRun param, const double *matU, double *red1,
              double *rho2);
int getNCINoPer(dataCube cube, dataRun param, const double *matU, double *red1,
                double *rho2);
#endif
