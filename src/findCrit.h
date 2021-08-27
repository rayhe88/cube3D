/**
 * @file   critical.h
 * @brief
 * @author Raymundo Hernández-Esparza
 * @date   August 2018.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "struct.h"

#ifndef _FINDCRIT_H_
#define _FINDCRIT_H_

//#define TOLFUN0  1.E-6 //original
//#define TOLFUN  1.E-5
//#define TOLFUN0 1.E-7 // original /1.E-7
//#define TOLFUN 1.E-7
//#define TOLGRD 1.E-5 // original
//#define TOLNRM 100.
//#define MAXITER1 30
//#define MAXITER2 1000
//#define PERCENT 1.5
#define ALPHA(X) (1.)
//#define TOLDIST 0.5E-1 //original
//#define TOLDIST 1.72E-1
//#define TOLDIST2  0.2 // Distance for hidrogens 0.06
//#define TOLDIST3  0.7 // Si hay pseudopotenciales esta
// distancia puede servir para
// quitarlos 0.66 u.a. == 0.35 A.
// 0.7 u.a. == 0.3704240743 A

#define XCP -5
#define NCP -4
#define BCP -3
#define RCP -2
#define CCP -1

typedef struct {
    int typ;
    int ran;
    int sig;
    double x;
    double y;
    double z;
    double fun;
} dataCritP;

int delRepCoor(int, double *, double *, dataRC);

int delCoorAtomic(int, double *, const double *, dataCube, dataRC);

int delCoorPseudo(int, double *, double, const double *, dataCube, dataRun,
                  dataRC);

int getTypeCrit(int, int);

void getSignature(double *, int *, int *);

int critPoints(dataCube, dataRun, dataRC, const double *, double, char *);

int rejectCube(double *, double, dataCube, const double *, dataRun, dataRC,
               double *, double *, double *, double *, double *, double *);

int refineCrit(int, double, dataCube, dataRun, dataRC, double *, const double *,
               char *);

int describeCrit(int, double, dataCube, dataRun, dataRC, double *,
                 const double *, char *);

int numCritical01(double, double, double, dataCube, dataRun, const double *,
                  double, double *);

int numCritical02(double, double, double, dataCube, dataRun, const double *,
                  double, double *);

int numCritical01Vec(double[3], dataCube, dataRun, const double *, double,
                     double *);
int numCritical02Vec(double[3], dataCube, dataRun, const double *, double,
                     double *);

int numCritical02Vec_exp(double[3], dataCube, dataRun, const double *, double,
                         double *);

double hessGrad(double *x, double *y, double *z, double *ngrad, double *val);

double getFunInCube(int i, int j, int k, int n1, int n2, double min, double *h,
                    double *field, double *grad);

int getLimits(double *, double *, double *, double);
int inCube(double, double, double, double[]);
int inMacroCube(double, double, double, double *, double *);

int getMaximumInRho(dataCube, double *coor2);

int createArrayCritP(int, dataCritP **, const char *);

int sortCritPoints(int, dataCritP *);

#endif
