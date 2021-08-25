#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mathTools.h"
#include "struct.h"

#ifndef _GEOM_DATA_H_
#define _GEOM_DATA_H_

#define NPUA0 75 // Number of points per atomic unit in 1D
#define NPUA1 50 // Number of points per atomic unit in 2D
#define NPUA2 8  // Number of points per atomic unit in 2D vec
#define NANG 36  // Number of points in the angular part
#define MAXITER 2000
#define EPS 0.0125

typedef struct {
    int i;
    int j;
    int sline;
    double f;
    double r[3];
    double s[2];
} dataSLine;

int getFieldInLine(double, dataCube, dataRun, const double *, char *);

int getFieldInPlane(double, dataCube, dataRun, const double *, char *);

int getGradLinesInPlane1(double, dataCube, dataRun, const double *, char *);
int getGradVectorsInPlane(double, dataCube, dataRun, const double *, char *);
int getStreamLinesInPlane(double, dataCube, dataRun, const double *, char *);

double getVecLMN(double[3], double[3], double[3], double[3], double[3],
                 double[3], double[3], double[3]);

void sortCoor(double[3], double[3], double[2]);

void reOrderAtoms(int *, int *, int *, dataCube, const double *);

void origin(double[3], double[3], double[3], double[2], double[3]);

void getCircumcentre(double[3], double[3], double[3], double[3]);
void getIncentre(double[3], double[3], double[3], double[3]);

void transformCube1(dataCube *, const double *);

void printInfoCube1(dataCube);

int firstSeed(int, dataSLine *);

int getStreamLine(int, int, int, double, dataCube, dataRun, const double *,
                  const double *, double[4][3], double[2], double[2],
                  dataSLine *, FILE *, FILE *);

int check(double, double, double[2], double[2], double[10]);
int checkMatrix(int, dataSLine *, double, double, int);

void printMatrixGeom(int, dataSLine *);

int checkFullMatrix(int, dataSLine *);

int findNeighbours(int, int, dataSLine *);

int getIdfromMatrix(int, int, int, dataSLine *);

void trans3Dto2D(double[3], double[3], double[3], double, double, double *,
                 double *);

void trans2Dto3D(double, double, double[4][3], double[3]);

double RungeKutta4(double[3], double[3], dataCube, dataRun, const double *,
                   const double *, double, int);
double RungeKutta4_fast(double[3], double[3], dataCube, dataRun, const double *,
                        const double *, double, int);
double EulerMethod(double[3], double[3], dataCube, dataRun, const double *,
                   const double *, double, int);

#endif
