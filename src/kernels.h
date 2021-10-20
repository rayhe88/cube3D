#include <stdio.h>
#include <stdlib.h>

#include "array.h"
#include "fields.h"
#include "file.h"
#include "struct.h"
#include "utils.h"

#ifndef _KERNELS_H_
#define _KERNELS_H_

typedef struct {
    int stat;
    int idx;
    int idy;
    int idz;
    double x;
    double y;
    double z;
} dataVoids;

void selectExec(dataCube, dataRun, dataRC, double *, char *);
int evalRGL(dataCube, dataRun, double *, char *);
int evalNCI(dataCube, dataRun, double *, char *);
int evalVoidVol(dataCube, dataRun, const double *, char *);
int evalRepCube(dataCube, dataRun, char *);

double getDenInCube(int, int, int, int, int, double *);

void getLogField(dataCube, double);

int getCubeIdx(int, int, int, int, int, int);

#endif
