#include <stdio.h>

#include "lagForm.h"
#include "transU.h"

#define SET_VOID (-5)
#define SET_NULL (-6)

#define TOLERANCE 1.E-2

typedef struct{
    int attr;
    int idx;
    double r[3];
    double fun0;
    double fun1;
    double fun2;
    double dx;
    double dy;
    double dz;
    double norm;
} dataCells;


//void copyCells    (dataCells , dataCells*);

void evalBasins   (dataCube cube, dataRun param, double *matU, char* name);

void createArrayCells( int, dataCells **, const char*);

int screening(int, dataCells *, dataCells*);

void getCellsPer  (dataCube cube, dataRun param, double *matU);

void getCellsNoPer(dataCube cube, dataRun param, double *matU);

void loadFieldBasins(int,int,int,int,int,double*, double*,double*);

void assignAttr(int,dataCells*, dataCube);

int getAttr(int,int,double, double[],double *);


void mergeSort(dataCells *, int, int );
void merge    (dataCells *, int, int, int);


