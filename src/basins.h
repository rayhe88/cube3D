#include <stdio.h>

#include "lagForm.h"
#include "transU.h"

#define SET_VOID (-5)
#define SET_NULL (-6)

#define TOLERANCE 1.E-10

typedef struct{
    int attr;
    int idx;
    int max;
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

int ascendingGLine(int, int*, double r[3], dataCells*,
                    dataCube, dataRun, double*);

void getHash(dataCells*, int*, int);

void evalBasins   (dataCube cube, dataRun param, double *matU, char* name);

void createArrayCells( int, dataCells **, const char*);

void getCellsPer  (dataCube cube, dataRun param, double *matU);

void getCellsNoPer(dataCube cube, dataRun param, double *matU);

int loadFieldBasins(int,int,int,int,int,double*, double*,double*);

void assignAttr(int,dataCells*, dataCube, int*,dataRun,double*);

int getAttr(int*,double[10][3],double[3],double[3],int,int,int,int,int,dataCube);


void mergeSort(dataCells *, int, int );
void merge    (dataCells *, int, int, int);


