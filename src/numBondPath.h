/**
 * @file   numBondPath.h
 * @brief 
 * @author Raymundo Hernández-Esparza.
 * @date   November 2018.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "struct.h"
#include "findCrit.h"

#ifndef _NUMBONDPATH_H_
 #define _NUMBONDPATH_H_

 #define BPATH_EPS 0.005  //0.02 // Bueno con 0.05
 #define NSTEP 5
 //#define TOL_DIST_ATM 0.35
 #define TOL_DIST_ATM 0.1
 #define MAXPTS 3000

 typedef struct{
    int bcp;
    double r[3];
    double dij;
 } dataBpath;

 int myIsNanInf( double);
 int myIsNanInf_V3(double*);

 void cpyVec3(double[3], double[3]);
 
 void getRiU(double[3], const double*, double[3]);

 void getKnRungeKuta( double[],double[]);

 int bondPath(int,dataCritP*,int,dataCritP*,int*,int*,dataCube,dataRun,double, const double*,char*);

 int perfectCube(int,double *r,double *min, double *max);

 int logFile( int bcp, int rcp, int ccp, int ncp, dataCritP *bondCrit,
              dataCritP *ringCrit, dataCritP *cageCrit, dataCritP *nnucCrit,
              dataCube cube, dataRun param, double min0, int*,int*,const double *matU, char *name);

 int logFileCSV( int bcp,dataCritP *bondCrit,
              dataCube cube, dataRun param, double min0,const double *matU, char *name,char *string);

 int axesCrit(int bcp, int rcp, int ccp, int ncp, dataCritP *bondCrit,
              dataCritP *ringCrit, dataCritP *cageCrit, dataCritP *nnucCrit,
              dataCube cube, dataRun param, double min0, const double *matU, char *name);


 void printLog1( int,int,int,int,int,FILE*);

 void centraMess(char*, FILE*);

 void getEnergies( double,double,double*,double*,double*);

 void printBpaths(int,int,FILE*,FILE*);

 void sortBpaths (int, int, int*, dataBpath*);

 void createArrayDataBpath(int,dataBpath**,char *);

 void mergeBPathBCP(dataBpath*,int, int);
 void sortBPathBCP(dataBpath*,int,int,int);

 void mergeBPathDist(dataBpath*,int, int);
 void sortBPathDist(dataBpath*,int,int,int);



#endif

