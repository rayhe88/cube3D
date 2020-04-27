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
 #define NSTEP 2
 //#define TOL_DIST_ATM 0.35
 #define TOL_DIST_ATM 0.1
 #define MAXPTS 4000

 int myIsNanInf( double);
 int myIsNanInf_V3(double*);

 void cpyVec3(double[3], double[3]);
 
 void getRiU(double[3], const double*, double[3]);

 void getKnRungeKuta( double[],double[]);

 int bondPath(int,dataCritP*,int,dataCritP*,int*,dataCube,dataRun,double, const double*,char*);

 int perfectCube(int,double *r,double *min, double *max);

 int logFile( int bcp, int rcp, int ccp, int ncp, dataCritP *bondCrit,
              dataCritP *ringCrit, dataCritP *cageCrit, dataCritP *nnucCrit,
              dataCube cube, dataRun param, double min0, int*,const double *matU, char *name);

 int logFileCSV( int bcp,dataCritP *bondCrit,
              dataCube cube, dataRun param, double min0,const double *matU, char *name);

 int axesCrit(int bcp, int rcp, int ccp, int ncp, dataCritP *bondCrit,
              dataCritP *ringCrit, dataCritP *cageCrit, dataCritP *nnucCrit,
              dataCube cube, dataRun param, double min0, const double *matU, char *name);


 void printLog1( int,int,int,int,int,FILE*);

 void centraMess(char*, FILE*);

 void getEnergies( double,double,double*,double*,double*);




#endif

