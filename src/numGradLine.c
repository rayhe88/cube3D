/**
 * @file  numGradLine.c
 * @brief
 * @author Raymundo Hernández-Esparza.
 * @date   December 2018.
 */

#include "file.h"
#include "jacobi.h"
#include "findCrit.h"
#include "mathTools.h"
#include "lagrange2.h"
#include "utils.h"
#include "transU.h"

#include <omp.h>

int ringGradLine( int rcp, dataCritP *ringCrit, int bcp, dataCritP *bondCrit, dataCube cube, dataRun param, double min0, const double *matU, char *name){

  char nameOut[128],tmpname[128];

  FILE *tmp;
  FILE *out;
  
  sprintf(nameOut,"%sRPath.xyz",name);
  tmpFile(&tmp,".c3dRpath_",tmpname,"w+");

  for(i=0;i<rcp;i++){
    x0 = ringCrit[i].x;
    y0 = ringCrit[i].y;
    z0 = ringCrit[i].z;

    numCritical02(x0,y0,z0,cube,param,matU,min0,val);

    matH[0] = val[4]; matH[1] = val[7]; matH[2] = val[8];
    matH[3] = val[7]; matH[4] = val[5]; matH[5] = val[9];
    matH[6] = val[8]; matH[7] = val[9]; matH[8] = val[6];

    JacobiNxN(matH,eval,evec);
    
    vecU[0] = evec[3]; vecU[1] = evec[4]; vecU[2] = evec[5];
    vecV[0] = evec[6]; vecV[1] = evec[7]; vecV[2] = evec[8];

  }


}
