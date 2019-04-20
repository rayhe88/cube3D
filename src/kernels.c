/**
 * @file   kernels.c
 * @brief
 * @author Raymundo Hernández-Esparza
 * @date   August 2018.
 */
#include "utils.h"
#include "graph.h"
#include "kernels.h"
#include "findCrit.h"
#include "mathTools.h"
#include "geomData.h"
#include "replicate.h"

/**
 * @brief 
 * @param
 * @param
 * @param
 * @param
 */
void selectExec(dataCube cube, dataRun param, double *matU, char *name){

  double min, max;
  double min2;
  int planeLine = 0;
  min = 0.;


  if( param.task == RED || param.task == GRA ||
      param.task == LAP || param.task == KIN ||
      param.task == VIR ){
      planeLine = 1;
  }

  if( param.geotask == LIN && planeLine == 1){
    fieldMinMax (cube,&min,&max);
    min2 = min - DELTA;
    getLogField (cube,min2);
    getFieldInLine(min,cube,param,matU,name);
  }

  if( param.geotask == PLA && planeLine == 1){
    fieldMinMax (cube,&min,&max);
    min2 = min - DELTA;
    getLogField (cube,min2);
    getFieldInPlane(min,cube,param,matU,name);
  }

  
  


  if( param.geotask != LIN && param.geotask != PLA || planeLine == 0 ){
    switch(param.task){
        case RED:
        case GRA:
        case LAP: 
        case KIN:
        case VIR: evalRGL     (cube,param,matU,name); break;
        case NCI: evalNCI     (cube,param,matU,name); break;
        case CRI: 
                  fieldMinMax (cube,&min,&max);
                  min2 = min - DELTA;
                  getLogField (cube,min2);
                  critPoints  (cube,param,matU,min,name); break;
        case VOI: evalVoidVol (cube,param,name);      break;
        case REP: evalRepCube (cube,param,name);      break;
      }
    }

  printBar(stdout);

}

/**
 * @brief 
 * @param
 * @param
 * @param
 * @param
 */
int evalRGL (dataCube cube, dataRun param, double *matU, char *name){
 
  int i;
  char text[64];
  char nameOut[128];
  double *field2;

  FILE *out;
  dataCube cube2;

  createArrayDou(cube.npt,&field2,"Field 2");

  for(i=0;i<cube.npt;i++)
    field2[i] = (double) 0;

  strcpy(nameOut,name);

  switch( param.task ){

    case RED: strcpy(text,"Reduced Gradient of field"); 
              strcat(nameOut,"Red.cube"); break;

    case GRA: strcpy(text,"Gradient of field"); 
              strcat(nameOut,"Grd.cube"); break;

    case LAP: strcpy(text,"Laplacian of field"); 
              strcat(nameOut,"Lap.cube"); break;

    case KIN: strcpy(text,"Abramov kinetic energy"); 
              strcat(nameOut,"Kin.cube"); break;

    case VIR: strcpy(text,"Virial field using Abramov kinetic energy"); 
              strcat(nameOut,"Vir.cube"); break;

  }


  //llamamos al field
  if( param.pbc == YES )
    getFieldPer  (cube,param,matU,field2);
  else
    getFieldNoPer(cube,param,matU,field2);
  
  openFile(&out,nameOut,"w+");

  cpyDataCube(cube,&cube2);
  cube2.field = NULL;
  cube2.field = field2;
  
  printCube(text,cube2,out);

  printBar(stdout);
  printf("  File %s was generated\n",nameOut);

  fclose(out);
  free(field2);

  return 0;
}

/**
 * @brief 
 * @param
 * @param
 * @param
 * @param
 */
int evalNCI (dataCube cube, dataRun param, double *matU, char *name){
 
  int i;
  char text1[64],text2[64];
  char nameOut1[128],nameOut2[128];
  double *fred1;
  double *frho2;

  FILE *out;

  dataCube cubeRed1;
  dataCube cubeRho2;

  createArrayDou(cube.npt,&fred1,"Field 1");
  createArrayDou(cube.npt,&frho2,"Field 2");

  for(i=0;i<cube.npt;i++){
    fred1[i] = (double) 0;
    frho2[i] = (double) 0;
  }

  strcpy(nameOut1,name);
  strcpy(nameOut2,name);

  strcpy(text1,"Reduced Gradient of field for NCI-Index");
  strcpy(text2," 100.*lambda_2*rho for NCI-Index");

  strcat(nameOut1,"NCIRed.cube");
  strcat(nameOut2,"NCIRho.cube");

  if( param.pbc == YES )
    getNCIPer  (cube,param,matU,fred1,frho2);
  else
    getNCINoPer(cube,param,matU,fred1,frho2);


  printDataNCI(cube.npt,param,fred1,frho2,name);

  openFile(&out,nameOut1,"w+");
  cpyDataCube(cube,&cubeRed1);
  cubeRed1.field = fred1;
  printCube(text1,cubeRed1,out);
  fclose(out);
  free(fred1);

  openFile(&out,nameOut2,"w+");
  cpyDataCube(cube,&cubeRho2);
  cubeRho2.field = NULL;
  cubeRho2.field = frho2;
  printCube(text2,cubeRho2,out);
  fclose(out);
  free(frho2);

  printBar(stdout);
  printf("  File %s was generated\n",nameOut1);
  printBar(stdout);
  printf("  FILE %s was generated\n",nameOut2);
  printBar(stdout);
  printf("  File %s.dat was generated\n",name);
  printBar(stdout);

#ifdef GRAPH
  plotNCIdat(name,param);
#endif
  getNCIgnu(name,param);
  printf("  File plot_NCI_%s.gp was generated\n",name);


  return 0;

}

int evalRepCube ( dataCube cube, dataRun param, char *name){

  char nameOut[128];
  sprintf(nameOut,"%s_%d%d%d.cube",name, 
          param.rep[0], param.rep[1], param.rep[2]);


  int tamanio = param.rep[0]*param.rep[1]*param.rep[2];
  printBanner(" Cube Info Rep ",stdout);
  printf(" Number of atoms      : % 10d\n",(cube.natm)*tamanio);
  printf(" Total field points   : % 10d\n",(cube.npt)*tamanio);
  printf(" Points in X Y Z axes : % 10d % 10d % 10d\n",
          (cube.pts[0])*(param.rep[0]) ,
          (cube.pts[1])*(param.rep[1]) ,
          (cube.pts[2])*(param.rep[2]) );

  replicate(cube,param.rep,nameOut);

  printBar(stdout);
  printf("  File %s was generated\n",nameOut);

  return 0;
}



int evalVoidVol(dataCube cube, dataRun param, char *name){

  unsigned int i,j,k;
  unsigned int n1,n2;
  unsigned int ncx,ncy,ncz,nct;
  unsigned int cbs;
  double vecA[3],vecB[3],vecC[3];
  double tmp[3],vol0,volt;
  double den,volvac;


  printBanner("Void volume",stdout);

  ncx = cube.pts[0] - 1;
  ncy = cube.pts[1] - 1;
  ncz = cube.pts[2] - 1;

  nct = ncx*ncy*ncz;
  
  vecA[0] = cube.mvec[0]; vecA[1] = cube.mvec[1]; vecA[2] = cube.mvec[2];
  vecB[0] = cube.mvec[3]; vecB[1] = cube.mvec[4]; vecB[2] = cube.mvec[5];
  vecC[0] = cube.mvec[6]; vecC[1] = cube.mvec[7]; vecC[2] = cube.mvec[8];

  crossProduct(vecA,vecB,tmp);
  vol0 = fabs( dotProduct(vecC,tmp) );

  volt = vol0*nct;

  n1 = cube.pts[1]*cube.pts[2];
  n2 = cube.pts[2];
  
  cbs=0;
  for( i = 0 ; i < ncx ; i++ ){
    for( j = 0 ; j < ncy ; j++ ){
      for( k = 0 ; k < ncz ; k++ ){
         den = getDenInCube(i,j,k,n1,n2,cube.field); 
         if( den <= param.vac )
           cbs++;
      }
    }
  }

  volvac = cbs*vol0;

  
   
  printf("   Cutoff in the density  : % 12.6lf a.u.\n",param.vac);
  printf("   Volume of the element  : % 12.6lf a.u.\n",vol0);
  printf("   Volume of the cell     : % 12.6lf a.u.\n",volt);
  printf("   Volume of the void     : % 12.6lf a.u.\n",volvac);
  printf("   Percentaje of the void : % 12.3lf \%\n",(volvac*100.0)/volt);

  char nameout[128];
  double x,y,z;
  FILE *out;

  sprintf(nameout,"%sVoid.xyz",name);

  openFile(&out,nameout,"w+");

  fprintf(out," %10d\n",cbs);
  fprintf(out," VOIDS for Cube3D project\n");

  for( i = 0 ; i < ncx ; i++ ){
    for( j = 0 ; j < ncy ; j++ ){
      for( k = 0 ; k < ncz ; k++ ){
         den = getDenInCube(i,j,k,n1,n2,cube.field); 
         if( den <= param.vac ){
           x = cube.min[0]  + ((double) i + 0.5 )*cube.hvec[0];
           y = cube.min[1]  + ((double) j + 0.5 )*cube.hvec[1];
           z = cube.min[2]  + ((double) k + 0.5 )*cube.hvec[2];

           x *= B2A;
           y *= B2A;
           z *= B2A;
           
           fprintf(out," Voi % 10.6lf % 10.6lf % 10.6lf\n",x,y,z);
         }
      }
    }
  }


  printBar(stdout);

  printf("  File %s was generated\n",nameout);


  fclose(out);
  
  return 0;
}

double getDenInCube(int i, int j, int k, int n1, int n2, double *field){
  
  int ip = i+1;
  int jp = j+1;
  int kp = k+1;

  double den;

  den  = field[i *n1 + j *n2 + k ];
  den += field[i *n1 + j *n2 + kp];
  den += field[i *n1 + jp*n2 + k ];
  den += field[i *n1 + jp*n2 + kp];
  
  den += field[ip*n1 + j * n2 + k ];
  den += field[ip*n1 + j * n2 + kp];
  den += field[ip*n1 + jp* n2 + k ];
  den += field[ip*n1 + jp* n2 + kp];

  return ( 0.125*den);
}

void getLogField (dataCube cube, double min2){

  unsigned int i;

  for( i = 0; i < cube.npt ; i++ )
    cube.field[i] = log(cube.field[i] - min2);


}