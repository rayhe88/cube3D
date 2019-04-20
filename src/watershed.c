#include "watershed.h"
#include "file.h"

#include "fields.h"

#include "merge.h"
#include "timing.h"

#define FRAC 0.01

void sortField(dataCube cube, dataRun param, int *fieldIdx, double *fieldSort){

  int i,npt;

  timeType ti,tf;
  double secs;

  npt = cube.npt;
  // Copy field from cube to fsort
  for(i=0; i < npt; i++){
    fieldIdx[i] = i;
    fieldSort[i] = cube.field[i];
  }

  gettimeofday(&ti,NULL);
  mergeSort(0,npt-1,fieldSort,fieldIdx);
  gettimeofday(&tf,NULL);

  secs = timing(&tf,&ti);

  printf(" Merge sort time : % 20.16g milliseconds\n",secs);

  // printing the field and idx
  FILE *out;
  
  openFile(&out,"sortField.dat","w+");

  for(i=0;i<npt;i++)
    fprintf(out," %12d  %12d  % 12.8E \n",i,
      fieldIdx[i],fieldSort[i]);

  fclose(out);
  
  double x,y,z;
  int j,k,idx = fieldIdx[0];
  int npy = cube.pts[1];
  int npz = cube.pts[2];

  i = idx / (npy * npz);
  j = (idx - i * (npy * npz)) / npz;
  k = idx % npz;

  x = cube.min[0] + i*cube.hvec[0];
  y = cube.min[1] + j*cube.hvec[1];
  z = cube.min[2] + k*cube.hvec[2];

  printf(" Indices : %5d %5d %5d\n",i,j,k);

  double xO1,yO1,zO1;
  double xH1,yH1,zH1;
  double xH2,yH2,zH2;
  xO1 = cube.coor[0];  yO1 = cube.coor[1];  zO1 = cube.coor[2];
  xH1 = cube.coor[3];  yH1 = cube.coor[4];  zH1 = cube.coor[5];
  xH2 = cube.coor[6];  yH2 = cube.coor[7];  zH2 = cube.coor[8];

  double r1 = sqrt( pow(x-xO1,2.) + pow(y-yO1,2.) + pow(z-zO1,2.));
  double r2 = sqrt( pow(x-xH1,2.) + pow(y-yH1,2.) + pow(z-zH1,2.));
  double r3 = sqrt( pow(x-xH2,2.) + pow(y-yH2,2.) + pow(z-zH2,2.));

  

  printf(" Comp = %5d %5d\n",idx,i*npy*npz + j*npz + k);  


  printf(" Coordinates of %12d : % 12.6lf % 12.6lf % 12.6lf\n",
          fieldIdx[0],x,y,z);

  printf(" Distancias \n");
  printf(" r - O1 : % 12.8lf\n",r1);
  printf(" r - H1 : % 12.8lf\n",r2);
  printf(" r - H2 : % 12.8lf\n",r3);

  printf(" h3 : % 12.8lf\n",sqrt(cube.hvec[0]*cube.hvec[0]+
                                 cube.hvec[1]*cube.hvec[1]+
                                 cube.hvec[2]*cube.hvec[2]));
    
}

int waterShed( dataCube cube, dataRun p, const double *matU,int *center, int *idxField,double *field){
  int t=0;
  int ny,nz;
  int n1,n2;
  int idxglb,iter;
  int index[3];
  int vecn[3];
  double grad[3];

  int (*fpol)(double,double,double,double*,double*);

  ny = cube.pts[1];
  nz = cube.pts[2];

  n1 = ny*nz;
  n2 = nz;

  switch( p.pol ){
    case  1: fpol = &gradPol01; break;
    case  2: fpol = &gradPol02; break;
    case  3: fpol = &gradPol03; break;
    case  4: fpol = &gradPol04; break;
    case  5: fpol = &gradPol05; break;
    case  6: fpol = &gradPol06; break;
    default: fpol = &gradPol02;

  }

  int npt = cube.npt;

  for (iter = 0; iter< npt; iter++)
    center[iter] = -1;

  center[idxField[0]] = 1;

  //printf(" %4d ",center[0]);

  int k = 1;

  for( iter=1; iter < 20; iter++ ){
    idxglb = idxField[iter];

    index[0] = idxglb/n1;
    index[1] = (idxglb - index[0]*n1)/n2;
    index[2] = idxglb%n2;    

    ///
    //
//  i = idx / (npy * npz);
//  j = (idx - i * (npy * npz)) / npz;
//  k = idx % npz;
    //
    //


    gradientVec(index,cube.pts,p,cube.hvec,cube.field,matU,fpol,grad);


    center[idxField[iter]] = asignaCenters( iter,idxField,center,index,cube.pts,cube.hvec,cube.min, grad );

    printf("Centro asignado: %4d \n",center[idxField[iter]]);
    /*k++;
    if( k == 10 || iter == npt - 1 ){
      printf("\n");
      k=0;
    }
    */
  }

  

  return 0;
}


int asignaCenters( int iter,int *idxField, int *center, int index[3], int npts[3], double hvec[3], 
                   double min[3], double grad[3] ){
    
    int idx[3];
    double r[3];
    static int atom=1;
    r[0] = min[0] + index[0]*hvec[0];
    r[1] = min[1] + index[1]*hvec[1];
    r[2] = min[2] + index[2]*hvec[2];

    grad[0] *= ( 1. + FRAC );
    grad[1] *= ( 1. + FRAC );
    grad[2] *= ( 1. + FRAC );

    r[0] += grad[0];
    r[1] += grad[1];
    r[2] += grad[2];

    getIndex3D(npts,r,min,hvec,idx);
    idx[0] += 1;
    idx[1] += 1;
    idx[2] += 1;

    int newIdx = idx[0]*npts[1]*npts[2] + idx[1]*npts[2] + idx[2];


    printf(" Entra     %5d %5d %5d :%8d ", index[0],index[1],index[2], index[0]*npts[1]*npts[2] + index[1]*npts[2] + index[2]);
    printf(" Sales %5d %5d %5d :%8d ", idx[0],idx[1],idx[2], newIdx);;

    printf(" Center : %4d ",center[newIdx]);

    if ( center[newIdx] > 0 )
      return center[newIdx];
    else
      return ++atom;

}

