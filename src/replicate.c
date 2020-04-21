#include "struct.h"
#include "replicate.h"
#include "array.h"
#include "lectura.h"
#include "file.h"
#include "utils.h"

int repInt(double rep){
    if ( rep >= 1. )
        return floor(rep) - 1;
    else
        return 0;
}

int repInt2(double rep){
    if ( rep >= 1. )
        return (int)ceil(rep);
    else
        return 1;
}
int replicate(int *npts,dataCube cubeOld, double rep[3],const double matT[9],
              const double matU[9],const char *name){

  char text[64];
  int i, npt2,check;
  int tamanio;
  int newxp = (int) ceil(cubeOld.pts[0]*rep[0]);
  int newyp = (int) ceil(cubeOld.pts[1]*rep[1]);
  int newzp = (int) ceil(cubeOld.pts[2]*rep[2]);
  int *zatm2;
  double *coor2;
  double *fieldNew;

  FILE *out;
  dataCube cubeNew;

  sprintf(text," Cube replicate %6lf x %6lf x %6lf times",rep[0],rep[1],rep[2]);
 

  openFile(&out,name,"w+");

  checkBoundaryCond(cubeOld,&check);

  printf(" Check value for boundary conditions: %d\n",check);

  newxp -= ( repInt( rep[0] )) * check;
  newyp -= ( repInt( rep[1] )) * check;
  newzp -= ( repInt( rep[2] )) * check;

  npt2 = newxp*newyp*newzp;

  createArrayDou ( npt2,&fieldNew,"Field2");

  cubeNew.pts[0] = newxp;
  cubeNew.pts[1] = newyp;
  cubeNew.pts[2] = newzp;
  cubeNew.npt    = npt2;

  cubeNew.natm = cubeOld.natm;
  for(i=0;i<3;i++){
    cubeNew.min[i]  = cubeOld.min[i];
    cubeNew.hvec[i] = cubeOld.hvec[i];
  }
  for(i=0;i<9;i++)
    cubeNew.mvec[i] = cubeOld.mvec[i];

  for(i=0;i<npt2;i++){
    fieldNew[i] = cubeOld.field[getIndexOld(check,i,cubeOld.pts,cubeNew.pts)];
  }

  tamanio  = cubeOld.natm;
  tamanio *= repInt2(rep[0]);
  tamanio *= repInt2(rep[1]);
  tamanio *= repInt2(rep[2]);

  createArrayInt(tamanio,&zatm2,"New zatm");
  createArrayDou(3*tamanio,&coor2,"New coor");

  tamanio = replicateCoor(check,cubeOld,zatm2,coor2,rep, matT,matU);

  cubeNew.natm = tamanio;

  cubeNew.zatm  = zatm2;
  cubeNew.coor  = coor2;
  cubeNew.field = fieldNew;

  (*npts) = tamanio;

  printCube(text,cubeNew,out);

  unloadData(&cubeNew,&zatm2,&coor2,&fieldNew);

  fclose(out);

  return 0;
}


int getIndexOld(int check, int idx, int oldpts[3], int newpts[3] ){

  int oldx,oldy,oldz;
  int newx,newy,newz;
  int oldnpx,oldnpy,oldnpz;
  int newnpy,newnpz;

  int oldIdx;

  oldnpx = oldpts[0];
  oldnpy = oldpts[1];  
  oldnpz = oldpts[2]; 

  newnpy = newpts[1];
  newnpz = newpts[2];


  newx = (idx/newnpz)/newnpy;
  newy = (idx/newnpz)%newnpy;
  newz = (idx%newnpz);

  // NOTA USAR -1 para cubos de Crystal
  //            0 para cubos de GPAW
  // esto viene dado por check
 
  oldx = newx%(oldnpx-check);
  oldy = newy%(oldnpy-check);
  oldz = newz%(oldnpz-check);

  oldIdx = oldx*oldnpz*oldnpy + oldy*oldnpz + oldz;
  

  return oldIdx;
}

void getRiT( double r[3], const double u[9], double q[3]){
    double rx,ry,rz;
    rx = r[0]; ry = r[1]; rz = r[2];

    q[0] =  rx - ry*(u[1]/u[4]);                                                  
    q[0] += rz*((u[1]*u[5] - u[2]*u[4])/(u[4]*u[8]));                             
    q[1]  = ry/u[4] - rz*u[5]/(u[4]*u[8]);                                        
    q[2]  = rz/u[8];  

}

int replicateCoor(int check, dataCube cubeOld, int *zatm2, double *coor2, 
                  double *rep,const double* matT,const double* matU){

  int i,j,k,l,n;
  int nrx,nry,nrz;

  int natm = cubeOld.natm;
  int npx = cubeOld.pts[0] - check;
  int npy = cubeOld.pts[1] - check;
  int npz = cubeOld.pts[2] - check;

  double hrx[3],hry[3],hrz[3];
  double hqx[3],hqy[3],hqz[3];
  double deltax,deltay,deltaz;

  hrx[0] = npx*(cubeOld.mvec[0]);
  hrx[1] = npx*(cubeOld.mvec[1]);
  hrx[2] = npx*(cubeOld.mvec[2]);
  getRiT(hrx, matT, hqx);

  hry[0] = npy*(cubeOld.mvec[3]);
  hry[1] = npy*(cubeOld.mvec[4]);
  hry[2] = npy*(cubeOld.mvec[5]);
  getRiT(hry, matT, hqy);

  hrz[0] = npz*(cubeOld.mvec[6]);
  hrz[1] = npz*(cubeOld.mvec[7]);
  hrz[2] = npz*(cubeOld.mvec[8]);
  getRiT(hrz, matT, hqz);

  deltax = getNormVec(hqx);
  deltay = getNormVec(hqy);
  deltaz = getNormVec(hqz);

  double max_x,max_y,max_z;
  max_x = deltax * rep[0];
  max_y = deltay * rep[1];
  max_z = deltaz * rep[2];

  nrx = (int) ceil ( rep[0] );
  nry = (int) ceil ( rep[1] );
  nrz = (int) ceil ( rep[2] );

  double ra[3],qa[3];
 
  l=0;
  for(i=0;i<nrx;i++){
    for(j=0;j<nry;j++){
      for(k=0;k<nrz;k++){
        for(n=0;n<natm;n++){
            ra[0] = cubeOld.coor[3*n];
            ra[1] = cubeOld.coor[3*n+1];
            ra[2] = cubeOld.coor[3*n+2];

            getRiT(ra, matT, qa);
            qa[0] += i*deltax;
            qa[1] += j*deltay;
            qa[2] += k*deltaz;

            if( qa[0] <= max_x && qa[1] <= max_y && qa[2] <= max_z){
                getRiT(qa, matU, ra);

                zatm2[l] = cubeOld.zatm[n];
                coor2[3*l]   = ra[0];
                coor2[3*l+1] = ra[1];
                coor2[3*l+2] = ra[2];
                l++;
            }
        }
      }
    }
  }
 
  return l-1;
}
