#include "file.h"
#include "fields.h"
#include "geomData.h"
#include "lagrange2.h"
#include "cubeIndex.h"
#include "utils.h"
#include "numBondPath.h"
#include "transU.h"

int getFieldInLine(double min0,dataCube cube, dataRun param, const double *matU,char *name){

  printBanner("Field-line",stdout);

  int at1,at2;
  int n;
  double dist;
  double h;
  double r1[3],r2[3];
  double ux,uy,uz,norm;

  double f[param.size];
  double xx[param.pol + 1];
  double yy[param.pol + 1];
  double zz[param.pol + 1];

  FILE *out;
  char nameOut[128];

  at1 = param.geom[0];
  at2 = param.geom[1];
  if( at1 > cube.natm-1 || at2 > cube.natm-1)
    return 1;

  sprintf(nameOut,"%s_line%d-%d.dat",name,at1,at2);
  openFile(&out,nameOut,"w+");

  at1--;
  at2--;

  // the coordinates are assigned to r1 and r2 vectors
  //

  r1[0] = cube.coor[3*at1];
  r1[1] = cube.coor[3*at1+1];
  r1[2] = cube.coor[3*at1+2];

  r2[0] = cube.coor[3*at2];
  r2[1] = cube.coor[3*at2+1];
  r2[2] = cube.coor[3*at2+2];

  fprintf(out,"# Coordinates atom %4d:  % 10.6lf % 10.6lf % 10.6lf\n",at1+1,r1[0],r1[1],r1[2]);
  fprintf(out,"# Coordinates atom %4d:  % 10.6lf % 10.6lf % 10.6lf\n",at2+1,r2[0],r2[1],r2[2]);
  dist = distance(r1,r2);

  ux = r2[0] - r1[0];
  uy = r2[1] - r1[1];
  uz = r2[2] - r1[2];

  norm = getNorm(ux,uy,uz);
  
  ux /= norm;
  uy /= norm;
  uz /= norm;

  n = ceil(dist*(NPUA1));
  h = dist/(double) (n-1);

  double (*gfun) (double*);
  int (*gNum) (double,double,double,double*,double*,double*,double*,int,int,const double*,double,double*);

  switch( param.task ){
    case RED: gfun = &getRed; 
              gNum = &gradient3DLog;
              break;
    case GRA: gfun = &getGrd; 
              gNum = &gradient3DLog;
              break;
    case LAP: gfun = &getLap;
              gNum = &getDerivatives3DLog;
              break;
    case KIN: gfun = &getKin;
              gNum = &getDerivatives3DLog;
              break;
    case VIR: gfun = &getVir;
              gNum = &getDerivatives3DLog;
              break;
    default : gfun = &getGrd;
              gNum = &gradient3DLog;
              break;
  }


  printf("  Line between the  atoms %3d and %3d\n",at1+1,at2+1);
  printf("  Distance in Angstrom  :  % 10.6lf\n",dist*B2A);
  printf("  Number of points      :  % 10d\n",n);
  printf("  Step in Angstrom      :  % 10.6lf\n",h*B2A);


  fprintf(out,"# Distance between atoms %3d and %3d is : %10.6lf Angstrom\n",at1+1,at2+2,dist*B2A);
  fprintf(out,"# Unitary vector between atoms : % 10.6lf i + % 10.6lf j + % 10.6lf k\n",ux,uy,uz);
  fprintf(out,"# Number of points  % 4d  step : % 10.6lf\n",n,h);
  ux *= h;
  uy *= h;
  uz *= h;
  fprintf(out,"#    h    vector between atoms : % 10.6lf i + % 10.6lf j + % 10.6lf k\n#\n",ux,uy,uz);
  fprintf(out,"#        r         field        field0\n");

  int i;
  double r[3];
  int idx[3];
  double val[10];
  double f0;

  for(i=0;i<n;i++){
    r[0] = r1[0] + i*ux;
    r[1] = r1[1] + i*uy;
    r[2] = r1[2] + i*uz;

    getIndex3D(cube.pts,r,cube.min,cube.hvec,idx);
    loadLocalField(idx,cube,param,xx,yy,zz,f);

    gNum(r[0],r[1],r[2],xx,yy,zz,f,param.pol,param.orth,matU,min0,val);
    f0 = gfun(val);
    
    fprintf(out," % 12.6lf % 12.6lf % 12.6lf\n",distance(r1,r)*B2A,f0,val[0]);

  }

  fclose(out);

  printBar(stdout);
  printf("  File %s was generated\n",nameOut);

  return 0;
}

int getFieldInPlane(double min0,dataCube cube, dataRun param, const double *matU,char *name){
  
  printBanner("Field-plane",stdout);

  int at1,at2,at3;
  int n;
  double dist,h;
  double r1[3],r2[3],r3[3];
  double r[3],q[3];

  double matT[9];

  getMatInv(matU,matT);

  double f[param.size];
  double xx[param.pol + 1];
  double yy[param.pol + 1];
  double zz[param.pol + 1];

  FILE *out;
  char nameOut[128];


  at1 = param.geom[0];
  at2 = param.geom[1];
  at3 = param.geom[2];

  if( at1 > cube.natm || at2 > cube.natm || at3 > cube.natm )
    return 1;
   

  sprintf(nameOut,"%s_plane%d-%d-%d.dat",name,at1,at2,at3);
  openFile(&out,nameOut,"w+");

  printf("  Plane conformed by the atoms %3d, %3d and %3d\n",at1,at2,at3);
  reOrderAtoms(&at1,&at2,&at3,cube);

  at1--;
  at2--;
  at3--;

  r1[0] = cube.coor[3*at1];
  r1[1] = cube.coor[3*at1+1];
  r1[2] = cube.coor[3*at1+2];

  r2[0] = cube.coor[3*at2];
  r2[1] = cube.coor[3*at2+1];
  r2[2] = cube.coor[3*at2+2];

  r3[0] = cube.coor[3*at3];
  r3[1] = cube.coor[3*at3+1];
  r3[2] = cube.coor[3*at3+2];

  double r0[3],vL[3],vM[3],vN[3],cg[3];

  dist = getVecLMN(r1,r2,r3,r0,cg,vL,vM,vN);

  //dist = distance(r0,cg);
  n = ceil(dist*(NPUA1));
  h = dist/(double) (n-1);

  printf("  Number of points      :  % 10d x %-10d\n",n,n);
  printf("  Distance in Angstrom  :        % 10.6lf\n",n*h*B2A);
  printf("  Step in Angstrom      :        % 10.6lf\n",h*B2A);

  double (*gfun) (double*);
  int (*gNum) (double,double,double,double*,double*,double*,double*,int,int,const double*,double,double*);

  switch( param.task ){
    case RED: gfun = &getRed; 
              gNum = &gradient3DLog;
              break;
    case GRA: gfun = &getGrd; 
              gNum = &gradient3DLog;
              break;
    case LAP: gfun = &getLap;
              gNum = &getDerivatives3DLog;
              break;
    case KIN: gfun = &getLap;
              gNum = &getDerivatives3DLog;
              break;
    case VIR: gfun = &getLap;
              gNum = &getDerivatives3DLog;
              break;
    default : gfun = &getGrd;
              gNum = &gradient3DLog;
              break;
  }

  int i,j;
  int idx[3];
  double val[10];
  double f0;
  double ejeL,ejeM;
  double cgL,cgM;

  fprintf(out,"# eje L    eje M   field    field0\n");
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      r[0] = r0[0] + i*h*vL[0] + j*h*vM[0];
      r[1] = r0[1] + i*h*vL[1] + j*h*vM[1];
      r[2] = r0[2] + i*h*vL[2] + j*h*vM[2];

      getRiU(r,matT,q);

      ejeL = dotProduct(r,vL);
      ejeM = dotProduct(r,vM);
      cgL  = dotProduct(cg,vL);
      cgM  = dotProduct(cg,vM);

      ejeL -= cgL;
      ejeM -= cgM;


      getIndex3D(cube.pts,q,cube.min,cube.hvec,idx);
      loadLocalField(idx,cube,param,xx,yy,zz,f);

      gNum(q[0],q[1],q[2],xx,yy,zz,f,param.pol,param.orth,matU,min0,val);
      f0 = gfun(val);

      fprintf(out," % 12.6lf % 12.6lf  % 14.8lf % 14.8lf \n",ejeL*B2A,ejeM*B2A,f0,val[0]);

    }
      fprintf(out," \n");
  }

  fclose(out);

  printBar(stdout);
  printf("  File %s was generated\n",nameOut);

  return 0;
}

double getVecLMN( double r1[3], double r2[3], double r3[3], double r0[3], double cg[3],
                  double vecL[3], double vecM[3], double vecN[3]){
  
  double tmp1[3],tmp2[3];
  double vecACm[3];
  double dist,norm;

  cg[0] = (r1[0] + r2[0] + r3[0])/(double) 3.;
  cg[1] = (r1[1] + r2[1] + r3[1])/(double) 3.;
  cg[2] = (r1[2] + r2[2] + r3[2])/(double) 3.;

  dist = 1.75*distance(r1,cg);

  tmp1[0] = (r1[0] - r2[0]);
  tmp1[1] = (r1[1] - r2[1]);
  tmp1[2] = (r1[2] - r2[2]);

  tmp2[0] = (r3[0] - r2[0]);
  tmp2[1] = (r3[1] - r2[1]);
  tmp2[2] = (r3[2] - r2[2]);

  vecACm[0] = 0.5*(r1[0] - r3[0]);
  vecACm[1] = 0.5*(r1[1] - r3[1]);
  vecACm[2] = 0.5*(r1[2] - r3[2]);

  crossProduct(tmp1,tmp2,vecN);
  norm = getNormVec(vecN);
  vecN[0] /= norm;
  vecN[1] /= norm;
  vecN[2] /= norm;

  crossProduct(vecN,vecACm,vecM);
  norm = getNormVec(vecM);
  vecM[0] /=(-1.*norm);
  vecM[1] /=(-1.*norm);
  vecM[2] /=(-1.*norm);

  crossProduct(vecN,vecM,vecL);
  norm = getNormVec(vecL);
  vecL[0] /= norm;
  vecL[1] /= norm;
  vecL[2] /= norm;
  
  r0[0] = cg[0] - dist*vecL[0] - dist*vecM[0];
  r0[1] = cg[1] - dist*vecL[1] - dist*vecM[1];
  r0[2] = cg[2] - dist*vecL[2] - dist*vecM[2];

  double p1[3];

  p1[0] = cg[0] - dist*vecL[0] + dist*vecM[0];
  p1[1] = cg[1] - dist*vecL[1] + dist*vecM[1];
  p1[2] = cg[2] - dist*vecL[2] + dist*vecM[2];

  double m[3],l[3],minmax[2];

  l[0] = dotProduct(r1,vecL);
  l[1] = dotProduct(r2,vecL);
  l[2] = dotProduct(r3,vecL);

  m[0] = dotProduct(r1,vecM);
  m[1] = dotProduct(r2,vecM);
  m[2] = dotProduct(r3,vecM);

  sortCoor(l,m,minmax);

  return distance(r0,p1);
}

void sortCoor(double l[3],double m[3],double minmax[2]){

  int i,j,tmp;
  double coorm[3],coorl[3];

  double extra = 1.5;

  for(i=0;i<3;i++){
    coorm[i] = m[i];
    coorl[i] = l[i];
  }

  for(i=0;i<3;i++)
    for(j=i+1;j<3;j++){
      if( coorm[i] > coorm[j]){
        tmp      = coorm[i];
        coorm[i] = coorm[j];
        coorm[j] = tmp;
      }   
    }   

  for(i=0;i<3;i++)
    for(j=i+1;j<3;j++){
      if( coorl[i] > coorl[j]){
        tmp      = coorl[i];
        coorl[i] = coorl[j];
        coorl[j] = tmp;
      }   
    }   

  coorm[0] -= extra;
  coorl[0] -= extra;

  coorm[2] += extra;
  coorl[2] += extra;

  if( coorm[0] < coorl[0] )
    minmax[0] = coorm[0];
  else
    minmax[0] = coorl[0];

  if( coorm[2] > coorl[2] )
    minmax[1] = coorm[2];
  else
    minmax[1] = coorl[2];


}

void origin( double vecL[3], double vecM[3], double vecN[3],double minmax[2], double r0[3]){
  double m0,n0,l0;
  double rx,ry,rz;
  double mx,my,mz;
  double nx,ny,nz;
  double lx,ly,lz;
  double det;

  l0 = minmax[0];
  m0 = minmax[0];
  n0 = 0.;

  lx = vecL[0]; ly = vecL[1]; lz = vecL[2];
  mx = vecM[0]; my = vecM[1]; mz = vecM[2];
  nx = vecN[0]; ny = vecN[1]; nz = vecN[2];

  det  = mx*( ny*lz - nz*ly );
  det += my*( nz*lx - nx*lz );
  det += mz*( nx*ly - ny*lx );

  rx = -lz*my*n0 + ly*mz*n0 + lz*m0*ny - l0*mz*ny - ly*m0*nz + l0*my*nz;
  ry =  lz*mx*n0 - lx*mz*n0 - lz*m0*nx + l0*mz*nx + lx*m0*nz - l0*mx*nz;
  rz = -ly*mx*n0 + lx*my*n0 + ly*m0*nx - l0*my*nx - lx*m0*ny + l0*mx*ny;

  r0[0] = rx/det;
  r0[1] = ry/det;
  r0[2] = rz/det;
 
}
void  reOrderAtoms(int *at1,int *at2, int *at3,dataCube cube){
  int i;
  int atA,atB,atC;

  double coorA[3],coorB[3],coorC[3];
  double vecAB[3],vecAC[3],vecBC[3];
  double nAB,nAC,nBC;

  atA = (*at1); atB = (*at2); atC = (*at3);

  atA--;
  atB--;
  atC--;

  for( i = 0; i < 3; i++){
    coorA[i] = cube.coor[3*atA+i];
    coorB[i] = cube.coor[3*atB+i];
    coorC[i] = cube.coor[3*atC+i];

    vecAB[i] = coorA[i] - coorB[i];
    vecAC[i] = coorA[i] - coorC[i];
    vecBC[i] = coorB[i] - coorC[i];
  }

  nAB = getNormVec(vecAB);
  nAC = getNormVec(vecAC);
  nBC = getNormVec(vecBC);


  if( nAB > nAC && nAB > nBC ){
    (*at2) = atC+1;
    (*at1) = atA+1;
    (*at3) = atB+1;
  }else{
    if( nAC > nAB && nAC > nBC){
      (*at2) = atB+1;
      (*at1) = atA+1;
      (*at3) = atC+1;
    }else{
      if( nBC > nAB && nBC > nAC){
        (*at2) = atA+1;
        (*at1) = atB+1;
        (*at3) = atC+1;
      }else{
        (*at1) = atA+1;
        (*at2) = atB+1;
        (*at3) = atC+1;

      }
    }
  }

}
