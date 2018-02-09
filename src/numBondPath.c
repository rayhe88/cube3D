#include <stdio.h>
#include <stdlib.h>

#include "file.h"
#include "array.h"
#include "analys.h"
#include "matvec.h"
#include "critical.h"
#include "mathTools.h"
#include "pruebaCoef.h"
#include "numBondPath.h"
#include "lagrangePol.h"
#include "tableP.h"

#define TOLERANCE 0.01
#define MAXPTS 800

int isInsideCube(double vecIn[3],double minmax[6]){


  if( vecIn[0] < minmax[0] ) return 1;
  if( vecIn[0] > minmax[1] ) return 1;
  if( vecIn[1] < minmax[2] ) return 1;
  if( vecIn[1] > minmax[3] ) return 1;
  if( vecIn[2] < minmax[4] ) return 1;
  if( vecIn[2] > minmax[5] ) return 1;

  return 0;
}

int myIsNanInf(double val){

  int ret=0;

  ret += isnan(val);
  ret += isinf(val);

  return ret;
}

int myIsNanInf_V3(double vec[3]){

  int ret = 0;
  
  ret += myIsNanInf(vec[0]);
  ret += myIsNanInf(vec[1]);
  ret += myIsNanInf(vec[2]);

  return ret;

}

int bondPath(int poly,int nna, int nato, int nCrit,int *type,double *coor, double *coorCrit,int *pts, 
             double *min,double *hvec, double *field,double *matT,const char* name){
  int i,j,k,id[3],itmp,iterap,iteran;
  int attractors, amico;
  int nucleo1,nucleo2;
  char tmpname[120];
  double *coorAttr;
  double *coef,val[10];
  double matHH[3][3];
  double eigenVec[3][3];
  double eigenVal[3],norm,dist,difmin,rij;
  double v1,v2,v3;
  double x,y,z;
  double xn,yn,zn;
  double x0,y0,z0;
  double k1x,k1y,k1z;
  double k2x,k2y,k2z;
  double ex0,ey0,ez0;
  double xatm,yatm,zatm;
  double xcrit,ycrit,zcrit;
  double fac = 0.52917;
  double vecIn[3],vecOut[3];
  double minmax[6];
  FILE *tmp;

  int size = pow(poly+1,3);


  minmax[0] = min[0]; minmax[1] = min[0] + pts[0]*hvec[0];
  minmax[2] = min[1]; minmax[3] = min[1] + pts[1]*hvec[1];
  minmax[4] = min[2]; minmax[5] = min[2] + pts[2]*hvec[2];

  createArrayDou(size,&coef,"Coef");


  tmpFile(&tmp,".bpath_",tmpname,"w+");

  attractors = nna + nato;

  createArrayDou(3*attractors,&coorAttr,"BP01");

  for(i=0;i<nato;i++){
    coorAttr[3*i]   = coor[3*i];
    coorAttr[3*i+1] = coor[3*i+1];
    coorAttr[3*i+2] = coor[3*i+2];
  }
  for(j=0;j<nna;j++){
    coorAttr[3*i]   = coorCrit[3*j];
    coorAttr[3*i+1] = coorCrit[3*j+1];
    coorAttr[3*i+2] = coorCrit[3*j+2];
    i++;
  }
  int bp=0;
  for(i=0; i < nCrit; i++){
    if( type[i] == -3 ){
      bp++;
      nucleo1 = nucleo2 = 0;


      xcrit = coorCrit[3*i];
      ycrit = coorCrit[3*i+1];
      zcrit = coorCrit[3*i+2];

      id[0] = 0; id[1] = 0; id[2] = 0;
      getData(poly,pts,id,xcrit,ycrit,zcrit,min,hvec,field,coef);
      NumericalCrit02(xcrit,ycrit,zcrit,coef,val);
      //evalPol2(poly,xcrit,ycrit,zcrit,coef,val);

      matHH[0][0] = val[4]; matHH[0][1] = val[5]; matHH[0][2] = val[6];
      matHH[1][0] = val[5]; matHH[1][1] = val[7]; matHH[1][2] = val[8];
      matHH[2][0] = val[6]; matHH[2][1] = val[8]; matHH[2][2] = val[9];

      jacobi(matHH,eigenVal,eigenVec);

      v1 = eigenVec[0][0];  v2 = eigenVec[0][1];  v3 = eigenVec[0][2]; 

      norm = getNorm(v1,v2,v3);

      xn = xcrit + (v1*0.02)/norm;
      yn = ycrit + (v2*0.02)/norm;
      zn = zcrit + (v3*0.02)/norm;

      iterap = 0; dist = 4.;

      id[0] = 0; id[1] = 0; id[2] = 0;

      while( iterap < MAXPTS && dist > 0 ){
        x0 = xn;
        y0 = yn;
        z0 = zn;

        getData(poly,pts,id,x0,y0,z0,min,hvec,field,coef);
        //evalPol2(poly,x0,y0,z0,coef,val);
        NumericalCrit02(x0,y0,z0,coef,val);

        ex0 = val[1]; 
        ey0 = val[2]; 
        ez0 = val[3];

        norm = getNorm(ex0,ey0,ez0);

        k1x = ex0/norm;
        k1y = ey0/norm;
        k1z = ez0/norm;

        x  = x0 + 0.01*k1x;
        y  = y0 + 0.01*k1y;
        z  = z0 + 0.01*k1z;

        getData(poly,pts,id,x,y,z,min,hvec,field,coef);
        //evalPol2(poly,x,y,z,coef,val);
        NumericalCrit02(x,y,z,coef,val);

        ex0 = val[1]; 
        ey0 = val[2]; 
        ez0 = val[3];

        norm = getNorm(ex0,ey0,ez0);

        k2x = ex0/norm;
        k2y = ey0/norm;
        k2z = ez0/norm;

        xn = x0 + 0.02*k2x;
        yn = y0 + 0.02*k2y;
        zn = z0 + 0.02*k2z;

        amico = 0;
        difmin = 1.E4;


        for(k=0; k < attractors; k++){
          xatm = coorAttr[3*k];
          yatm = coorAttr[3*k+1];
          zatm = coorAttr[3*k+2];

          rij = distance(xn,xatm,yn,yatm,zn,zatm);

          if( rij < difmin){
            difmin = rij;
            nucleo1 = k;
          }
          if( rij <= TOLERANCE){
            amico = 1;
            nucleo1 = k;
          }

        }// end for k

        if( amico == 0 ){
          iterap++;
          dist = 4.;
          vecIn[0] = xn*fac; vecIn[1] = yn*fac; vecIn[2] = zn*fac;

          if( isInsideCube(vecIn,minmax) == 0 ){
            matVecProduct(vecIn,matT,vecOut);
          
            if( myIsNanInf_V3(vecOut) == 0 )
              fprintf(tmp," BP%04d  % 10.6lf   % 10.6lf  % 10.6lf\n",i+1-nna,vecOut[0],vecOut[1],vecOut[2]);
          }

        }else{
          dist = -1.;
        }
      }//end while


 // At this moment we change the direction
      norm = getNorm(v1,v2,v3);

      xn = xcrit - (v1*0.02)/norm;
      yn = ycrit - (v2*0.02)/norm;
      zn = zcrit - (v3*0.02)/norm;

      iteran = 0; dist = 4.;
   
      id[0] = 0; id[1] = 0; id[2] = 0;
      while( iteran < MAXPTS && dist > 0 ){
        x0 = xn;
        y0 = yn;
        z0 = zn;

        getData(poly,pts,id,x0,y0,z0,min,hvec,field,coef);

        NumericalCrit02(x0,y0,z0,coef,val);
        //evalPol2(poly,x0,y0,z0,coef,val);

        ex0 = val[1]; 
        ey0 = val[2]; 
        ez0 = val[3];

        norm = getNorm(ex0,ey0,ez0);

        k1x = ex0/norm;
        k1y = ey0/norm;
        k1z = ez0/norm;

        x  = x0 - 0.01*k1x;
        y  = y0 - 0.01*k1y;
        z  = z0 - 0.01*k1z;

        getData(poly,pts,id,x,y,z,min,hvec,field,coef);
        NumericalCrit02(x,y,z,coef,val);
        //evalPol2(poly,x,y,z,coef,val);

        ex0 = val[1]; 
        ey0 = val[2]; 
        ez0 = val[3];

        norm = getNorm(ex0,ey0,ez0);


        k2x = ex0/norm;
        k2y = ey0/norm;
        k2z = ez0/norm;

        xn = x0 + 0.02*k2x;
        yn = y0 + 0.02*k2y;
        zn = z0 + 0.02*k2z;

        amico = 0;
        difmin = 1.E4;

        for(k=0; k < attractors; k++){
          xatm = coorAttr[3*k];
          yatm = coorAttr[3*k+1];
          zatm = coorAttr[3*k+2];

          rij = distance(xn,xatm,yn,yatm,zn,zatm);

          if( rij < difmin){
            difmin = rij;
            nucleo2 = k;
          }
          if( rij <= TOLERANCE){
            amico = 1;
            nucleo2 = k;
          }

        }// end for k

        if( amico == 0 ){
          iteran++;
          dist = 4.;
          vecIn[0] = xn*fac; vecIn[1] = yn*fac; vecIn[2] = zn*fac;

          if( isInsideCube(vecIn,minmax) == 0 ){
            matVecProduct(vecIn,matT,vecOut);
            if( myIsNanInf_V3(vecOut) == 0 )
              fprintf(tmp," BP%04d  % 10.6lf   % 10.6lf  % 10.6lf\n",i+1-nna,vecOut[0],vecOut[1],vecOut[2]);
          }
        }else{
          dist = -1;
        }
      }//end while



      if( nucleo1 > nucleo2){
        itmp = nucleo1;
        nucleo1 = nucleo2;
        nucleo2 = itmp;
      }

      printf(" --> Punto critico %3d  Une : %2d %2d\n",bp,nucleo1+1,nucleo2+1);

    }


  }//fin for


  free(coef);

  rewind(tmp);

  int npoints=0;
  char c;

  while( (c = fgetc(tmp)) != EOF){
    if( c == '\n' )
      npoints++;
  }
  rewind(tmp);

  char namebpath[128];
  FILE *out;

  sprintf(namebpath,"%sBPath.xyz",name);

  openFile(&out,namebpath,"w+");

  fprintf(out," %10d\n",npoints);
  fprintf(out," Bond path file in Aangstrom\n");

  while( (c = fgetc(tmp)) != EOF ){
    fprintf(out,"%c",c);
  }


  fclose(tmp);
  fclose(out);

  remove(tmpname);
  printf(" Terminamos de imprimir los bondPaths\n");
  printf(" En el archivo [%s]\n",namebpath);
  return 0;
}

int SortCoordinates(int ncrit, int *type, double *coorcrit){
  int i,j,k;
  int *typeOrd;
  double *sort;

  createArrayInt( ncrit, &typeOrd, "tipo");
  createArrayDou( 3*ncrit, &sort, "coordenadas");

  for( i=0, j=0; i < ncrit; i++){
    if ( type[i] == -4 ) {// alias para NNACP
      typeOrd[j]  = type[i];
      sort[3*j]   = coorcrit[3*i];
      sort[3*j+1] = coorcrit[3*i+1];
      sort[3*j+2] = coorcrit[3*i+2];
      j++;
    }

  }

  for( i=0; i < ncrit; i++){
    if ( type[i] == -3 ) {// alias para BCP
      typeOrd[j]  = type[i];
      sort[3*j]   = coorcrit[3*i];
      sort[3*j+1] = coorcrit[3*i+1];
      sort[3*j+2] = coorcrit[3*i+2];
      j++;
    }
  }

  for( i=0; i < ncrit; i++){
    if ( type[i] == -2 ) {// alias para RCP
      typeOrd[j]  = type[i];
      sort[3*j]   = coorcrit[3*i];
      sort[3*j+1] = coorcrit[3*i+1];
      sort[3*j+2] = coorcrit[3*i+2];
      j++;
    }
  }

  for( i=0; i < ncrit; i++){
    if ( type[i] == -1 ) {// alias para CCP
      typeOrd[j]  = type[i];
      sort[3*j]   = coorcrit[3*i];
      sort[3*j+1] = coorcrit[3*i+1];
      sort[3*j+2] = coorcrit[3*i+2];
      j++;
    }
  }

  for( i=0; i < ncrit; i++){
    if ( type[i] == 0 ) {// alias para dummy
      typeOrd[j]  = type[i];
      sort[3*j]   = coorcrit[3*i];
      sort[3*j+1] = coorcrit[3*i+1];
      sort[3*j+2] = coorcrit[3*i+2];
      j++;
    }
  }
 
  k=0;
  for( i=0; i < ncrit; i++){
    type[i] = typeOrd[i];
    if( type[i] > -5 && type[i] < 0)
      k++;
    j = 3*i;
    coorcrit[j] = sort[j];
    coorcrit[j+1] = sort[j+1];
    coorcrit[j+2] = sort[j+2];

  }
 

 free(sort);
 free(typeOrd);
 return k;
}

int getData(int poly,int *pts, int *id,double x,double y, double z,double *min, double *hvec,double *field, double *coef){
  
  int indi,indj,indk;
  int npx,npy,npz,flag;
  double x1,y1,z1;
  double *fun;

  int size =  pow(poly+1,3);

  createArrayDou(size,&fun,"func");
  npx = pts[0];
  npy = pts[1];
  npz = pts[2];

  getCube(pts,x,y,z,min,hvec,&indi,&indj,&indk);
  flag  = fabs(indi-id[0]);
  flag += fabs(indj-id[1]);
  flag += fabs(indk-id[2]);


  if( flag != 0){
    x1 = min[0] + indi*hvec[0];
    y1 = min[1] + indj*hvec[1];
    z1 = min[2] + indk*hvec[2];
    loadField(poly,indi,indj,indk,npx,npy,npz,field,fun);

    coefficients(hvec[0],hvec[1],hvec[2],x1,y1,z1,fun,coef);
    //getCoeff(poly,hvec,x1,y1,z1,fun,coef);

    id[0] = indi; id[1] = indj; id[2] = indk;
  }


  free(fun);
  return 0;
}

  
int logFile(int poly,int nato, int nCrit,int *type,double *coor, double *coorCrit,int *pts, 
             double *min,double *hvec, double *field,double *matT,const char* name){

  int i,j;
  char logname[128];
  double factor = 0.52917;
  FILE *log;
  int ncc,nrc,nbc,nna;
  int ntyp[5];

  double *coef;
  double vecIn[3],vecOut[3];



  sprintf(logname,"%sCrit.log",name);
  openFile(&log,logname,"w+");
   
  ncc = nrc = nbc = nna = 0;

  for(i=0;i<nCrit;i++){
    if( type[i] == -1 ) ncc++;
    if( type[i] == -2 ) nrc++;
    if( type[i] == -3 ) nbc++;
    if( type[i] == -4 ) nna++;
  }
  ntyp[0] = ncc;
  ntyp[1] = nrc;
  ntyp[2] = nbc;
  ntyp[3] = nna;
  ntyp[4] = nato;


  int size = pow(poly+1,3);
  double x,y,z;
  int id[3];
  double val[10];
  char symb[4];


  createArrayDou(size,&coef,"Coef");

  printLog1(log,ntyp);


  id[0] = 0; id[1] = 0; id[2] = 0;
  for(i=0; i < nato; i++){
    x = coor[3*i];
    y = coor[3*i+1];
    z = coor[3*i+2];

    getData(poly,pts,id,x,y,z,min,hvec,field,coef);
    NumericalCrit02(x,y,z,coef,val);

    vecIn[0] = x*factor;
    vecIn[1] = y*factor;
    vecIn[2] = z*factor;
    matVecProduct(vecIn,matT,vecOut);

    getAtomicSymbol(zatm[i],4,symb);



    fprintf(log,"%5d  %5s",i+1,symb);

    fprintf(log," % 10.6lf % 10.6lf % 10.6lf % 14.8lf % 16.7lf\n",vecOut[0],vecOut[1],vecOut[2],val[0],val[4]+val[7]+val[9]);

  }


  if( nna > 0 ) lineDiv(log);

  for(i=0; i < nCrit; i++){

    if( type[i] == -4 ){
      x = coorCrit[3*i];
      y = coorCrit[3*i+1];
      z = coorCrit[3*i+2];

      id[0] = 0; id[1] = 0; id[2] = 0;
      getData(poly,pts,id,x,y,z,min,hvec,field,coef);
      NumericalCrit02(x,y,z,coef,val);

      vecIn[0] = x*factor;
      vecIn[1] = y*factor;
      vecIn[2] = z*factor;
      matVecProduct(vecIn,matT,vecOut);

      fprintf(log,"%5d  %5s",nato+i+1,"NNACP");
      fprintf(log," % 10.6lf % 10.6lf % 10.6lf % 14.8lf % 16.7lf\n",vecOut[0],vecOut[1],vecOut[2],val[0],val[4]+val[7]+val[9]);
    }
  }

  lineDiv(log);
  centraMess("Bond, Ring and Cage critical points",log);

  if( nbc > 0 ) lineDiv(log);
  for(i=0,j=0; i < nCrit; i++){

    if( type[i] == -3 ){
      x = coorCrit[3*i];
      y = coorCrit[3*i+1];
      z = coorCrit[3*i+2];

      id[0] = 0; id[1] = 0; id[2] = 0;
      getData(poly,pts,id,x,y,z,min,hvec,field,coef);
      NumericalCrit02(x,y,z,coef,val);

      vecIn[0] = x*factor;
      vecIn[1] = y*factor;
      vecIn[2] = z*factor;
      matVecProduct(vecIn,matT,vecOut);

      fprintf(log,"%5d  %5s",j+1,"BCP  ");
      fprintf(log," % 10.6lf % 10.6lf % 10.6lf % 14.8lf % 16.7lf\n",vecOut[0],vecOut[1],vecOut[2],val[0],val[4]+val[7]+val[9]);
      j++;
    }
  }

  if( nrc > 0 ) lineDiv(log);
  for(i=0,j=0; i < nCrit; i++){

    if( type[i] == -2 ){
      x = coorCrit[3*i];
      y = coorCrit[3*i+1];
      z = coorCrit[3*i+2];

      id[0] = 0; id[1] = 0; id[2] = 0;
      getData(poly,pts,id,x,y,z,min,hvec,field,coef);
      NumericalCrit02(x,y,z,coef,val);

      vecIn[0] = x*factor;
      vecIn[1] = y*factor;
      vecIn[2] = z*factor;
      matVecProduct(vecIn,matT,vecOut);

      fprintf(log,"%5d  %5s",j+1,"RCP  ");
      fprintf(log," % 10.6lf % 10.6lf % 10.6lf % 14.8lf % 16.7lf\n",vecOut[0],vecOut[1],vecOut[2],val[0],val[4]+val[7]+val[9]);
      j++;
    }
  }

  if( ncc > 0 ) lineDiv(log);
  for(i=0,j=0; i < nCrit; i++){

    if( type[i] == -1 ){
      x = coorCrit[3*i];
      y = coorCrit[3*i+1];
      z = coorCrit[3*i+2];

      id[0] = 0; id[1] = 0; id[2] = 0;
      getData(poly,pts,id,x,y,z,min,hvec,field,coef);
      NumericalCrit02(x,y,z,coef,val);

      vecIn[0] = x*factor;
      vecIn[1] = y*factor;
      vecIn[2] = z*factor;
      matVecProduct(vecIn,matT,vecOut);

      fprintf(log,"%5d  %5s",j+1,"CCP  ");
      fprintf(log," % 10.6lf % 10.6lf % 10.6lf % 14.8lf % 16.7lf\n",vecOut[0],vecOut[1],vecOut[2],val[0],val[4]+val[7]+val[9]);
      j++;
    }
  }
  lineDiv(log);

  double rho,lap,ngrad;
  double kinG,kinK,virial,energyH;
  double l1,l2,l3, matHH[3][3];
  double eigenVec[3][3],eigenVal[3];
 

  //////////////////////////////////////////////////////////////////////////////////////
  if( nna > 0 ){
    centraMess("No-nuclear attractors",log);
    
    for(i=0,j=0; i < nCrit; i++){
      if( type[i] == -4 ){

        lineDiv(log);
        x = coorCrit[3*i];
        y = coorCrit[3*i+1];
        z = coorCrit[3*i+2];

        vecIn[0] = x*factor;
        vecIn[1] = y*factor;
        vecIn[2] = z*factor;
        matVecProduct(vecIn,matT,vecOut);

        id[0] = 0; id[1] = 0; id[2] = 0;
        getData(poly,pts,id,x,y,z,min,hvec,field,coef);
        NumericalCrit02(x,y,z,coef,val);

	     rho   = val[0];
	     ngrad = getNorm(val[1],val[2],val[3]);
	     lap   = val[4] + val[7] + val[9];

	     getEnergies(rho,lap,&kinG,&kinK,&virial);
	     energyH = kinG + virial;
	
        matHH[0][0] = val[4]; matHH[0][1] = val[5]; matHH[0][2] = val[6];
        matHH[1][0] = val[5]; matHH[1][1] = val[7]; matHH[1][2] = val[8];
        matHH[2][0] = val[6]; matHH[2][1] = val[8]; matHH[2][2] = val[9];

        jacobi(matHH,eigenVal,eigenVec);
	     sortJacobi(eigenVal);
	     l1 = eigenVal[0]; l2 = eigenVal[1]; l3 = eigenVal[2];


        fprintf(log,"\n   Critical Point  %7d of type (3,-3)\n\n",j+1);
        fprintf(log,"   Coordinates [A]    % 12.8E  % 12.8E  % 12.8E  \n\n",
               vecOut[0],vecOut[1],vecOut[2]);

	     fprintf(log,"   Density            % 12.8E\n",rho);
	     fprintf(log,"   NormGrad           % 12.8E\n",ngrad);
	     fprintf(log,"   Laplacian          % 12.8E\n\n",lap);

	     fprintf(log,"   Kinetic Energy G   % 12.8E\n",kinG);
	     fprintf(log,"   Kinetic Energy K   % 12.8E\n",kinK);
        fprintf(log,"   Virial field V     % 12.8E\n",virial);
	     fprintf(log,"   Total energy H     % 12.8E\n\n",energyH);

        fprintf(log,"                     | % 12.8E  % 12.8E  %12.8E |\n",matHH[0][0],matHH[0][1],matHH[0][2]);
        fprintf(log,"   Hessian Matrix  = | % 12.8E  % 12.8E  %12.8E |\n",matHH[1][0],matHH[1][1],matHH[1][2]);
        fprintf(log,"                     | % 12.8E  % 12.8E  %12.8E |\n\n",matHH[2][0],matHH[2][1],matHH[2][2]);


	     fprintf(log,"   Eigenvalues        % 12.8E  % 12.8E  % 12.8E\n\n",l1,l2,l3);

        jacobi(matHH,eigenVal,eigenVec);
        matVecProduct(eigenVal,matT,vecOut);
        sortJacobi(vecOut);
        l1 = vecOut[0]; l2 = vecOut[1]; l3 = vecOut[2];
	     fprintf(log,"   Eigenvalues        % 12.8E  % 12.8E  % 12.8E\n\n",l1,l2,l3);


        j++;
      }
    }
  }
  //////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////
  if( nbc > 0 ){
    lineDiv(log);
    centraMess("Bond Critical Points",log);
    
    for(i=0,j=0; i < nCrit; i++){
      if( type[i] == -3 ){
        lineDiv(log);
        x = coorCrit[3*i];
        y = coorCrit[3*i+1];
        z = coorCrit[3*i+2];

        vecIn[0] = x*factor;
        vecIn[1] = y*factor;
        vecIn[2] = z*factor;
        matVecProduct(vecIn,matT,vecOut);

        id[0] = 0; id[1] = 0; id[2] = 0;
        getData(poly,pts,id,x,y,z,min,hvec,field,coef);
        NumericalCrit02(x,y,z,coef,val);

	     rho   = val[0];
	     ngrad = getNorm(val[1],val[2],val[3]);
	     lap   = val[4] + val[7] + val[9];

	     getEnergies(rho,lap,&kinG,&kinK,&virial);
	     energyH = kinG + virial;
	
        matHH[0][0] = val[4]; matHH[0][1] = val[5]; matHH[0][2] = val[6];
        matHH[1][0] = val[5]; matHH[1][1] = val[7]; matHH[1][2] = val[8];
        matHH[2][0] = val[6]; matHH[2][1] = val[8]; matHH[2][2] = val[9];

        jacobi(matHH,eigenVal,eigenVec);
	     sortJacobi(eigenVal);
	     l1 = eigenVal[0]; l2 = eigenVal[1]; l3 = eigenVal[2];


        fprintf(log,"\n   Critical Point  %7d of type (3,-1)\n\n",j+1);
        fprintf(log,"   Coordinates [A]    % 12.8E  % 12.8E  % 12.8E \n\n",
                vecOut[0],vecOut[1],vecOut[2]);

	     fprintf(log,"   Density            % 12.8E\n",rho);
	     fprintf(log,"   NormGrad           % 12.8E\n",ngrad);
	     fprintf(log,"   Laplacian          % 12.8E\n\n",lap);

	     fprintf(log,"   Kinetic Energy G   % 12.8E\n",kinG);
	     fprintf(log,"   Kinetic Energy K   % 12.8E\n",kinK);
        fprintf(log,"   Virial field V     % 12.8E\n",virial);
        fprintf(log,"   Total energy H     % 12.8E\n\n",energyH);

        fprintf(log,"                     | % 12.8E  % 12.8E  %12.8E |\n",matHH[0][0],matHH[0][1],matHH[0][2]);
        fprintf(log,"   Hessian Matrix  = | % 12.8E  % 12.8E  %12.8E |\n",matHH[1][0],matHH[1][1],matHH[1][2]);
        fprintf(log,"                     | % 12.8E  % 12.8E  %12.8E |\n\n",matHH[2][0],matHH[2][1],matHH[2][2]);

	     fprintf(log,"   Eigenvalues        % 12.8E  % 12.8E  % 12.8E\n\n",l1,l2,l3);

	     fprintf(log,"   Bond ellipticity   %12.8lf\n",(l1/l2)-1.);
	     fprintf(log,"   Eta index          % 12.8E\n\n",fabs(l1)/l3);

        jacobi(matHH,eigenVal,eigenVec);
        matVecProduct(eigenVal,matT,vecOut);
        sortJacobi(vecOut);
        l1 = vecOut[0]; l2 = vecOut[1]; l3 = vecOut[2];
	     fprintf(log,"   Eigenvalues        % 12.8E  % 12.8E  % 12.8E\n\n",l1,l2,l3);

	     fprintf(log,"   Bond ellipticity   %12.8lf\n",(l1/l2)-1.);
	     fprintf(log,"   Eta index          % 12.8E\n\n",fabs(l1)/l3);

        j++;
      }
    }
  }
  //////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////
  if( nrc > 0 ){
    lineDiv(log);
    centraMess("Ring Critical Points",log);
    
    for(i=0,j=0; i < nCrit; i++){
      if( type[i] == -2 ){

        lineDiv(log);
        x = coorCrit[3*i];
        y = coorCrit[3*i+1];
        z = coorCrit[3*i+2];

        vecIn[0] = x*factor;
        vecIn[1] = y*factor;
        vecIn[2] = z*factor;
        matVecProduct(vecIn,matT,vecOut);

        id[0] = 0; id[1] = 0; id[2] = 0;
        getData(poly,pts,id,x,y,z,min,hvec,field,coef);
        NumericalCrit02(x,y,z,coef,val);

	     rho   = val[0];
	     ngrad = getNorm(val[1],val[2],val[3]);
	     lap   = val[4] + val[7] + val[9];

	     getEnergies(rho,lap,&kinG,&kinK,&virial);
	     energyH = kinG + virial;
	
        matHH[0][0] = val[4]; matHH[0][1] = val[5]; matHH[0][2] = val[6];
        matHH[1][0] = val[5]; matHH[1][1] = val[7]; matHH[1][2] = val[8];
        matHH[2][0] = val[6]; matHH[2][1] = val[8]; matHH[2][2] = val[9];

        jacobi(matHH,eigenVal,eigenVec);
	     sortJacobi(eigenVal);
        l1 = eigenVal[0]; l2 = eigenVal[1]; l3 = eigenVal[2];


        fprintf(log,"\n   Critical Point  %7d of type (3,1)\n\n",j+1);
        fprintf(log,"   Coordinates [A]    % 12.8E  % 12.8E  % 12.8E  \n\n",
                vecOut[0],vecOut[1],vecOut[2]);

	     fprintf(log,"   Density            % 12.8E\n",rho);
	     fprintf(log,"   NormGrad           % 12.8E\n",ngrad);
	     fprintf(log,"   Laplacian          % 12.8E\n\n",lap);

	     fprintf(log,"   Kinetic Energy G   % 12.8E\n",kinG);
	     fprintf(log,"   Kinetic Energy K   % 12.8E\n",kinK);
	     fprintf(log,"   Virial field V     % 12.8E\n",virial);
	     fprintf(log,"   Total energy H     % 12.8E\n\n",energyH);

        fprintf(log,"                     | % 12.8E  % 12.8E  %12.8E |\n",matHH[0][0],matHH[0][1],matHH[0][2]);
        fprintf(log,"   Hessian Matrix  = | % 12.8E  % 12.8E  %12.8E |\n",matHH[1][0],matHH[1][1],matHH[1][2]);
        fprintf(log,"                     | % 12.8E  % 12.8E  %12.8E |\n\n",matHH[2][0],matHH[2][1],matHH[2][2]);

	     fprintf(log,"   Eigenvalues        % 12.8E  % 12.8E  % 12.8E\n\n",l1,l2,l3);


        jacobi(matHH,eigenVal,eigenVec);
        matVecProduct(eigenVal,matT,vecOut);
        sortJacobi(vecOut);
        l1 = vecOut[0]; l2 = vecOut[1]; l3 = vecOut[2];
	     fprintf(log,"   Eigenvalues        % 12.8E  % 12.8E  % 12.8E\n\n",l1,l2,l3);

        j++;
      }
    }
  }
  //////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////
  if( ncc > 0 ){
    lineDiv(log);
    centraMess("Cage Critical Points",log);
    
    for(i=0,j=0; i < nCrit; i++){
      if( type[i] == -1 ){

        lineDiv(log);
        x = coorCrit[3*i];
        y = coorCrit[3*i+1];
        z = coorCrit[3*i+2];

        vecIn[0] = x*factor;
        vecIn[1] = y*factor;
        vecIn[2] = z*factor;
        matVecProduct(vecIn,matT,vecOut);

        id[0] = 0; id[1] = 0; id[2] = 0;
        getData(poly,pts,id,x,y,z,min,hvec,field,coef);
        NumericalCrit02(x,y,z,coef,val);

	     rho   = val[0];
	     ngrad = getNorm(val[1],val[2],val[3]);
	     lap   = val[4] + val[7] + val[9];

	     getEnergies(rho,lap,&kinG,&kinK,&virial);
	     energyH = kinG + virial;
	
        matHH[0][0] = val[4]; matHH[0][1] = val[5]; matHH[0][2] = val[6];
        matHH[1][0] = val[5]; matHH[1][1] = val[7]; matHH[1][2] = val[8];
        matHH[2][0] = val[6]; matHH[2][1] = val[8]; matHH[2][2] = val[9];

        jacobi(matHH,eigenVal,eigenVec);
	     sortJacobi(eigenVal);
        l1 = eigenVal[0]; l2 = eigenVal[1]; l3 = eigenVal[2];


        fprintf(log,"\n   Critical Point  %7d of type (3,3)\n\n",j+1);
        fprintf(log,"   Coordinates [A]    % 12.8E  % 12.8E  % 12.8E  \n\n",
                vecOut[0],vecOut[1],vecOut[2]);

	     fprintf(log,"   Density            % 12.8E\n",rho);
	     fprintf(log,"   NormGrad           % 12.8E\n",ngrad);
	     fprintf(log,"   Laplacian          % 12.8E\n\n",lap);

	     fprintf(log,"   Kinetic Energy G   % 12.8E\n",kinG);
	     fprintf(log,"   Kinetic Energy K   % 12.8E\n",kinK);
	     fprintf(log,"   Virial field V     % 12.8E\n",virial);
	     fprintf(log,"   Total energy H     % 12.8E\n\n",energyH);

        fprintf(log,"                     | % 12.8E  % 12.8E  %12.8E |\n",matHH[0][0],matHH[0][1],matHH[0][2]);
        fprintf(log,"   Hessian Matrix  = | % 12.8E  % 12.8E  %12.8E |\n",matHH[1][0],matHH[1][1],matHH[1][2]);
        fprintf(log,"                     | % 12.8E  % 12.8E  %12.8E |\n\n",matHH[2][0],matHH[2][1],matHH[2][2]);

	     fprintf(log,"   Eigenvalues        % 12.8E  % 12.8E  % 12.8E\n\n",l1,l2,l3);

        j++;
      }
    }
  }
  //////////////////////////////////////////////////////////////////////////////////////
  
  lineDiv(log);
  centraMess("End of File",log);
  lineDiv(log);
  
  fclose(log);

  return 0;
}

void sortJacobi(double vec[3]){

  int i,j;
  double tmp;

  for(i=0; i < 3; i++){
    for(j=i+1; j < 3; j++){
      if(vec[i] > vec[j] ){
        tmp = vec[i];
	vec[i] = vec[j];
	vec[j] = tmp;
      }
    }
  }



}

void getEnergies(double rho, double lap, double *kinG,double *kinK,double *vir){

  double kG,kK;

  kG = (3./10.)*pow(3.*M_PI*M_PI,2./3.)*pow(rho,5./3.);
  kG += (1./6.)*lap;

  kK = kG - 0.25*lap;


  (*vir) = 0.25*lap - 2.*kG;

  (*kinG) = kG;
  (*kinK) = kK;


}


void printLog1(FILE *log,int ntyp[5]){
  

 
  char tchar[128];

  time_t t= time(NULL);
  struct tm *tlocal = localtime(&t);

  int ncc,nrc,nbc,nna,nnu;

  ncc = ntyp[0];
  nrc = ntyp[1];
  nbc = ntyp[2];
  nna = ntyp[3];
  nnu = ntyp[4];

  strftime(tchar,128,"%H:%M:%S  %d/%m/%y",tlocal);


  lineDiv(log);

  centraMess("File log for cube3d",log);
  lineDiv(log);
  fprintf(log,"%82s\n",tchar);
  lineDiv(log);
  
  fprintf(log," This file contains information about critical points and properties, these\n");
  fprintf(log," information are determineted by code cube3d v.0\n");
  fprintf(log," The units in this  file for  distance are  Angstrom and  the units for fields are\n");
  fprintf(log," atomic units.\n\n");
  fprintf(log,"  NNACP               Non-nuclear attractor critical point\n");
  fprintf(log,"  BCP                 Bond critical point\n");
  fprintf(log,"  RCP                 Ring critical point\n");
  fprintf(log,"  CCP                 Cage critical point\n");
  fprintf(log,"  Kinetic Energy G    Kinetic Energy Density in the form Lagrangian\n");
  fprintf(log,"  Kinetic Energy K    Kinetic Energy Density in the form Hamiltonian\n");
  fprintf(log,"  Virial Field        Or Potential Energy Density  V = -K - G\n");
  //fprintf(log,"  Lagrangian Density  -(1/4) Laplacian or K - G\n");
 // fprintf(log,"  KenergyG/Density    Kinetic Energy Density per electron\n");
  
  fprintf(log,"  Eigenvalues         Eigenvalues for Hessian matrix: lambda_1<lambda_2<lambda_3\n");
  fprintf(log,"  Bond Ellipiticy     lambda_1/lambda_2 - 1\n");
  fprintf(log,"  Eta index           |lambda_1|/lambda_3\n\n");
  lineDiv(log);
  centraMess("Information about the system",log);
  lineDiv(log);
  fprintf(log,"    Nuclei                : %15d\n",nnu);
  fprintf(log,"    NNACP   (3,-3)        : %15d\n",nna);
  fprintf(log,"    BCP     (3,-1)        : %15d\n",nbc);
  fprintf(log,"    RCP     (3,+1)        : %15d\n",nrc);
  fprintf(log,"    CCP     (3,+3)        : %15d\n",ncc);
  fprintf(log,"    Total Critical Points : %15d\n",nna+nbc+nrc+ncc);
  lineDiv(log);
  fprintf(log,"    Poincare-Hopf rule\n");
  fprintf(log,"    nuclei + NNACP + RCP - BCP - CCP = %5d\n",nnu+nna+nrc-nbc-ncc);
  lineDiv(log);

  centraMess("Nuclear and Non-nuclear Attractors",log);
  lineDiv(log);
  fprintf(log,"  Atractor                Coordinates [A]            Density         Laplacian\n");
  lineDiv(log);



}

void centraMess(char *mess,FILE *out){
  int nmax = 82;
  int i, n = strlen(mess);
  int spacio;
  if( n%2 == 1)
    n++;
  if(nmax%2 == 1)
    nmax++;

  spacio = nmax/2 - n/2;

  if( spacio > 0 )
    for(i=0;i<spacio;i++)
      fprintf(out," ");
  fprintf(out,"%s\n",mess);


}


void lineDiv(FILE *out){
  int i;
  
  for(i=0;i<82;i++)
    fprintf(out,"=");
  
  fprintf(out,"\n");


}

