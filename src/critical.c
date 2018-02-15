#include "file.h"
#include "array.h"
#include "analys.h"
#include "tableP.h"
#include "matvec.h"
#include "critical.h"
#include "pruebaCoef.h"
#include "numBondPath.h"
#include "lagrangePol.h"

int getIndex(char c,int ptq, double q, double q0, double hq){
  int i;
  int index;
  double qi;
  double qmin;
  index = -1;
  qmin = 1.E10;
  for(i=0; i < ptq; i++){
    qi = q0 + i*hq;
    if( fabs( q - qi) < qmin){
      qmin = fabs( q - qi);
      index = i;
    }
  }

  return  index;
}


int getCube(int *pts,double x, double y, double z, double *min,double *h,int *i,int *j, int *k){

  int itmp,jtmp,ktmp;

  itmp = getIndex('x',pts[0],x,min[0],h[0]);
  jtmp = getIndex('y',pts[1],y,min[1],h[1]);
  ktmp = getIndex('z',pts[2],z,min[2],h[2]);

  (*i) = itmp;
  (*j) = jtmp;
  (*k) = ktmp;

  
  return 0;
}


int findCriticalPoints(int pol,dataCube data1,int *zato, double *coor, double *field, double *matT, const char *name){

  int i,j,k;
  int npx,npy,npz;
  int mu, cubos;
  int *idx;
  double fun2[27];
  double xi,yi,zi;
  double x1,y1,z1;
  double hx,hy,hz;
  double *xyz;

  FILE *out;
  char nameCrit[128];
  sprintf(nameCrit,"%sCrit.xyz",name);
  openFile(&out,nameCrit,"w+");

  xi = data1.min[0];  yi = data1.min[1];  zi = data1.min[2];
  hx = data1.hvec[0]; hy = data1.hvec[1]; hz = data1.hvec[2];
 
  npx = data1.pts[0];
  npy = data1.pts[1];
  npz = data1.pts[2];

  cubos = (int) (npx-1)*(npy-1)*(npz-1);

  createArrayInt(cubos,&idx,"idxCube");
  createArrayDou(3*cubos,&xyz,"xyzCritico");
  for(mu=0 ; mu<cubos ; mu++)
    idx[mu] = -1;

  for( i = 1; i < npx-1 ; i+=2){
    x1 = xi + i*hx;
    for( j = 1; j < npy-1 ; j+=2){
      y1 = yi + j*hy;
      for( k = 1; k < npz-1 ; k+=2){
        z1 = zi + k*hz;

        loadField(2,i,j,k,npx,npy,npz,field,fun2);

        findCritical(i,j,k,data1.pts,idx,hx,hy,hz,x1,y1,z1,xyz,fun2);
      }
    }
  }  

  refineCritical(data1,pol,idx,coor,field,xyz,matT,out,name);

  return 0;
}

double getDist2(double x1, double y1, double z1, double x2, double y2, double z2){

  return (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1);

}

double getNorm(double vx, double vy, double vz){

  return sqrt( vx*vx + vy*vy + vz*vz);

}

int findCritical(int i, int j, int k, int *pts, int *idx,
                 double hx, double hy, double hz,
                 double x1, double y1, double z1, double *xyz,double *f){

   int ncy,ncz;
   double c[27];
   double hmax;
   double x0,x2,y0,y2,z0,z2;
   double limits[6];
   //double hvec[3];

   //hvec[0] = hx;
   //hvec[1] = hy;
   //hvec[2] = hz;

   ncy = pts[1] - 1;
   ncz = pts[2] - 1;

   coefficients(hx,hy,hz,x1,y1,z1,f,c);
   //getCoeff(2,hvec,x1,y1,z1,f,c);

   limits[0] = x1 - hx;
   limits[1] = x1 + hx;
   limits[2] = y1 - hy;
   limits[3] = y1 + hy;
   limits[4] = z1 - hz;
   limits[5] = z1 + hz;

   hx *= 0.5;
   hy *= 0.5;
   hz *= 0.5;
   hmax = mayor3(hx,hy,hz);

   x0 = x1 - hx; x2 = x1 + hx;
   y0 = y1 - hy; y2 = y1 + hy;
   z0 = z1 - hz; z2 = z1 + hz;

   cubeFail(i-1,j-1,k-1,ncy,ncz,idx,x0,y0,z0,hmax,xyz,limits,c);
   cubeFail(i-1,j-1,k+1,ncy,ncz,idx,x0,y0,z2,hmax,xyz,limits,c);
   cubeFail(i-1,j+1,k-1,ncy,ncz,idx,x0,y2,z0,hmax,xyz,limits,c);
   cubeFail(i-1,j+1,k+1,ncy,ncz,idx,x0,y2,z2,hmax,xyz,limits,c);
   cubeFail(i+1,j-1,k-1,ncy,ncz,idx,x2,y0,z0,hmax,xyz,limits,c);
   cubeFail(i+1,j-1,k+1,ncy,ncz,idx,x2,y0,z2,hmax,xyz,limits,c);
   cubeFail(i+1,j+1,k-1,ncy,ncz,idx,x2,y2,z0,hmax,xyz,limits,c);
   cubeFail(i+1,j+1,k+1,ncy,ncz,idx,x2,y2,z2,hmax,xyz,limits,c);
   

   return 0;
}

int cubeFail( int i, int j, int k, int ncy, int ncz, int *idx,
              double xi,double yi, double zi,double rmax,
              double *xyz,double *limits,double *c){
  
  int indice,iter;
  int difcoor;
  double val[10];
  double f0,norm;
  double x,y,z;
  double xn,yn,zn;

  
  indice = i*ncy*ncz + j*ncz + k;

  if( idx[indice] != 1 ){

  //NumericalCrit01(xi,yi,zi,c,val);
  evalPol1(2,xi,yi,zi,c,val);

  f0 = val[0];
  norm = getNorm(val[1],val[2],val[3]);

  if( fabs(f0) < TOLFUN && norm < 1.E-4){
  //if( fabs(f0) < TOLFUN ){
    idx[indice] = 0;
  }else{
  
    idx[indice] = 0;
    iter = 0;
    x = xi; y = yi; z = zi;
    norm = 1.E3;
    difcoor = 0;

    while( (iter < MAXITER) && (norm > TOLGRD) && (difcoor == 0) ){
      //NumericalCrit02(x,y,z,c,val);
      evalPol2 (2,x,y,z,c,val);
      f0 = val[0];
      norm = getNorm(val[1],val[2],val[3]);

      if( fabs(f0) < TOLFUN){
        iter = 2*MAXITER;
        x = xi + 1.E3;
        y = yi + 1.E3;
        z = zi + 1.E3;
        idx[indice] = 0;
      }else{

        hessGrad(&xn,&yn,&zn,&norm,val);
        x += xn;
        y += yn;
        z += zn;
        iter++;
      }
      difcoor = inCube(x,y,z,limits);
    } // fin while

    if( difcoor == 0 && norm <= TOLNRM ){
      idx[indice] = 1;
      //printf(" >> indice uno %2d %2d %10d %10.6E \n",inCube(x,y,z,limits),idx[indice],indice,norm);
      //printf(" He  % 10.6lf % 10.6lf % 10.6lf\n",xi,yi,zi);
      xyz[3*indice]   = x;
      xyz[3*indice+1] = y;
      xyz[3*indice+2] = z;
    }


    }// fin del if-else
  }  
  
  return 0;
}

int signature(double *val, int *r, int *s){
  int i;
  int rint,sint;
  double tmp;

  rint = sint = 0;
  for(i = 0; i < 3; i++){
    tmp = val[i];
    if( tmp != 0. ){
      rint++;
      if( tmp < 0.) 
        sint --; 
      else 
        sint++;
    }
  }
  
  (*r) = rint;
  (*s) = sint;


  return 0;
}

int inCube(double x, double y, double z, double limits[6]){
  
  // devolvera 0 si estamos dentro del macrocubo
  // devolvera 1 si estamos fuera del macrocubo
  // las distancias del macrocubo son dadas por limits[6];

  int flag = 0;

  if( limits[0] < x && limits[1] > x )
    flag += 1;
  if( limits[2] < y && limits[3] > y )
    flag += 2;
  if( limits[4] < z && limits[5] > z )
    flag += 4;

  if( flag == 7 )
    flag = 0;
  else 
    flag = 1;

  return flag;


}

double mayor3( double a, double b, double c){
  double max;

  if( a > b && a > c)
    max = a;
  else{

    if( b > a && b > c)
      max = b;
    else
      max = c;
  }

  return max;

}

double distance(double x1, double x2, double y1, double y2, double z1, double z2){

  double r;
  r  = (x2-x1)*(x2-x1);
  r += (y2-y1)*(y2-y1);
  r += (z2-z1)*(z2-z1);

  return sqrt(r);
}

double hessGrad(double *x, double *y, double *z,double *ngrad,double *val){
  double vx,vy,vz;
  double a,b,c,d,e,g;
  double dis,dis1;
  double t1,t2,t3,t4,t5,t6;

  /****************************************************************************
   *         | a  b  c |
   * mat H = | b  d  e |
   *         | c  e  g |
   *
   *           | t1 t2 t3 |
   * mat H-1 = | t2 t4 t5 |
   *           | t3 t5 t6 |
   *
   * (x,y,z) = mat H-1 . vecGrad(vx,vy,vz)
   *
   *****************************************************************************/
  
  
  vx = val[1];
  vy = val[2];
  vz = val[3];

  (*ngrad) = getNorm(vx,vy,vz);

  a = val[4];
  b = val[5];
  c = val[6];
  d = val[7];
  e = val[8];
  g = val[9];

  dis = 2.*b*c*e - c*c*d - a*e*e - b*b*g + a*d*g;
  dis1 = -1./dis;

  t1 = ( d*g - e*e)*dis1;
  t2 = ( c*e - b*g)*dis1;
  t3 = ( b*e - c*d)*dis1;
  t4 = ( a*g - c*c)*dis1;
  t5 = ( b*c - a*e)*dis1;
  t6 = ( a*d - b*b)*dis1;
  
  (*x) = t1*vx + t2*vy + t3*vz;
  (*y) = t2*vx + t4*vy + t5*vz;
  (*z) = t3*vx + t5*vy + t6*vz;

  
  
  return dis;
}



int refineCritical( dataCube data1, int pol,int *idx,
                    double *coor,double *field,double *xyz,double *matT,FILE *out,const char *name){

  int j,k,flag;
  int npx,npy,npz;
  int indi,indj,indk;
  int mu,nu,cubos;
  int size;

  double xat,yat,zat;
  double fac = 0.52917;
  double hvec[3];
  double *fun,*coef;
  double x1,y1,z1;
  double x,y,z;
  double xn,yn,zn;
  double limits[6];
  double rij;
  int iter=0; 
  double norm = 1.E3;
  int difcoor=0;
  double matA[6],val[10];
  double eigval[3];
  int r,s,count;
  int ncp,bcp,rcp,ccp;
  int *type;
  double *coorCrit;

  npx = data1.pts[0];
  npy = data1.pts[1];
  npz = data1.pts[2];

  hvec[0] = data1.hvec[0];
  hvec[1] = data1.hvec[1];
  hvec[2] = data1.hvec[2];
  
  size = (int) (pol+1)*(pol+1)*(pol+1);
  createArrayDou(size,&fun,"funct");
  createArrayDou(size,&coef,"coef");

  cubos = (int) (npx-1)*(npy-1)*(npz-1);
  
  count = 0;
  for(mu=0; mu < cubos; mu++){
    if(idx[mu] == 1 ){
      count++;
    }
  }


  // Vamos a quitar los valores repetidos
  double *aux;
  double *xyzOut;
  double *xyzInp;
  createArrayDou(count,&aux,"Auxiliar");
  createArrayDou(3*count,&xyzOut,"Auxlilar");
  createArrayDou(3*count,&xyzInp,"Auxlilar");
  
  nu=0;
  for(mu = 0; mu < cubos; mu++){
    if(idx[mu] == 1){
      xyzInp[3*nu]   = xyz[3*mu];
      xyzInp[3*nu+1] = xyz[3*mu+1];
      xyzInp[3*nu+2] = xyz[3*mu+2];
      nu++;
    }
  }
  printf(" Possible critical points : %10d\n",count);

  //int contador;
  int i;
  //double num;
  printf("=====================================================\n");
  printf(" %4d\n",count);
  printf(" Posibles puntos \n");
  for(i=0;i<count;i++)
    printf(" POS  % 10.6lf % 10.6lf % 10.6lf\n",xyzInp[3*i]*fac,
                                                xyzInp[3*i+1]*fac,
                                                xyzInp[3*i+2]*fac);
  printf("=====================================================\n");

  /*nu = 0; i = 0;

  for(mu=0;mu<count;mu++){
    contador = 0;
    num = getDist2(xyzInp[3*mu],xyzInp[3*mu+1],xyzInp[3*mu+2],5.,-10.,-5.);
    //printf(" mu % 2d   num %20.10E\n",mu,num);

    aux[nu] = num;
    nu++;
    for(k=0;k<count;k++){
      if( fabs(aux[k]-num) < 1.E-7)
        contador++;
    }

    if( contador == 1){
      xyzOut[3*i]   = xyzInp[3*mu];
      xyzOut[3*i+1] = xyzInp[3*mu+1];
      xyzOut[3*i+2] = xyzInp[3*mu+2];
      i++;
    }

    

  }

  */
  i = deleteRepeated(count, xyzInp,xyzOut);
  printf(" Despues : %d\n",i);
  count = i;


  printf("=====================================================\n");
  printf(" %4d\n",count);
  printf(" Posibles puntos Despues\n");
  for(i=0;i<count;i++)
    printf(" POS  % 10.6lf % 10.6lf % 10.6lf\n",xyzOut[3*i]*fac,
                                                xyzOut[3*i+1]*fac,
                                                xyzOut[3*i+2]*fac);
  printf("=====================================================\n");

  free(aux);
  free(xyzInp);
  //exit(EXIT_FAILURE);
  // DEspues de quitar los repetidos

  createArrayInt(count,&type,"typeCritical");
  createArrayDou(3*count,&coorCrit,"CoorCrit");



  ncp = bcp = rcp = ccp = 0;
  k= 0;
  printf("=====================================================\n");
  for(mu=0; mu < count; mu++){
    flag = 1;
    x = xyzOut[3*mu];
    y = xyzOut[3*mu+1];
    z = xyzOut[3*mu+2];

    getCube(data1.pts,x,y,z,data1.min,hvec,&indi,&indj,&indk);
    loadField(pol,indi,indj,indk,npx,npy,npz,field,fun);

    x1 = data1.min[0] + indi*hvec[0];
    y1 = data1.min[1] + indj*hvec[1];
    z1 = data1.min[2] + indk*hvec[2];

    //if(indi==0) indi++;
    //if(indj==0) indj++;
    //if(indk==0) indk++;

    //limits[0] = data1.min[0] + (indi-1)*hvec[0];
    //limits[1] = data1.min[0] + (indi+1)*hvec[0];
    //limits[2] = data1.min[1] + (indj-1)*hvec[1];
    //limits[3] = data1.min[1] + (indj+1)*hvec[1];
    //limits[4] = data1.min[2] + (indk-1)*hvec[2];
    //limits[5] = data1.min[2] + (indk+1)*hvec[2];

    limitsMacro(pol,indi,indj,indk,data1,limits);

    coefficients(hvec[0],hvec[1],hvec[2],x1,y1,z1,fun,coef);
//    getCoeff(pol,hvec,x1,y1,z1,fun,coef);

    iter=0;  norm = 1.E3;  difcoor=0;

    while( iter < 3*MAXITER && norm > TOLGRD2 && difcoor == 0){
      NumericalCrit02(x,y,z,coef,val);
      //evalPol2(pol,x,y,z,coef,val);
      hessGrad(&xn,&yn,&zn,&norm,val);
      x += xn;
      y += yn;
      z += zn;
      iter++;
      difcoor = inCube(x,y,z,limits);
    }

    flag = !inCube(x,y,z,limits);

    if( flag == 1){
      type[k]    = 0;
      coorCrit[3*k]   = x;
      coorCrit[3*k+1] = y;
      coorCrit[3*k+2] = z;

      evalPol2(pol,x,y,z,coef,val);

      matA[0] = val[4];
      matA[1] = val[5];
      matA[2] = val[6];
      matA[3] = val[7];
      matA[4] = val[8];
      matA[5] = val[9];

      valoresPropios3x3(matA,eigval);

      signature(eigval,&r,&s);

      for(j=0; j < data1.natm; j++){
        xat = coor[3*j];
        yat = coor[3*j+1];
        zat = coor[3*j+2];
        rij = distance(x,xat,y,yat,z,zat);
        if( rij < TOLDIS ){
          rij = -10;
          j = data1.natm+1;
        }
          
      }
  
     double ngrad = getNorm(val[1],val[2],val[3]);

     if( (fabs(val[0]) > CUTRHO) && ( ngrad < TOLNRM) && (rij > 0) ){
     //if( (fabs(val[0]) > CUTRHO) && ( ngrad < 1.E-5) && (rij > 0) ){

      if( s == -3 ){
        type[k] = -4;
        ncp++;
      }
      if( s == -1 ){
        type[k] = -3;
        bcp++;
      }
      if( s == 1 ){
        type[k] = -2;
        rcp++;
      }
      if( s == 3 ){
        type[k] = -1;
        ccp++;
      }
      k++;
    }
  }

}

  
  int ncrit;
  ncrit = SortCoordinates(count,type,coorCrit);
  printf(" --> Posibles puntos criticos: % 5d\n",count);
  printf(" --> Reales puntos criticos  : % 5d\n",ncrit);
  printf(" --> Puntos criticos No nuc  : % 5d\n",ncp);
  printf(" --> Puntos criticos enlace  : % 5d\n",bcp);
  printf(" --> Puntos criticos anillo  : % 5d\n",rcp);
  printf(" --> Puntos criticos caja    : % 5d\n",ccp);

  printf("=====================================================\n");

  fprintf(out," % 3d\n",data1.natm+ncrit);
  fprintf(out," File crit in format XYZ\n");
  char symb[4];
  extern int *zatm;
  double vecIn[3],vecTmp[3];
  for( i = 0; i < data1.natm; i++){
    getAtomicSymbol(zatm[i],4,symb);
    vecIn[0] = coor[3*i  ]*fac;
    vecIn[1] = coor[3*i+1]*fac;
    vecIn[2] = coor[3*i+2]*fac;
    matVecProduct(vecIn,matT,vecTmp);
    fprintf(out,"   %-8s % 10.6lf % 10.6lf % 10.6lf\n",symb,vecTmp[0],
                                                            vecTmp[1],
                                                            vecTmp[2]);

  }


  for( i = 0; i < ncrit; i++ ){
    switch(type[i]){
      case  -1: fprintf(out,"   %-8s","CCP"); break;
      case  -2: fprintf(out,"   %-8s","RCP"); break;
      case  -3: fprintf(out,"   %-8s","BCP"); break;
      case  -4: fprintf(out,"   %-8s","NNACP"); break;
    }
    vecIn[0] = coorCrit[3*i  ]*fac;
    vecIn[1] = coorCrit[3*i+1]*fac;
    vecIn[2] = coorCrit[3*i+2]*fac;
    matVecProduct(vecIn,matT,vecTmp);

    fprintf(out," % 10.6lf % 10.6lf % 10.6lf\n",vecTmp[0],vecTmp[1],vecTmp[2]);

  }

  fclose(out);
  
  int *bonding;
  createArrayInt(2*bcp,&bonding,"bonding");

  bondPath(pol,ncp,data1.natm,ncrit,type,bonding,coor,coorCrit,data1.pts,data1.min,hvec,field,matT,name);

  logFile (pol,data1.natm,ncrit,type,bonding,coor,coorCrit,data1.pts,data1.min,hvec,field,matT,name);


  printf("=====================================================\n");

  free(coef);
  free(fun);
  free(type);
  
  
  return 0;
}

int deleteRepeated(int n, double *arrayInp, double *arrayOut){

  int i,j,k,m;
  int cont;
  double vecAux[3*n];
  double tmp[3];

  j=m=0;
  for(i=0;i<n;i++){
    cont = 0;
    tmp[0]     = arrayInp[3*i];
    tmp[1]     = arrayInp[3*i+1];
    tmp[2]     = arrayInp[3*i+2];
    vecAux[3*j]   = tmp[0];
    vecAux[3*j+1] = tmp[1];
    vecAux[3*j+2] = tmp[2];
    j++;
    for(k=0;k<j;k++)
      if( fabs(vecAux[3*k]   - tmp[0]) < 1.E-6 && 
          fabs(vecAux[3*k+1] - tmp[1]) < 1.E-6 && 
          fabs(vecAux[3*k+2] - tmp[2]) < 1.E-6 ){
        cont++;
      }   
    if( cont == 1){ 
      arrayOut[3*m  ] = tmp[0];
      arrayOut[3*m+1] = tmp[1];
      arrayOut[3*m+2] = tmp[2];
      m++;
    }   


  }
  printf(" Entraron : %d \n",n);
  printf(" Salieron : %d \n",m);

  return m;

}
