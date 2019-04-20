/**
 * @file   lagrange.c
 * @brief 
 * @author Raymundo Hernández-Esparza.
 * @date   August 2018.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lagrange2.h"
#include "transU.h"

#define YES 1
#define NOT 0

int getData(double x,int n,int pol,double *xcon, double *fun,double *xpol, double *fpol){

  int i,dat;
  int index0;

  dat = pol + 1;

  for(i=0;i<n-1;i++){
    if( x >= xcon[i] && x < xcon[i+1])
      index0 = i;
  }

  printf(" GetData function  Polinomio : %4d  Ndata %4d  Index0 = %4d\n",pol,dat,index0);
  int min,max;
  switch (dat){
    case  2: min =  0; max =  1; break;
    case  3: min = -1; max =  1; break;

    case  4: min = -1; max =  2; break;
    case  5: min = -2; max =  2; break;

    case  6: min = -2; max =  3; break;
    case  7: min = -3; max =  3; break;

    case  8: min = -3; max =  4; break;
    case  9: min = -4; max =  4; break;

    case 10: min = -4; max =  5; break;
    case 11: min = -5; max =  5; break;

    case 12: min = -5; max =  6; break;
    case 13: min = -6; max =  6; break;

    case 14: min = -6; max =  7; break;
    case 15: min = -7; max =  7; break;

    case 16: min = -7; max =  8; break;
    case 17: min = -8; max =  8; break;

    case 18: min = -8; max =  9; break;
    case 19: min = -9; max =  9; break;

    case 20: min = - 9; max =10; break;
    case 21: min = -10; max =10; break;
  }

  min = index0 + min;
  max = index0 + max;

  int k = 0;
  for(i=min;i<=max;i++){
    xpol[k] = xcon[i];
    fpol[k] = fun [i];
    printf(" % 10.6lf ",xpol[k]);
    k++;
  }
  printf("\n");
  return 0;
}

int getLagParameters(double x0, double *xcon, double *xxmu, double *delta, int n){
  
  int i,j;

  for( i = 0 ; i <=n ; i++ )
    xxmu[i] = x0 - xcon[i];

  for( i = 0 ; i <=n ; i++ ){
    delta[i] = (double) 1.;
    for( j = 0 ; j <=n ; j++ )
      if( j != i )
        delta[i] *= (xcon[i] - xcon[j]);
  }


  return 0;
}
/*****************************************************************************/
double lagrange(double *xxmu, double delta, int k, int n){

  int mu;
  double valor;

  valor = (double) 1.;
  
  for( mu = 0 ; mu <= n ; mu++ ){
    if( mu != k){
      valor *= xxmu[mu];
    }
  }

  return valor/delta;
}

double lagrange1(double *xxmu, double delta, int k, int n){

  int nu,mu;
  double sum,prod;

  sum = (double) 0.;
  
  for( nu = 0; nu <= n ; nu++){
    prod = (double) 0.;
    if(nu != k){
      prod = (double) 1.;
      for( mu = 0; mu <= n ; mu++){
        if( mu != k && mu != nu){
          prod *= xxmu[mu];
        }
      }
    }
    sum += prod;
  }
  

  return sum/delta;
}

double lagrange2(double *xxmu, double delta, int k, int n){

  int nu,mu,la;
  double sum,prod;

  sum = (double) 0.;
  
  for( nu = 0; nu <= n ; nu++){
   for( la = 0; la <= n ; la++){
    if( la != k && la != nu){
    prod = (double) 0.;
    if(nu != k){
      prod = (double) 1.;
      for( mu = 0; mu <= n ; mu++){
        if( mu != k && mu != nu && mu != la){
          prod *= xxmu[mu];
        }
      }
     }
    sum += prod;
    }
   }
  }
  

  return sum/delta;
}
/*****************************************************************************/
double interpolador1D(double x0, double *xcon, double *f,int n){
  
  int i; 
  double valor = (double) 0.;
  double xxmu[n+1];
  double delta[n+1];

  printf(" Valor inicial: % 10.8lf\n",xcon[0]);
 
  getLagParameters(x0,xcon,xxmu,delta, n);
  
  for( i = 0 ; i <=n ; i++ )
    valor += f[i]*lagrange(xxmu,delta[i],i,n);


  return valor;
}
double interpolador1Dder1(double x0, double *xcon, double *f,int n){

  int i; 
  double valor = (double) 0.;
  double xxmu[n+1];
  double delta[n+1];

  getLagParameters(x0,xcon,xxmu,delta, n);

  for( i = 0 ; i <=n ; i++ )
    valor += f[i]*lagrange1(xxmu,delta[i],i,n);

  return valor;
}

double interpolador1Dder2(double x0, double *xcon, double *f,int n){

  int i; 
  double valor = (double) 0.;
  double xxmu[n+1];
  double delta[n+1];

  getLagParameters(x0,xcon,xxmu,delta, n);

  for( i = 0 ; i <=n ; i++ )
    valor += f[i]*lagrange2(xxmu,delta[i],i,n);


  return valor;
}

int getDerivatives1D(double x0, double *xcon, double *f,int n,double *val){

  int i; 
  double valor0,valor1,valor2;
  double xxmu[n+1];
  double delta[n+1];

  valor0 = valor1 = valor2 = 0.;


  getLagParameters(x0,xcon,xxmu,delta, n);

  for( i = 0 ; i <=n ; i++ ){
    valor0 += f[i]*lagrange (xxmu,delta[i],i,n);
    valor1 += f[i]*lagrange1(xxmu,delta[i],i,n);
    valor2 += f[i]*lagrange2(xxmu,delta[i],i,n);
  }

  val[0] = valor0;
  val[1] = valor1;
  val[2] = valor2;

  return 0;
}
/*****************************************************************************/

double interpolador3D(double    x0, double    y0, double    z0, 
                      double *xcon, double *ycon, double *zcon,
                      double  *fun, int n){
  int index;
  int m = n+1;
  int i,j,k;
  double lx[n+1],ly[n+1],lz[n+1];
  double valor = (double) 0.;
  double xxmu[n+1];
  double yymu[n+1];
  double zzmu[n+1];
  double deltax[n+1];
  double deltay[n+1];
  double deltaz[n+1];

  getLagParameters(x0,xcon,xxmu,deltax,n);
  getLagParameters(y0,ycon,yymu,deltay,n);
  getLagParameters(z0,zcon,zzmu,deltaz,n);

  for( i = 0 ; i <= n ; i++ ){
    lx[i] = lagrange(xxmu,deltax[i],i,n);
    ly[i] = lagrange(yymu,deltay[i],i,n);
    lz[i] = lagrange(zzmu,deltaz[i],i,n);
  }

  valor = (double) 0.;
  for( i = 0 ; i <=n ; i++ )
    for( j = 0 ; j <=n ; j++ )
      for( k = 0 ; k <=n ; k++ ){
        index = i*m*m + j*m + k;
        
        valor += fun[index]*lx[i]*ly[j]*lz[k];
      }
 
  return valor;
}

int gradient3D(double    x0, double    y0, double    z0, 
                     double *xcon, double *ycon, double *zcon,
                     double  *fun, int n ,int trans, const double *matU, double *val){
  int index;
  int m = n+1;
  int i,j,k;
  double l0x[n+1],l0y[n+1],l0z[n+1];
  double l1x[n+1],l1y[n+1],l1z[n+1];

  double f,val0;
  double valx,valy,valz;
  double xxmu[n+1];
  double yymu[n+1];
  double zzmu[n+1];
  double deltax[n+1];
  double deltay[n+1];
  double deltaz[n+1];


  getLagParameters(x0,xcon,xxmu,deltax,n);
  getLagParameters(y0,ycon,yymu,deltay,n);
  getLagParameters(z0,zcon,zzmu,deltaz,n);

  for( i = 0 ; i <= n ; i++ ){
    l0x[i] = lagrange (xxmu,deltax[i],i,n);
    l0y[i] = lagrange (yymu,deltay[i],i,n);
    l0z[i] = lagrange (zzmu,deltaz[i],i,n);

    l1x[i] = lagrange1(xxmu,deltax[i],i,n);
    l1y[i] = lagrange1(yymu,deltay[i],i,n);
    l1z[i] = lagrange1(zzmu,deltaz[i],i,n);
  }

  val0  = (double) 0.;
  valx  = valy  = valz  = (double) 0.;

  for( i = 0 ; i <=n ; i++ )
    for( j = 0 ; j <=n ; j++ )
      for( k = 0 ; k <=n ; k++ ){
        index = i*m*m + j*m + k;
        f = fun[index];
        
        val0  += f*l0x[i]*l0y[j]*l0z[k];

        valx  += f*l1x[i]*l0y[j]*l0z[k];
        valy  += f*l0x[i]*l1y[j]*l0z[k];
        valz  += f*l0x[i]*l0y[j]*l1z[k];

      }

  val[0] = val0;
  val[1] = valx;
  val[2] = valy;
  val[3] = valz;

  if(trans != YES )
    trans01(val,matU);
 
 
  return 0;
}

int getDerivatives3D(double    x0, double    y0, double    z0, 
                     double *xcon, double *ycon, double *zcon,
                     double  *fun, int n ,int trans, const double *matU, double *val){
  int index;
  int m=n+1;
  int i,j,k;
  double l0x[n+1],l0y[n+1],l0z[n+1];
  double l1x[n+1],l1y[n+1],l1z[n+1];
  double l2x[n+1],l2y[n+1],l2z[n+1];

  double f;
  double val0,valx,valy,valz;
  double valxx,valyy,valzz;
  double valxy,valxz,valyz;
  double xxmu[n+1];
  double yymu[n+1];
  double zzmu[n+1];
  double deltax[n+1];
  double deltay[n+1];
  double deltaz[n+1];


  getLagParameters(x0,xcon,xxmu,deltax,n);
  getLagParameters(y0,ycon,yymu,deltay,n);
  getLagParameters(z0,zcon,zzmu,deltaz,n);

  for( i = 0 ; i <= n ; i++ ){
    l0x[i] = lagrange (xxmu,deltax[i],i,n);
    l0y[i] = lagrange (yymu,deltay[i],i,n);
    l0z[i] = lagrange (zzmu,deltaz[i],i,n);

    l1x[i] = lagrange1(xxmu,deltax[i],i,n);
    l1y[i] = lagrange1(yymu,deltay[i],i,n);
    l1z[i] = lagrange1(zzmu,deltaz[i],i,n);

    l2x[i] = lagrange2(xxmu,deltax[i],i,n);
    l2y[i] = lagrange2(yymu,deltay[i],i,n);
    l2z[i] = lagrange2(zzmu,deltaz[i],i,n);
  }

  val0  = (double) 0.;
  valx  = valy  = valz  = (double) 0.;
  valxx = valxy = valxz = (double) 0.;
  valyy = valyz = valzz = (double) 0.;

  for( i = 0 ; i <=n ; i++ )
    for( j = 0 ; j <=n ; j++ )
      for( k = 0 ; k <=n ; k++ ){
        index = i*m*m + j*m + k;
        f = fun[index];
        
        val0  += f*l0x[i]*l0y[j]*l0z[k];

        valx  += f*l1x[i]*l0y[j]*l0z[k];
        valy  += f*l0x[i]*l1y[j]*l0z[k];
        valz  += f*l0x[i]*l0y[j]*l1z[k];

        valxx += f*l2x[i]*l0y[j]*l0z[k];
        valyy += f*l0x[i]*l2y[j]*l0z[k];
        valzz += f*l0x[i]*l0y[j]*l2z[k];

        valxy += f*l1x[i]*l1y[j]*l0z[k];
        valxz += f*l1x[i]*l0y[j]*l1z[k];
        valyz += f*l0x[i]*l1y[j]*l1z[k];
      }

  val[0] = val0;
  val[1] = valx;
  val[2] = valy;
  val[3] = valz;
  val[4] = valxx;
  val[5] = valyy;
  val[6] = valzz;
  val[7] = valxy;
  val[8] = valxz;
  val[9] = valyz;

  if(trans != YES )
    trans02(val,matU);
 
  return 0;
}
/*************************************************************************/
int gradient3DLog(double    x0, double    y0, double    z0, 
                  double *xcon, double *ycon, double *zcon,
                  double  *fun, int n ,int trans, const double *matU, double min, double *val){
  int index;
  int m = n+1;
  int i,j,k;
  double l0x[n+1],l0y[n+1],l0z[n+1];
  double l1x[n+1],l1y[n+1],l1z[n+1];

  double f,val0;
  double efx;
  double valx,valy,valz;
  double xxmu[n+1];
  double yymu[n+1];
  double zzmu[n+1];
  double deltax[n+1];
  double deltay[n+1];
  double deltaz[n+1];


  getLagParameters(x0,xcon,xxmu,deltax,n);
  getLagParameters(y0,ycon,yymu,deltay,n);
  getLagParameters(z0,zcon,zzmu,deltaz,n);

  for( i = 0 ; i <= n ; i++ ){
    l0x[i] = lagrange (xxmu,deltax[i],i,n);
    l0y[i] = lagrange (yymu,deltay[i],i,n);
    l0z[i] = lagrange (zzmu,deltaz[i],i,n);

    l1x[i] = lagrange1(xxmu,deltax[i],i,n);
    l1y[i] = lagrange1(yymu,deltay[i],i,n);
    l1z[i] = lagrange1(zzmu,deltaz[i],i,n);
  }

  val0  = (double) 0.;
  valx  = valy  = valz  = (double) 0.;

  for( i = 0 ; i <=n ; i++ )
    for( j = 0 ; j <=n ; j++ )
      for( k = 0 ; k <=n ; k++ ){
        index = i*m*m + j*m + k;
        f = fun[index];
        
        val0  += f*l0x[i]*l0y[j]*l0z[k];

        valx  += f*l1x[i]*l0y[j]*l0z[k];
        valy  += f*l0x[i]*l1y[j]*l0z[k];
        valz  += f*l0x[i]*l0y[j]*l1z[k];

      }

  efx = exp(val0);

  val[0] = efx + min - DELTA;
  val[1] = efx*valx;
  val[2] = efx*valy;
  val[3] = efx*valz;

  if(trans != YES )
    trans01(val,matU);
 
 
  return 0;
}

int getDerivatives3DLog( double    x0, double    y0, double    z0, 
                         double *xcon, double *ycon, double *zcon,
                         double  *fun, int n ,int trans, const double *matU, double min, double *val){
  int index;
  int m=n+1;
  int i,j,k;
  double l0x[n+1],l0y[n+1],l0z[n+1];
  double l1x[n+1],l1y[n+1],l1z[n+1];
  double l2x[n+1],l2y[n+1],l2z[n+1];

  double f,efx;
  double val0,valx,valy,valz;
  double valxx,valyy,valzz;
  double valxy,valxz,valyz;
  double xxmu[n+1];
  double yymu[n+1];
  double zzmu[n+1];
  double deltax[n+1];
  double deltay[n+1];
  double deltaz[n+1];


  getLagParameters(x0,xcon,xxmu,deltax,n);
  getLagParameters(y0,ycon,yymu,deltay,n);
  getLagParameters(z0,zcon,zzmu,deltaz,n);

  for( i = 0 ; i <= n ; i++ ){
    l0x[i] = lagrange (xxmu,deltax[i],i,n);
    l0y[i] = lagrange (yymu,deltay[i],i,n);
    l0z[i] = lagrange (zzmu,deltaz[i],i,n);

    l1x[i] = lagrange1(xxmu,deltax[i],i,n);
    l1y[i] = lagrange1(yymu,deltay[i],i,n);
    l1z[i] = lagrange1(zzmu,deltaz[i],i,n);

    l2x[i] = lagrange2(xxmu,deltax[i],i,n);
    l2y[i] = lagrange2(yymu,deltay[i],i,n);
    l2z[i] = lagrange2(zzmu,deltaz[i],i,n);
  }

  val0  = (double) 0.;
  valx  = valy  = valz  = (double) 0.;
  valxx = valxy = valxz = (double) 0.;
  valyy = valyz = valzz = (double) 0.;

  for( i = 0 ; i <=n ; i++ )
    for( j = 0 ; j <=n ; j++ )
      for( k = 0 ; k <=n ; k++ ){
        index = i*m*m + j*m + k;
        f = fun[index];
        
        val0  += f*l0x[i]*l0y[j]*l0z[k];

        valx  += f*l1x[i]*l0y[j]*l0z[k];
        valy  += f*l0x[i]*l1y[j]*l0z[k];
        valz  += f*l0x[i]*l0y[j]*l1z[k];

        valxx += f*l2x[i]*l0y[j]*l0z[k];
        valyy += f*l0x[i]*l2y[j]*l0z[k];
        valzz += f*l0x[i]*l0y[j]*l2z[k];

        valxy += f*l1x[i]*l1y[j]*l0z[k];
        valxz += f*l1x[i]*l0y[j]*l1z[k];
        valyz += f*l0x[i]*l1y[j]*l1z[k];
      }

  efx = exp(val0);

  val[0] = efx + min - DELTA;
  val[1] = efx*valx;
  val[2] = efx*valy;
  val[3] = efx*valz;

  val[4] = efx*(valx*valx + valxx);
  val[5] = efx*(valy*valy + valyy);
  val[6] = efx*(valz*valz + valzz);

  val[7] = efx*(valx*valy + valxy);
  val[8] = efx*(valx*valz + valxz);
  val[9] = efx*(valy*valz + valyz);

  if(trans != YES )
    trans02(val,matU);
 
  return 0;
}
