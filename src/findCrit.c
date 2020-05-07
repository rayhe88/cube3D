/**
 * @file  findCrit.c
 * @brief 
 * @author Raymundo Hernández-Esparza.
 * @date   August 2018.
 */
#include "file.h"
#include "array.h"
#include "utils.h"
#include "fields.h"
#include "tableP.h"
#include "version.h"
#include "findCrit.h"
#include "cubeIndex.h"
#include "mathTools.h"
#include "numBondPath.h"

#include "transU.h"

#include <omp.h>

void printPoss(int npc, double* coorIn, const double *matU){
    int i;

    double q[3];
    double r[3];

    FILE *f;

    f = fopen("posibles.xyz","w+");

    fprintf(f," %5d\n",npc);
    fprintf(f,"    \n");
    for(i=0;i<npc;i++){
        q[0] = coorIn[3*i];
        q[1] = coorIn[3*i+1];
        q[2] = coorIn[3*i+2];
        getRiU(q,matU,r);

        fprintf(f,"  PC  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);


    }

    fclose(f);


}

int createArrayCritP( int n, dataCritP **ptr, const char *mess){
  
  if ( n < 1 ){
    printf(" There is an error in the size for [%s]\n",mess);
    exit(EXIT_FAILURE);
  }

  *ptr = (dataCritP *) malloc( n*sizeof(dataCritP));
  if( *ptr == NULL ){
    printf(" Failed to allocate memory : [%s]\n",mess);
    exit(EXIT_FAILURE);
  }

  return 0;
}

int critPoints ( dataCube cube, dataRun param,const double *matU, double min0, char *name ){

  int    i,j,k;
  int    ncx = cube.pts[0] - 1;
  int    ncy = cube.pts[1] - 1;
  int    ncz = cube.pts[2] - 1;
  double  hx = cube.hvec[0];
  double  hy = cube.hvec[1];
  double  hz = cube.hvec[2];
  double x,y,z;
  double rout[8][3],rin[3];
  double fun;
  int    *ctrl  = NULL;
  int     nct = ncx*ncy*ncz;
  int      n1 = cube.pts[1]*cube.pts[2];
  int      n2 = cube.pts[2];
  int      mu,idx[3],npc=0;

  int nc1 = ncy*ncz;
  int nc2 = ncz;

  char tmpname[128];
  FILE *tmp;
  tmpFile(&tmp,".c3dFindC",tmpname,"w+");

  cube.max[0] = cube.min[0] + (cube.pts[0] - 1 ) * cube.hvec[0];
  cube.max[1] = cube.min[1] + (cube.pts[1] - 1 ) * cube.hvec[1];
  cube.max[2] = cube.min[2] + (cube.pts[2] - 1 ) * cube.hvec[2];


  double norm;

  /*************************************************************/
  for( i = 0; i < cube.natm; i++ ){
    rin[0] = cube.coor[3*i];
    rin[1] = cube.coor[3*i+1];
    rin[2] = cube.coor[3*i+2];
    trans00(rin,matU);

    cube.coor[3*i]   = rin[0]; 
    cube.coor[3*i+1] = rin[1]; 
    cube.coor[3*i+2] = rin[2];
  }
  rin[0] = cube.min[0];
  rin[1] = cube.min[1];
  rin[2] = cube.min[2];
  trans00(rin,matU);
  cube.min[0] = rin[0]; 
  cube.min[1] = rin[1]; 
  cube.min[2] = rin[2];

  /*
  rin[0] = cube.max[0];
  rin[1] = cube.max[1];
  rin[2] = cube.max[2];
  trans00(rin,matU);
  cube.max[0] = rin[0]; 
  cube.max[1] = rin[1]; 
  cube.max[2] = rin[2];
  */

  /*************************************************************/

  double  x0 = cube.min[0];
  double  y0 = cube.min[1];
  double  z0 = cube.min[2];

  printBanner(" Searching of critical points ",stdout);
  printf("  Cutoff for the function         : % 10.6E\n",TOLFUN);
  printf("  Tolerance in the gradient       : % 10.6E\n",TOLGRD);
  printf("  Cube overlap percentage         : % 10.2lf%c \n",PERCENT*100.,'%');
  printf("  Maximum iterations step 1       : % 10d\n",MAXITER1);
  printf("  Maximum iterations step 2       : % 10d\n",MAXITER2);
  printBar(stdout);

  createArrayInt(nct,&ctrl,"Index Array");
  for( i = 0; i < nct ; i++ )
    ctrl[i] = -1;

  // Primer discriminante Elimamos los puntos donde el campo es menor a 
  // TOLFUN0
  npc = 0;
  for( i = 0; i < ncx ; i++ ){
    for( j = 0 ; j < ncy ; j++ ){
      for( k = 0 ; k < ncz ; k++ ){
        mu = i*ncy*ncz + j*ncz + k;
        fun  = getFunInCube(i,j,k,n1,n2,min0,cube.hvec,cube.field,&norm);
        if( fun >= TOLFUN0  ){
          npc++;
          ctrl [mu] = 0;
        }
      }
    }
  }
/**************************************/

  printf("  Search space %7d cubes\n",npc);
  printf("  Corresponds to %6.2lf%c  total\n",(double) (100.*npc)/(ncx*ncy*ncz),'%');

/**************************************/

  double *xyz1=NULL;
  double xp,yp,zp;
  double limits[6];
  double xpp,ypp,zpp;

  double f[param.size];
  double xx[param.pol + 1];
  double yy[param.pol + 1];
  double zz[param.pol + 1];

  int ret[8];
  int nu = 0;

#pragma omp parallel private(i,j,k,nu,mu,idx,x,y,z,rin,rout,xp,yp,zp,xpp,ypp,zpp,ret,xx,yy,zz,f,limits)\
                     shared(nc1,nc2,nct,ctrl,x0,y0,z0,hx,hy,hz,cube,matU,param,min0,tmp)
{
#pragma omp single
{
  printBar(stdout);
  printf("  Number of threads for searching : %6d\n",omp_get_num_threads());
  printBar(stdout);
}
#pragma omp barrier

#pragma omp for schedule (dynamic)
  for( mu = 0; mu < nct; mu++){
    i = mu/nc1;
    j = (mu - i*nc1)/nc2;
    k = mu%nc2;

    if( ctrl[mu] == 0 ) {
       idx[0] = i;
       idx[1] = j;
       idx[2] = k;
       x = x0 + (double) i*hx;
       y = y0 + (double) j*hy;
       z = z0 + (double) k*hz;

       rin[0] = x + 0.5*hx;
       rin[1] = y + 0.5*hy;
       rin[2] = z + 0.5*hz;

       loadLocalField(idx,cube,param,xx,yy,zz,f);

       getLimits(rin,cube.hvec,limits);
       
       xp = x + 0.25*hx;
       yp = y + 0.25*hy;
       zp = z + 0.25*hz;

       xpp = x + 0.75*hx;
       ypp = y + 0.75*hy;
       zpp = z + 0.75*hz;
       

       rin[0] = xp; rin[1] = yp; rin[2] = zp;
       ret[0] = rejectCube(rin,min0,cube,matU,param,limits,xx,yy,zz,f,rout[0]);

       rin[0] = xp; rin[1] = yp; rin[2] = zpp;
       ret[1] = rejectCube(rin,min0,cube,matU,param,limits,xx,yy,zz,f,rout[1]);

       rin[0] = xp; rin[1] = ypp; rin[2] = zp;
       ret[2] = rejectCube(rin,min0,cube,matU,param,limits,xx,yy,zz,f,rout[2]);

       rin[0] = xp; rin[1] = ypp; rin[2] = zpp;
       ret[3] = rejectCube(rin,min0,cube,matU,param,limits,xx,yy,zz,f,rout[3]);

       rin[0] = xpp; rin[1] = yp; rin[2] = zp;
       ret[4] = rejectCube(rin,min0,cube,matU,param,limits,xx,yy,zz,f,rout[4]);

       rin[0] = xpp; rin[1] = yp; rin[2] = zpp;
       ret[5] = rejectCube(rin,min0,cube,matU,param,limits,xx,yy,zz,f,rout[5]);

       rin[0] = xpp; rin[1] = ypp; rin[2] = zp;
       ret[6] = rejectCube(rin,min0,cube,matU,param,limits,xx,yy,zz,f,rout[6]);

       rin[0] = xpp; rin[1] = ypp; rin[2] = zpp;
       ret[7] = rejectCube(rin,min0,cube,matU,param,limits,xx,yy,zz,f,rout[7]);
       
#pragma omp critical
{
       for( nu = 0; nu < 8 ; nu++ ){
         if( ret[nu] == 1 )
             fprintf(tmp," % 20.12lf % 20.12lf % 20.12lf \n",rout[nu][0],rout[nu][1],rout[nu][2]); 
       }

} // fin omp critical


    } // fin if
           
  }//fin for

}//fin omp
  free(ctrl);

  char c;
  rewind(tmp);
  nu = 0;
  while( (c = fgetc(tmp)) != EOF){
    if( c == '\n')
      nu++;
  }
  rewind(tmp);
  npc = nu;
  createArrayDou(3*npc,&xyz1,"Coordenadas 1");
  for(nu=0;nu<npc;nu++)
    fscanf(tmp,"%lf %lf %lf\n",&xyz1[3*nu],&xyz1[3*nu+1],&xyz1[3*nu+2]);

  fclose(tmp);
  remove(tmpname);
  
  printf("  Possible critical points        : %6d \n",npc);

  refineCrit(npc,min0,cube,param,xyz1,matU,name);
  free(xyz1);
  return 0;
}

int rejectCube(double *rin, double min0, dataCube cube, const double *matU, dataRun param,
               double *limits,double *xx, double *yy, double *zz, double *f, double *rout){

  int ret= -1;
  int iter,difcoor;
  double xi = rin[0];
  double yi = rin[1];
  double zi = rin[2];
  double fun, norm;
  double val[10];
  double x,y,z;
  double xn,yn,zn;
  double alphax = ALPHA(cube.hvec[0]);
  double alphay = ALPHA(cube.hvec[1]);
  double alphaz = ALPHA(cube.hvec[2]);

  
  gradient3DLog(xi,yi,zi,xx,yy,zz,f,param.pol,param.orth,matU,min0,val);
  //gradient3D(xi,yi,zi,xx,yy,zz,f,param.pol,param.orth,matU,val);

  fun  = fabs(val[0]);
  norm = getGrd(val);
  rout[0] = 1.E5;
  rout[1] = 1.E5;
  rout[2] = 1.E5;

  if( fun < TOLFUN && norm > TOLNRM ){
    ret = -1;
  }else{
    iter = 0;
    norm = 1.E4;
    difcoor = 0;
    x = xi; y = yi; z = zi;
    
    while( (iter < MAXITER1) && (norm > TOLGRD) && (difcoor == 0) ){

      getDerivatives3DLog(x,y,z,xx,yy,zz,f,param.pol,param.orth,matU,min0,val);
     // getDerivatives3D(x,y,z,xx,yy,zz,f,param.pol,param.orth,matU,val);

      fun  = fabs(val[0]);
      if( fun < TOLFUN ){
        iter = 2*MAXITER1;
        x = 1.E5;
        y = 1.E5;
        z = 1.E5;
      }else{
        hessGrad(&xn,&yn,&zn,&norm,val);
        x += alphax*xn;
        y += alphay*yn;
        z += alphaz*zn;
        iter++;
      }
      difcoor = inCube(x,y,z,limits);
    }// end of while

      if( difcoor == 0 && norm <= TOLNRM){
        gradient3DLog(x,y,z,xx,yy,zz,f,param.pol,param.orth,matU,min0,val);
       // gradient3D(x,y,z,xx,yy,zz,f,param.pol,param.orth,matU,val);
       
        fun  = fabs(val[0]);
        norm = getGrd(val);
        ret = 1;
        rout[0] = x;
        rout[1] = y;
        rout[2] = z;
        
        if( inMacroCube(x,y,z,cube.min,cube.max) == 1 )
          ret = -1;


      }


  }//fin del if-else




  return ret;
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

  a = val[4];
  d = val[5];
  g = val[6];

  b = val[7];
  c = val[8];
  e = val[9];

  (*ngrad) = sqrt(vx*vx + vy*vy + vz*vz);

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

int getLimits(double rin[3],double *h,double *limits){


  limits[0] = rin[0] - PERCENT*h[0]*0.5;
  limits[1] = rin[0] + PERCENT*h[0]*0.5;

  limits[2] = rin[1] - PERCENT*h[1]*0.5;
  limits[3] = rin[1] + PERCENT*h[1]*0.5;

  limits[4] = rin[2] - PERCENT*h[2]*0.5;
  limits[5] = rin[2] + PERCENT*h[2]*0.5;


  return 0;
}

int inCube(double x, double y, double z, double limits[6]){

// esta función devolverá 0 si estamos dentro del cubo mas una tolerancia
// si no devolvera 1 si estamos fuera
// las distancias limite estan dadas por limits

  int flag = 0;
  
  if( limits[0] < x && limits[1] > x )
    flag += 1;

  if( limits[2] < y && limits[3] > y )
    flag += 2;

  if( limits[4] < z && limits[5] > z )
    flag += 4;
  
  if( flag == 7 )
    return 0;
  else
    return 1;
}

int inMacroCube(double x,double y,double z,double *min, double *max) {
  int flag = 0;

  if( min[0] < x && max[0] > x )
    flag += 1;

  if( min[1] < y && max[1] > y )
    flag += 2;

  if( min[2] < z && max[2] > z )
    flag += 4;

  if (flag == 7)
    return 0;
  else
    return 1;

}

int numCritical01Vec(double r[3], dataCube cube, dataRun param, const double *matU, double min0, double *val){

  int ret = numCritical01(r[0],r[1],r[2],cube, param, matU, min0, val);

  return ret;
}

int numCritical01(double x, double y, double z, dataCube cube, 
                  dataRun param, const double *matU,double min0,double *val){

  int idx[3];
  double f[param.size];
  double xx[param.pol + 1];
  double yy[param.pol + 1];
  double zz[param.pol + 1];

  idx[0] = getIndex(cube.pts[0],x,cube.min[0],cube.hvec[0]);
  idx[1] = getIndex(cube.pts[1],y,cube.min[1],cube.hvec[1]);
  idx[2] = getIndex(cube.pts[2],z,cube.min[2],cube.hvec[2]);

  loadLocalField(idx,cube,param,xx,yy,zz,f);
  gradient3DLog(x,y,z,xx,yy,zz,f,param.pol,param.orth,matU,min0,val);
  //gradient3D(x,y,z,xx,yy,zz,f,param.pol,param.orth,matU,val);


  return 0;
}

int numCritical02Vec(double r[3], dataCube cube, dataRun param, const double *matU, double min0, double *val){

  int ret = numCritical02(r[0],r[1],r[2],cube, param, matU, min0, val);

  return ret;
}

int numCritical02(double x, double y, double z, dataCube cube, 
                  dataRun param, const double *matU,double min0,double *val){

  int idx[3];
  double f[param.size];
  double xx[param.pol + 1];
  double yy[param.pol + 1];
  double zz[param.pol + 1];

  idx[0] = getIndex(cube.pts[0],x,cube.min[0],cube.hvec[0]);
  idx[1] = getIndex(cube.pts[1],y,cube.min[1],cube.hvec[1]);
  idx[2] = getIndex(cube.pts[2],z,cube.min[2],cube.hvec[2]);

  loadLocalField(idx,cube,param,xx,yy,zz,f);
  getDerivatives3DLog(x,y,z,xx,yy,zz,f,param.pol,param.orth,matU,min0,val);
  //getDerivatives3D(x,y,z,xx,yy,zz,f,param.pol,param.orth,matU,val);

  return 0;
}
int numCritical02Vec_exp(double r[3], dataCube cube, dataRun param, const double *matU, double min0, double *val){
  int i,idx[3];
  double f[param.size];
  double xx[param.pol + 1];
  double yy[param.pol + 1];
  double zz[param.pol + 1];


  idx[0] = getIndex(cube.pts[0],r[0],cube.min[0],cube.hvec[0]);
  idx[1] = getIndex(cube.pts[1],r[1],cube.min[1],cube.hvec[1]);
  idx[2] = getIndex(cube.pts[2],r[2],cube.min[2],cube.hvec[2]);

  loadLocalField(idx,cube,param,xx,yy,zz,f);

  for(i=0;i<param.size;i++){
    f[i] = exp(f[i]) + min0 - DELTA;
  }

  getDerivatives3D(r[0],r[1],r[2],xx,yy,zz,f,param.pol,param.orth,matU,val);

  return 0;
}

double getFunInCube(int i, int j, int k, int n1, int n2, double min0, double *h,double *field,double *norm){
  
  int ip = i+1;
  int jp = j+1;
  int kp = k+1;

  double den = 0;
  double dx,dy,dz;
  double val[8];
  double grad[3];
  double efx,f;

  val[0] = field[i *n1 + j *n2 + k ];
  val[1] = field[i *n1 + j *n2 + kp];
  val[2] = field[i *n1 + jp*n2 + k ];
  val[3] = field[i *n1 + jp*n2 + kp];
  
  val[4] = field[ip*n1 + j *n2 + k ];
  val[5] = field[ip*n1 + j *n2 + kp];
  val[6] = field[ip*n1 + jp*n2 + k ];
  val[7] = field[ip*n1 + jp*n2 + kp];


   f  =  val[0] + val[1] + val[2] + val[3] + val[4] + val[5] + val[6] + val[7];

  dx  = -val[0] - val[1] - val[2] - val[3] + val[4] + val[5] + val[6] + val[7];

  dy  = -val[0] - val[1] + val[2] + val[3] - val[4] - val[5] + val[6] + val[7];

  dz  = -val[0] + val[1] - val[2] + val[3] - val[4] + val[5] - val[6] + val[7];

  f  *= 0.125;
  
  efx = exp(f);;

  den = efx + min0 - DELTA;


  dx /= (4.*h[0]);
  dy /= (4.*h[1]);
  dz /= (4.*h[2]);

  grad[0] = efx*dx;
  grad[1] = efx*dy;
  grad[2] = efx*dz;


  (*norm) = sqrt( grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2]);


  return fabs(den);
}


int  refineCrit(int npc, double min0, dataCube cube, dataRun param, double *coorIn, const double *matU, char *name){
  
  int mu,nu;
  int ncrit,ncrit2,ncrit3;
  int iter;
  double *xyz= NULL;
  double *coorOut=NULL;
  double x,y,z;
  double xn,yn,zn;
  double norm;
  double val[10];


  createArrayDou(3*npc,&xyz,"Coordinates Intermediate");

  printPoss(npc,coorIn,matU); //  XXX- Delete after test

  nu=0;
  for( mu = 0; mu < npc; mu++){
    x = coorIn[3*mu];
    y = coorIn[3*mu+1];
    z = coorIn[3*mu+2];

    
    iter=0; norm = 1.E5;

    while( iter < MAXITER2 && norm > TOLGRD ){
      numCritical02(x,y,z,cube,param,matU,min0,val);
      hessGrad(&xn,&yn,&zn,&norm,val);
      x += xn;
      y += yn;
      z += zn;
      iter++;
    } // end of WHILE

    if( iter == MAXITER2 ){
        while( iter < 2 * MAXITER2 && norm > TOLGRD ){
            x = (x + coorIn[3*mu]  )/2.;
            y = (y + coorIn[3*mu+1])/2.;
            z = (z + coorIn[3*mu+2])/2.;
            numCritical02(x,y,z,cube,param,matU,min0,val);
            hessGrad(&xn,&yn,&zn,&norm,val);
            x += 0.01*xn;
            y += 0.01*yn;
            z += 0.01*zn;
            iter++;
        } // end of WHILE
    }

    if( norm < 1.E-3 && inMacroCube(x,y,z,cube.min,cube.max) == 0 ){
    //if( norm < 1.E-4 && inMacroCube(x,y,z,cube.min,cube.max) == 0 ){ original
      xyz[ 3*nu    ] = x;
      xyz[ 3*nu+ 1 ] = y;
      xyz[ 3*nu+ 2 ] = z;
      nu++;
    }
  }// end of FOR 

  createArrayDou(3*npc,&coorOut,"Coordinates Finale");

  ncrit = delRepCoor(nu,xyz,coorOut);

  free(xyz);

  printf("  Critical points (delRepCoor)    : %6d \n",ncrit);
  ncrit2 = delCoorAtomic( ncrit, coorOut, matU,cube);
  printf("  Critical points (delCoorAtomic) : %6d\n",ncrit2);
  ncrit3 = delCoorPseudo( ncrit2, coorOut,min0,matU,cube,param);
  printf("  Critical points (delCoorPseudo) : %6d\n",ncrit3);

  //ncrit3 = ncrit2;

  // Caraterizamos los puntos criticos
  describeCrit(ncrit3,min0,cube,param,coorOut,matU,name);
  
  free(coorOut);

    
  return 0;
}

int describeCrit(int ncrit, double min0, dataCube cube,
    dataRun param, double *coor,const double *matU, char *name){

  int i,j,k,l,m;
  int mu,tipo;
  int poincare;
  int ran,sig;
  int ncp = 0;
  int bcp = 0;
  int rcp = 0;
  int ccp = 0;
  int xcp = 0;
  double x,y,z;
  double val[10];
  double matH[6];
  double eval[3];
  double fun;

  dataCritP *pointC = NULL; 

  dataCritP *bondCrit = NULL;
  dataCritP *ringCrit = NULL;
  dataCritP *cageCrit = NULL;
  dataCritP *nnucCrit = NULL;
  dataCritP *degeCrit = NULL;

  createArrayCritP(ncrit, &pointC,"Critical 1");
   
  for( mu = 0; mu < ncrit; mu++ ){
    x = coor[3*mu];
    y = coor[3*mu+1];
    z = coor[3*mu+2];

    numCritical02(x,y,z,cube,param,matU,min0,val);
    matH[0] = val[4]; matH[1] = val[7]; matH[2] = val[8];
                      matH[3] = val[5]; matH[4] = val[9];
                                        matH[5] = val[6];

    //JacobiNxN (matH,eval,evec);
    valoresPropios3x3_v0(matH,eval);
    getSignature( eval, &ran, &sig);
    tipo = getTypeCrit(ran,sig);

    pointC[mu].x   = x;
    pointC[mu].y   = y;
    pointC[mu].z   = z;

    pointC[mu].ran = ran;
    pointC[mu].sig = sig;
    pointC[mu].fun = val[0];

    pointC[mu].typ = tipo;

    switch( tipo ){
      case NCP : ncp++; break;
      case BCP : bcp++; break;
      case RCP : rcp++; break;
      case CCP : ccp++; break;
      default  : xcp++; 
    }

  } // end of for


  int *bonding;
  int *cells;
  if ( bcp != 0 ) {
    createArrayCritP(bcp, &bondCrit,"Critical BCP");
    createArrayInt (2*bcp,&bonding,"Bonding array");
    createArrayInt (2*bcp,&cells,"cells array");
  }
  if ( rcp != 0 ) createArrayCritP(rcp, &ringCrit,"Critical RCP");
  if ( ccp != 0 ) createArrayCritP(ccp, &cageCrit,"Critical CCP");
  if ( ncp != 0 ) createArrayCritP(ncp, &nnucCrit,"Critical NCP");
  if ( xcp != 0 ) createArrayCritP(xcp, &degeCrit,"Critical XCP");

  i = j = k = l = m = 0;
  for( mu = 0; mu < ncrit; mu++ ){
    x = pointC[mu].x;
    y = pointC[mu].y;
    z = pointC[mu].z;
    ran = pointC[mu].ran;
    sig = pointC[mu].sig;
    fun = pointC[mu].fun;
    tipo = pointC[mu].typ;

    switch( tipo ){
      case NCP :  nnucCrit[i].x   =   x ;  nnucCrit[i].y   =   y;  nnucCrit[i].z   =   z;
                  nnucCrit[i].ran = ran ;  nnucCrit[i].sig = sig;  nnucCrit[i].fun = fun;
                  nnucCrit[i].typ = tipo;  i++;
                  break;
      case BCP :  bondCrit[j].x   =   x ;  bondCrit[j].y   =   y;  bondCrit[j].z   =   z;
                  bondCrit[j].ran = ran ;  bondCrit[j].sig = sig;  bondCrit[j].fun = fun;
                  bondCrit[j].typ = tipo;  j++;
                  break;
      case RCP :  ringCrit[k].x   =   x ;  ringCrit[k].y   =   y;  ringCrit[k].z   =   z;
                  ringCrit[k].ran = ran ;  ringCrit[k].sig = sig;  ringCrit[k].fun = fun;
                  ringCrit[k].typ = tipo;  k++;
                  break;
      case CCP :  cageCrit[l].x   =   x ;  cageCrit[l].y   =   y;  cageCrit[l].z   =   z;
                  cageCrit[l].ran = ran ;  cageCrit[l].sig = sig;  cageCrit[l].fun = fun;
                  cageCrit[l].typ = tipo;  l++;
                  break;
      default  :  degeCrit[m].x   =   x ;  degeCrit[m].y   =   y;  degeCrit[m].z   =   z;
                  degeCrit[m].ran = ran ;  degeCrit[m].sig = sig;  degeCrit[m].fun = fun;
                  degeCrit[m].typ = tipo;  m++;
                  break;
    }

  }


  if ( bcp != 0 ) sortCritPoints( bcp, bondCrit);
  if ( rcp != 0 ) sortCritPoints( rcp, ringCrit);
  if ( ccp != 0 ) sortCritPoints( ccp, cageCrit);
  if ( ncp != 0 ) sortCritPoints( ncp, nnucCrit);
  if ( xcp != 0 ) sortCritPoints( xcp, degeCrit);

  poincare = ncp+cube.natm - bcp + rcp - ccp; 

  printBar(stdout);
  printBanner(" Critical points Info ",stdout);
  printf("    Bond Critical Points  (BCP)   : %6d\n",bcp);
  printf("    Ring Critical Points  (RCP)   : %6d\n",rcp);
  printf("    Cage Critical Points  (CCP)   : %6d\n",ccp);
  printf("    Nuclear Attractors    (NACP)  : %6d\n",cube.natm);
  printf("    No-nuclear Attractors (NNACP) : %6d\n",ncp);
  if(xcp != 0 )
    printf("    Degenerated Critical Points   : %6d\n",xcp);
  printBar(stdout);
  printf("       Poincare-Hopf relationship \n");
  printf("    (N)NACP - BCP + RCP - CCP     = %6d\n",poincare);


  FILE *out;
  char nameOut[128];
  char symb[4];
  double ri[3];

  sprintf(nameOut,"%sCritP.xyz",name);
  openFile(&out,nameOut,"w+");
 
  /**** Imprimimos el xyz ****/
  /**** Geometria del sistema ****/
  fprintf(out,"  %5d\n",cube.natm + ncrit);
  fprintf(out,"  File created by Cube3D-%s for critical points\n",VERSION);
  for(i=0; i < cube.natm; i++ ){
    getAtomicSymbol(cube.zatm[i],4,symb);

    ri[0] = cube.coor[3*i  ]*B2A;
    ri[1] = cube.coor[3*i+1]*B2A;
    ri[2] = cube.coor[3*i+2]*B2A;
    itrans00(ri,matU);
    fprintf(out,"  %-8s % 10.6lf % 10.6lf % 10.6lf\n",symb,
                ri[0], ri[1], ri[2]);
  }
  /***** NNACP *****/
  for(i=0; i < ncp ; i++ ){
    ri[0] = nnucCrit[i].x*B2A;
    ri[1] = nnucCrit[i].y*B2A;
    ri[2] = nnucCrit[i].z*B2A;
    itrans00(ri,matU);
    fprintf(out,"  %-8s % 10.6lf % 10.6lf % 10.6lf\n","NNACP",
                ri[0], ri[1], ri[2]);
  }
  /***** BCP *****/
  for(i=0; i < bcp ; i++ ){
    ri[0] = bondCrit[i].x*B2A;
    ri[1] = bondCrit[i].y*B2A;
    ri[2] = bondCrit[i].z*B2A;
    itrans00(ri,matU);
    fprintf(out,"  BCP%04d % 10.6lf % 10.6lf % 10.6lf\n",
                i+1,ri[0], ri[1], ri[2]);
  }
  /***** RCP *****/
  for(i=0; i < rcp ; i++ ){
    ri[0] = ringCrit[i].x*B2A;
    ri[1] = ringCrit[i].y*B2A;
    ri[2] = ringCrit[i].z*B2A;
    itrans00(ri,matU);
    fprintf(out,"  %-8s % 10.6lf % 10.6lf % 10.6lf\n","RCP",
                ri[0], ri[1], ri[2]);
  }
  /***** CCP *****/
  for(i=0; i < ccp ; i++ ){
    ri[0] = cageCrit[i].x*B2A;
    ri[1] = cageCrit[i].y*B2A;
    ri[2] = cageCrit[i].z*B2A;
    itrans00(ri,matU);
    fprintf(out,"  %-8s % 10.6lf % 10.6lf % 10.6lf\n","CCP",
                ri[0], ri[1], ri[2]);
  }
  /***** XCP *****/
  for(i=0; i < xcp ; i++ ){
    ri[0] = degeCrit[i].x*B2A;
    ri[1] = degeCrit[i].y*B2A;
    ri[2] = degeCrit[i].z*B2A;
    itrans00(ri,matU);
    fprintf(out,"  %-8s % 10.6lf % 10.6lf % 10.6lf\n","XCP",
                ri[0], ri[1], ri[2]);
  }


 bondPath(bcp,bondCrit,ncp,nnucCrit,bonding,cells,cube,param,min0,matU,name);

  printBar(stdout);
  printf("  File %s was generated\n",nameOut);
  logFile (bcp,rcp,ccp,ncp,bondCrit,ringCrit,cageCrit,nnucCrit,cube,param, min0,
       bonding,cells,matU,name);
  
  logFileCSV(bcp,bondCrit,cube,param,min0,matU,name,"BCP");
  logFileCSV(rcp,ringCrit,cube,param,min0,matU,name,"RCP");
  logFileCSV(ccp,cageCrit,cube,param,min0,matU,name,"CCP");
  logFileCSV(xcp,degeCrit,cube,param,min0,matU,name,"XCP");

  //axesCrit (bcp,rcp,ccp,ncp,bondCrit,ringCrit,cageCrit,nnucCrit,cube,param, min0,
  //     matU,name);

  if ( bcp != 0 ) {
    free(bondCrit);
    free(bonding);
  }
  if ( rcp != 0 ) free(ringCrit);
  if ( ccp != 0 ) free(cageCrit);
  if ( ncp != 0 ) free(nnucCrit);
  if ( xcp != 0 ) free(degeCrit);

  free(pointC);
  return 0;
}

int delRepCoor(int n, double *xyzIn, double *xyzOut){

  //keys for control
  // 'N' is similar to NULL, initialize value
  // 'U' if the value is unique and not delete
  // 'D' if the value will be deleted
  int i,j;
  char control[n];
  double ri[3],rj[3];
  double dist;
  for( i = 0; i < n; i++ )
    control[i] = 'N';

  control[0] = 'U';

  for( i = 0; i < n-1; i++ ){
    ri[0] = xyzIn[3*i];
    ri[1] = xyzIn[3*i + 1];
    ri[2] = xyzIn[3*i + 2];

    for( j = i+1; j < n; j++ ){
      if( control[j] != 'D' ){
        rj[0] = xyzIn[3*j];
        rj[1] = xyzIn[3*j + 1];
        rj[2] = xyzIn[3*j + 2];

        dist = distance(ri,rj);
        if( dist <= TOLDIST )
          control[j] = 'D';
        else
          control[j] = 'U';
      }
    }
    
  }

  int npc2=0;
  for( i = 0; i < n; i++ )
    if( control[i] == 'U' )
      npc2++;


  j=0;
  for( i = 0; i < n; i++ )
    if( control[i] == 'U' ){
      xyzOut[3*j  ] = xyzIn[3*i  ];
      xyzOut[3*j+1] = xyzIn[3*i+1];
      xyzOut[3*j+2] = xyzIn[3*i+2];
      j++;
    }


  return npc2;
}

void getSignature( double *eval, int *ran, int *sig){
  int i;
  int rint=0,sint=0;
  double tmp;
  
  for( i = 0; i < 3; i++ ){
    tmp = eval[i];
    if ( tmp != 0. ){
      rint++;
      if( tmp < 0. ) sint--;
      if( tmp > 0. ) sint++;
    }
  }

  (*ran) = rint;
  (*sig) = sint;

}

int getTypeCrit( int ran, int sig){
  
  int ret = 0;
  
  if ( ran != 3 ) {
    ret = XCP;
  } else {
    switch (sig) {
      case -3: ret = NCP; break;
      case -1: ret = BCP; break;
      case  1: ret = RCP; break;
      case  3: ret = CCP; break;
    }
  }
  
  return ret;
}

/**
 * This function sort the Critical points by the value's function
 */
int sortCritPoints( int n, dataCritP *cp){

  int i,j;
  int itmp;
  double ftmp;

  for ( i = 0; i < n-1; i++ ){
    for ( j = i + 1; j < n; j++ ){
      if ( cp[i].fun < cp[j].fun ){

        itmp = cp[i].ran;
        cp[i].ran = cp[j].ran;
        cp[j].ran = itmp;

        itmp = cp[i].sig;
        cp[i].sig = cp[j].sig;
        cp[j].sig = itmp;

        itmp = cp[i].typ;
        cp[i].typ = cp[j].typ;
        cp[j].typ = itmp;

        ftmp = cp[i].fun;
        cp[i].fun = cp[j].fun;
        cp[j].fun = ftmp;

        ftmp = cp[i].x;
        cp[i].x = cp[j].x;
        cp[j].x = ftmp;

        ftmp = cp[i].y;
        cp[i].y = cp[j].y;
        cp[j].y = ftmp;

        ftmp = cp[i].z;
        cp[i].z = cp[j].z;
        cp[j].z = ftmp;


      }
    
    }
  }



  return 0;
}

int delCoorAtomic(int n,double *coor, const double *matU, dataCube cube){

  //keys for control
  // 'N' is similar to NULL, initialize value
  // 'U' if the value is unique and not delete
  // 'D' if the value will be deleted
  int i,j;
  char control[n];
  double ri[3],rj[3];
  double tmp[3*n];
  double dist;

  for( i = 0; i < n; i++ )
    control[i] = 'N';

  for( i = 0; i < 3*n; i++ ){
    tmp[i] = coor[i];
    coor[i] = 0.;
  }

  for( i = 0; i < cube.natm; i++ ){
    ri[0] = cube.coor[3*i];
    ri[1] = cube.coor[3*i + 1];
    ri[2] = cube.coor[3*i + 2];

    for( j = 0; j < n; j++ ){
      if( control[j] != 'D' ){
        rj[0] = tmp[3*j];
        rj[1] = tmp[3*j + 1];
        rj[2] = tmp[3*j + 2];

        dist = distance(ri,rj);

        if( dist <= TOLDIST2 )
          control[j] = 'D';
        else
          control[j] = 'U';
      }
    }
    
  }

  int npc2=0;
  for( i = 0; i < n; i++ )
    if( control[i] == 'U' )
      npc2++;


  j=0;
  for( i = 0; i < n; i++ )
    if( control[i] == 'U' ){
      coor[3*j  ] = tmp[3*i  ];
      coor[3*j+1] = tmp[3*i+1];
      coor[3*j+2] = tmp[3*i+2];
      j++;
    }


  return npc2;
}

int delCoorPseudo(int n,double *coor1, double min0,const double *matU,
                  dataCube cube,dataRun param){

  //keys for control
  // 'N' is similar to NULL, initialize value
  // 'U' if the value is unique and not delete
  // 'D' if the value will be deleted
  int i,j;
  char control[n];
  double ri[3],rj[3];
  double tmp[3*n];
  double dist;
  double val[10],fun,lap;

  for( i = 0; i < n; i++ )
    control[i] = 'U';

  for( i = 0; i < 3*n; i++ ){
    tmp[i] = coor1[i];
    coor1[i] = 0.;
  }

  for( i = 0; i < cube.natm; i++){
    ri[0] = cube.coor[3*i];
    ri[1] = cube.coor[3*i + 1];
    ri[2] = cube.coor[3*i + 2];

    numCritical02(ri[0],ri[1],ri[2],cube,param,matU,min0,val);

    fun = fabs(val[0]);
    lap = getLap(val);

    if ( fun < 0.5 && lap > 0. ){ // If these condition is true, I have  pseudopotential 
                                  // in the nuclei i;
      for( j = 0; j < n; j++ ){
        if( control[j] != 'D' ){
          rj[0] = tmp[3*j];
          rj[1] = tmp[3*j + 1];
          rj[2] = tmp[3*j + 2];

          dist = distance(ri,rj);

          if( dist <= TOLDIST3 )
            control[j] = 'D';
          else
            control[j] = 'U';
        }
        
      }
    }
  }



  int npc2=0;
  for( i = 0; i < n; i++ )
    if( control[i] == 'U' )
      npc2++;


  j=0;
  for( i = 0; i < n; i++ )
    if( control[i] == 'U' ){
      coor1[3*j  ] = tmp[3*i  ];
      coor1[3*j+1] = tmp[3*i+1];
      coor1[3*j+2] = tmp[3*i+2];
      j++;
    }


  return npc2;
}

