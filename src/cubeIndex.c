#include "cubeIndex.h"

int getLeftRight( int pol ,int *min, int *max){
  int tmin,tmax;
  int dat = pol + 1;

  switch (dat){
    case  2: tmin =   0; tmax =  1; break;
    case  3: tmin =  -1; tmax =  1; break;

    case  4: tmin =  -1; tmax =  2; break;
    case  5: tmin =  -2; tmax =  2; break;

    case  6: tmin =  -2; tmax =  3; break;
    case  7: tmin =  -3; tmax =  3; break;

    case  8: tmin =  -3; tmax =  4; break;
    case  9: tmin =  -4; tmax =  4; break;

    case 10: tmin =  -4; tmax =  5; break;
    case 11: tmin =  -5; tmax =  5; break;

    case 12: tmin =  -5; tmax =  6; break;
    case 13: tmin =  -6; tmax =  6; break;

    case 14: tmin =  -6; tmax =  7; break;
    case 15: tmin =  -7; tmax =  7; break;

    case 16: tmin =  -7; tmax =  8; break;
    case 17: tmin =  -8; tmax =  8; break;

    case 18: tmin =  -8; tmax =  9; break;
    case 19: tmin =  -9; tmax =  9; break;

    case 20: tmin = - 9; tmax = 10; break;
    case 21: tmin = -10; tmax = 10; break;
  }
  
  (*min) = tmin;
  (*max) = tmax;
  
  return dat*dat*dat;
}

int checkIndex(int* index,int *vecn, dataRun param){
  int i;
  int id1,id2;

  /*int p,q,r;
  p = index[0];
  q = index[1];
  r = index[2];*/


  for(i = 0; i < 3; i++ ){
    id1 = index[i] + param.izq; 
    id2 = index[i] + param.der;

    if( id1 < 0 )
      index[i] = index[i] - id1;
    if( id2 >= vecn[i] )
      index[i] = index[i] -( id2 - (vecn[i] - 1));

  }
  //printf(" CHECK INDEX: %3d -> %3d | %3d -> %3d | %3d -> %3d\n",
  //               p,index[0],q,index[1],r,index[2]);

}


int getIndex( int ptq, double q, double q0, double hq ){
  int seguir;
  int i,index;
  double qi,qj;

  i = 0;
  qi = q0;
  seguir = YES;

  //while( i < ptq-1 && seguir == YES ){
  while( i < ptq && seguir == YES ){
    qj = qi + hq;
    if( q > qi && q < qj ){
      index = i;
      seguir = !YES;
    }
    qi = qj;
    i++;
  }

  return index;

}

int getPeriodicIndex(int p, int np ){
  
  np = np -1;
  p = p%np; // con esto meto p al intervalo
            //  [0 , np-1)
            //  funciona muy bien para cubos de 
            //  crystal 09

  //if (p >= 0 && p < np - 1 )
  if (p >= 0 && p <= np - 1 ) //original 
    return p;

  if( p < 0)
    return ( p + np );

  //if( p >= np - 1)
  if( p > np - 1) // original
    return ( p - np - 1);
}


int getIndex3D( int *pts, double *r, double *min, double *h, int *index){
  int i;
  for( i = 0; i < 3 ; i++ )
    index[i] = getIndex(pts[i],r[i],min[i],h[i]);

  return 0;
}

int loadLocalField(int *index,dataCube cube, dataRun param, double *xx, double *yy, double *zz, double *f){
 
  int n1 = cube.pts[1]*cube.pts[2];
  int n2 = cube.pts[2];
  int i,j,k;
  int p,q,r,mu;
  int p2,q2,r2;

  if( param.pbc == NOT) {
  // In this part  the field is load for evaluating the lagrange interpolators
  // without boundary conditions periodic
    checkIndex(index,cube.pts,param);
  
    i = index[0];
    j = index[1];
    k = index[2];

    mu = 0;
    for( p = i + param.izq ; p <= i + param.der ; p++ ){
      xx[mu] = cube.min[0] + p*cube.hvec[0];
      mu++;
    }
  
    mu = 0;
    for( p = j + param.izq ; p <= j + param.der ; p++ ){
      yy[mu] = cube.min[1] + p*cube.hvec[1];
      mu++;
    }
  
    mu = 0;
    for( p = k + param.izq ; p <= k + param.der ; p++ ){
      zz[mu] = cube.min[2] + p*cube.hvec[2];
      mu++;
    }

    mu = 0;
    for( p = i + param.izq ; p <= i + param.der ; p++ ){
      for( q = j + param.izq ; q <= j + param.der ; q++ ){
        for( r = k + param.izq ; r <= k + param.der ; r++ ){
          f[mu] = cube.field[ p*n1 + q*n2 + r];
          mu++;

        } 
      }
    }

  }else{
  // In this part  the field is load for evaluating the lagrange interpolators
  // with boundary conditions periodic
    i = index[0];
    j = index[1];
    k = index[2];

    mu = 0;
    for( p = i + param.izq ; p <= i + param.der ; p++ ){
      xx[mu] = cube.min[0] + p*cube.hvec[0];
      mu++;
    }
  
    mu = 0;
    for( p = j + param.izq ; p <= j + param.der ; p++ ){
      yy[mu] = cube.min[1] + p*cube.hvec[1];
      mu++;
    }
  
    mu = 0;
    for( p = k + param.izq ; p <= k + param.der ; p++ ){
      zz[mu] = cube.min[2] + p*cube.hvec[2];
      mu++;
    }

    mu = 0;
    for( p = i + param.izq ; p <= i + param.der ; p++ ){
      p2 = getPeriodicIndex(p,cube.pts[0]);

      for( q = j + param.izq ; q <= j + param.der ; q++ ){
        q2 = getPeriodicIndex(q,cube.pts[1]);

        for( r = k + param.izq ; r <= k + param.der ; r++ ){
          r2 = getPeriodicIndex(r,cube.pts[2]);

          f[mu] = cube.field[ p2*n1 + q2*n2 + r2];
          mu++;

        }
      }
    }


  }
 
  return 0;
}
