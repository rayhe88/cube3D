
#include "pruebaCoef.h"

int getCoeff(int poly,double *hvec, double x0,double y0, double z0,double *f, double *c){

  int size;
  int size2;
  int i,j,k,p;
  int l,m,n,mu;
  double x,y,z;
  double min,max;
  double facx,facy,facz;
  double *matX,*vecC;
  int *vecI;

  p = poly + 1;
  int nrhs = 1;

  size = pow(p,3);
  size2 = size*size;

  vecC = (double*) malloc(size*sizeof(double));
  vecI = (int*) malloc(size*sizeof(int));
  matX = (double*) malloc(size2*sizeof(double));

  for(i=0;i<size;i++)
    vecC[i] = f[i];

  switch(poly){
    case 1: min = -1;
            max =  0; break;
    case 2: min = -1;
            max =  1; break;
    case 3: min = -2;
            max =  1; break;
    case 4: min = -2;
            max =  2; break;
    case 5: min = -3;
            max =  2; break;
    case 6: min = -3;
            max =  3; break;
    case 7: min = -4;
            max =  3; break;
    case 8: min = -4;
            max =  4; break;
    case 9: min = -5;
            max =  4; break;
    case 10: min = -5;
            max =  5; break;
  }
  double hx,hy,hz;

  hx = hvec[0];
  hy = hvec[1];
  hz = hvec[2];

  mu = 0;
  for( i = min; i <= max; i++){
    x = x0 + i*hx;
    for( j = min; j <= max; j++){
      y = y0 + j*hy;
      for( k = min; k <= max; k++){
        z = z0 + k*hz;

        for(m=0; m<p ; m++){
          facx = pow(x,m);
          for(n=0; n<p ; n++){
            facy = pow(y,n);
            for(l=0; l<p ; l++){
              facz = pow(z,l);
              if( m == 0 && n == 0 && l == 0)
                matX[mu] = 1.;
              else
                matX[mu] = facx*facy*facz;

              mu++;
            }
          }
        }

      }
    }
  }
    


  int info;

  info = LAPACKE_dgesv(LAPACK_ROW_MAJOR,size,nrhs,matX,size,vecI,vecC,1);

  if( info != 0)
    printf(" Error en la diagonalizacion\n");

  for(i=0;i<size;i++){
    c[i] = vecC[i];
  }


  free(matX);
  free(vecI);
  free(vecC);

  

  return 0;

}

int loadField(int poly,int i, int j, int k,int npx, int npy, int npz, double *field, double *fun){

  int ip,jp,kp;
  int nu,idx;
  int tmpi,tmpj,tmpk;

  int min, max;

  switch(poly){
    case 1: min = 1;
            max = 0; break;
    case 2: min = 1;
            max = 1; break;
    case 3: min = 2;
            max = 1; break;
    case 4: min = 2;
            max = 2; break;
    case 5: min = 3;
            max = 2; break;
    case 6: min = 3;
            max = 3; break;
    case 7: min = 4;
            max = 3; break;
    case 8: min = 4;
            max = 4; break;
    case 9: min = 5;
            max = 4; break;
    case 10: min = 5;
             max = 5; break;
  }

  nu = 0;
  for(ip = i - min; ip <= i + max; ip++)
    for(jp = j - min; jp <= j + max; jp++)
      for(kp = k - min; kp <= k + max; kp++){

        tmpi = ip;
        tmpj = jp;
        tmpk = kp;
        // De esta forma meto  las condiciones periodicas del sistema
	     if(ip < 0) tmpi += npx;  
	     if(jp < 0) tmpj += npy;
	     if(kp < 0) tmpk += npz;

	     if(ip >= npx ) tmpi -= npx;
	     if(jp >= npy ) tmpj -= npy;
	     if(kp >= npz ) tmpk -= npz;
        // finalizo las condiciones periodicas 

        idx = tmpi*npy*npz + tmpj*npz + tmpk; 
        fun[nu] = field[idx];
        nu++;
      }

  

  return 0;
}

