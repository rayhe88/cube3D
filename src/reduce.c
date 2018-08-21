#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>

#define ARRAYSIZE 1000

void   printXYZ(int n, double *array, const char*,const char *name);
void   genArray(int n, double,double,double *array);
void   deleteRepeated(int n,double* xyzInp, double *xyzOut, int *m);
double randomGen(double min, double max);


int main(int argc, char *argv[]){

  srand48((unsigned) time(NULL));

  int i,n,m;
  double min,max;

  n = atoi(argv[1]);
  min = (double) atof(argv[2]);
  max = (double) atof(argv[3]);


  double xyzInp[9*ARRAYSIZE];
  double xyzOut[9*ARRAYSIZE];

  if( n > ARRAYSIZE ){
     printf(" El tamaño es muy grande para esta prueba\n");
     exit(EXIT_FAILURE);
  }
  printf(" FLT_EPSILON % 20.16E\n",FLT_EPSILON);
  printf(" DBL_EPSILON % 20.16E\n",DBL_EPSILON);
  printf(" DBL_EPSILON % 20.16E\n",DBL_EPSILON/2.);
  if( 1. == 1.+DBL_EPSILON)
    printf(" Esto esta raro no?\n");
  genArray( n,min,max,xyzInp);

  printXYZ(3*n,xyzInp,"H","INPUT.xyz");

  deleteRepeated(n,xyzInp,xyzOut,&m);


  printXYZ(m,xyzOut,"He","OUTPUT.xyz");
  

  exit(EXIT_FAILURE);
}



double randomGen(double min, double max){

  double val;

  val = (max-min)*drand48() + min;

  return val;
}


void   genArray(int n, double min, double max, double *array){
  int i;
  double x,y,z;

  for(i=0;i<n;i++){
    x = randomGen(min,max);
    y = randomGen(min,max);
    z = randomGen(min,max);
#ifdef SPHERE
    double theta,phi,rr;
    rr    = (max - min)/2.;
    theta = randomGen(0,M_PI);
    phi   = randomGen(0,2*M_PI);
    x     = rr*sin(theta)*cos(phi)+min;
    y     = rr*sin(theta)*sin(phi)+min;
    z     = rr*cos(theta)+min;
#endif
    array[9*i]   = x;
    array[9*i+1] = y;
    array[9*i+2] = z;
    array[9*i+3] = x;
    array[9*i+4] = y;
    array[9*i+5] = z;
    array[9*i+6] = x;
    array[9*i+7] = y;
    array[9*i+8] = z;
  }

}

void printXYZ(int n, double *array, const char *sym,const char *name){
  
  int i;
  FILE *out;

  out = fopen(name,"w+");
  if(out == NULL ) {
    printf(" No se pudo abrir el archivo  [%s]\n",name);
    exit(EXIT_FAILURE);
  }

  fprintf(out,"  %3d\n",n);
  fprintf(out," archivo con posibles puntos [%s]\n",name);
  for(i=0;i<n;i++)
    fprintf(out,"  %4s   % 10.6lf % 10.6lf % 10.6lf\n",sym,array[3*i],array[3*i+1],array[3*i+2]);


  fclose(out);
}

void deleteRepeated(int n,double *arrayInp,double *arrayOut,int *m){
  
  int i,j,k;
  n *= 3;
  double vecAux[9*n];
  double vectmp[3];
  /*
  (*m) = n;
  for(i=0;i<n;i++){
    arrayOut[3*i  ] = arrayInp[3*i  ];
    arrayOut[3*i+1] = arrayInp[3*i+1];
    arrayOut[3*i+2] = arrayInp[3*i+2];
  }
  printf(" Entraron : %d \n",n);
  printf(" Salieron : %d \n",(*m));
  */
  int cont;
  int z;
  j=z=0;
  for(i=0;i<n;i++){
    cont = 0;
    vectmp[0]     = arrayInp[3*i];
    vectmp[1]     = arrayInp[3*i+1];
    vectmp[2]     = arrayInp[3*i+2];
    vecAux[3*j]   = vectmp[0];
    vecAux[3*j+1] = vectmp[1];
    vecAux[3*j+2] = vectmp[2];
    j++;
    for(k=0;k<j;k++)
      if( vecAux[3*k] == vectmp[0] && vecAux[3*k+1] == vectmp[1] && vecAux[3*k+2] == vectmp[2]){
        cont++;
      }
    if( cont == 1){
      arrayOut[3*z  ] = vectmp[0];
      arrayOut[3*z+1] = vectmp[1];
      arrayOut[3*z+2] = vectmp[2];
      z++;
    }


  }
  (*m) = z;
  printf(" Entraron : %d \n",n);
  printf(" Salieron : %d \n",(*m));

  // Termina
}
