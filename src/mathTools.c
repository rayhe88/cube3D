#include "mathTools.h"


double mayor( double matAA[][N], int *p, int *q){
  int i, j;

  *p = 1;
  *q = 0;
  double temp = fabs( matAA[*p][*q]);

  for( i = 2; i < N; i++)
    for( j = 0; j < i; j++)
      if( fabs(matAA[i][j]) > temp){
        temp = fabs(matAA[i][j]);
        *p = i;
        *q = j;
      }

  return temp;

}

int jacobi(double matAA[N][N], double valores[N], double eigenvectors[N][N]){
  
  int i, j, p, q, n = 0;
  double cot, sen, tan, cos, max, ip, iq, temp;

  for( i = 0; i < N; i++){
    for( j = 0; j < N; j++){
      eigenvectors[i][j] = 0.0;
      eigenvectors[i][i] = 1.0;
    }
  }
  max = mayor(matAA, &p, &q);
  while( n < NMAX && max > TOL ){
    cot = (matAA[q][q] - matAA[p][p]) / (2.*matAA[p][q]);

    if (cot >= TOL)
      tan = -cot + sqrt(1. + cot*cot);
    else
      tan = -cot - sqrt(1. + cot*cot);

    cos = 1./sqrt(1. + tan*tan);
    sen = tan*cos;

    for( j = 0; j < N; j++){
      if ( j != p && j != q ) {
        temp = matAA[p][j];
        matAA[p][j] = matAA[p][j]*cos - matAA[q][j]*sen;
        matAA[j][p] = matAA[p][j];
        matAA[q][j] = temp*sen + matAA[q][j]*cos;
        matAA[j][q] = matAA[q][j];
      }
    }//For  
    
    matAA[q][q] = matAA[q][q] + tan * matAA[p][q];
    matAA[p][p] = matAA[p][p] - tan * matAA[p][q];
    matAA[p][q] = 0.;
    matAA[q][p] = 0.;

    for( j = 0; j < N; j++){
      ip = cos * eigenvectors[j][p] + sen * eigenvectors[j][q];
      iq = - sen * eigenvectors[j][p] + cos * eigenvectors[j][q];
      eigenvectors[j][p] = ip;
      eigenvectors[j][q] = iq;
    }
    n++;
    max = mayor(matAA, &p, &q);

  }//While

  if ( n == NMAX) printf(" Maximum number of iterations\n");
    j = 0;
    for (i = 0; i < N; i++){
      valores[i] = matAA[i][i];
      if ( valores[i] > 0 )
        j = i;
    }
      if( j != 0){
        sen = eigenvectors[0][0]; //x
        cos = eigenvectors[0][1]; //y
        cot = eigenvectors[0][2]; //z
        eigenvectors[0][0] = eigenvectors[j][0];
        eigenvectors[0][1] = eigenvectors[j][1];
        eigenvectors[0][2] = eigenvectors[j][2];
        eigenvectors[j][0] = sen;
        eigenvectors[j][1] = cos;
        eigenvectors[j][2] = cot;

        sen = valores[0];
        valores[0] = valores[j];
        valores[j] = sen;
      }

  return 0;
}

