#include "jacobi.h"

int JacobiNxN (double *matH, double *eval,double *evec){

  int i,j;
  int k,l,count;
  double matA[9], matP[9], matT[9];
  double tmpA[9], matQ[9];
  double max,tmp,tol,x,y,z,s,c;
  double sumaq, sign;

  IdtMat(N,matQ);

  //matrizH es la del input
  for(i=0;i<9;i++)
    matA[i] = matH[i];
  count=0;
  do{

    k = 0; l = 1;

    max = fabs(matA[IDX(N,k,l)]);

    for(i=0;i<N-1;i++){
      for(j=i+1;j<N;j++){
        tmp = fabs( matA[IDX(N,i,j)] );
        if( tmp > max ){
          k = i;  l = j;
          max = tmp;
        }
      }
    }
    sumaq= 0.;
    for(i=0;i<N;i++)
      sumaq += (matA[IDX(N,i,i)]*matA[IDX(N,i,i)]);

    tol = (EPS*sqrt(sumaq)/(double) N);
    if (max < tol ) break;



    y = matA[IDX(N,k,k)] - matA[IDX(N,l,l)];
    if( fabs(y) < ZERO ){
      c = s = sin(M_PI/4.);
    }else{
      x = 2.*matA[IDX(N,k,l)];
      z = sqrt( x*x + y*y );

      sign = (x/y);
      sign /= fabs(x/y);

      c =      sqrt( (z + y)/(2.*z));
      s = sign*sqrt( (z - y)/(2.*z));
    }

    IdtMat(N,matP);
    matP[IDX(N,k,k)] =  c;  matP[IDX(N,k,l)] =  s;
    matP[IDX(N,l,k)] = -s;  matP[IDX(N,l,l)] =  c;


  
    TrsMat(N,matP,matT);

    Mul3Mat(N,matP,matA,matT,tmpA);
    CpyMat (N,tmpA,matA);
    MulMat (N,matQ,matT,tmpA);
    CpyMat (N,tmpA,matQ);


    count++;

  }while( count < MAXIT);

  eval[0] = matA[0];
  eval[1] = matA[4];
  eval[2] = matA[8];

  sortEigen(N,eval,matQ);

  for(i=0;i<9;i++)
    evec[i] = matQ[i];


  return 0;
}

void sortEigen(int n, double *eval, double *matQ){
  int i,j,k;
  double tmp;

  for(i=0;i<n-1;i++){
    for(j=i+1;j<n;j++){
      if( eval[i] > eval[j] ){
        tmp = eval[i];
        eval[i] = eval[j];
        eval[j] = tmp;

        for(k=0;k<n;k++){
          tmp = matQ[3*i+k];
          matQ[3*i+k] = matQ[3*j+k];
          matQ[3*j+k] = tmp;
        }


      }

    }

  }

}

void IdtMat(int n, double *out){
  int i,j;

  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      if( i == j ){
        out[IDX(n,i,j)] = (double) 1.;
      }else{
        out[IDX(n,i,j)] = (double) 0.;
      }


}

void CpyMat( int n, double *in, double *out){
  int i;

  for(i=0;i<n*n;i++)
    out[i] = in[i];

}

void TrsMat( int n, double *in, double *out){
  int i,j;

  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      out[IDX(n,j,i)] = in[IDX(n,i,j)];


}

void MulMat( int n,double *in1, double *in2, double *out){

  int i,j,k;
  double tmp;

  for(i=0;i<n;i++)
    for(j=0;j<n;j++){
      tmp = (double) 0.;
      for(k=0;k<n;k++)
        tmp += in1[IDX(n,i,k)]*in2[IDX(n,k,j)];

      out[IDX(n,i,j)] = tmp;

    }
     

}

void Mul3Mat( int n, double *in1, double *in2, double *in3, double *out){

  double tmp[n*n];

  MulMat(n,in2,in3,tmp);
  MulMat(n,in1,tmp,out);
}

