#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef _JACOBI_H_
#define _JACOBI_H_

#define N 3
#define ZERO 1.E-8
#define EPS 1.E-8
#define MAXIT 500
#define IDX(X, Y, N) ((N) * (X) + (Y))

void IdtMat(int n, double *out);
void CpyMat(int n, double *in, double *out);
void TrsMat(int n, double *in, double *out);
void MulMat(int n, double *in1, double *in2, double *out);
void Mul3Mat(int n, double *in1, double *in2, double *in3, double *out);
void sortEigen(int n, double *vec, double *mat);

int JacobiNxN(double *matH, double *eval, double *evec);
;

#endif
