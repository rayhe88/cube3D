#include "matvec.h"


// Return  vector's norm for a 3D-vector
double getNormVec(double vec[3]){

  double valor;

  valor = vec[0]*vec[0];
  valor += vec[1]*vec[1];
  valor += vec[2]*vec[2];

  valor = sqrt(valor);

  return valor;

}

// Return  dot product between 2 3D-vectors 
double dotProduct(double vecA[3],double vecB[3]){

  double valor;

  valor = vecA[0]*vecB[0];
  valor += vecA[1]*vecB[1];
  valor += vecA[2]*vecB[2];

  return valor;

}

// This function calculates the cross product between 2 3D-vectors
void crossProduct(double vecA[3],double vecB[3],double vecOut[3]){

  vecOut[0] = vecA[1]*vecB[2]-vecA[2]*vecB[1];
  vecOut[1] = vecA[2]*vecB[0]-vecA[0]*vecB[2];
  vecOut[2] = vecA[0]*vecB[1]-vecA[1]*vecB[0];

}

// This function calculates the product between vector and matrix 3x3.
void matVecProduct(double v[3],double m[9],double vecOut[3]){

  vecOut[0] = v[0]*m[0] + v[1]*m[1] + v[2]*m[2];
  vecOut[1] = v[0]*m[3] + v[1]*m[4] + v[2]*m[5];
  vecOut[2] = v[0]*m[6] + v[1]*m[7] + v[2]*m[8];

  
}

double detMat(double m[9]){

  double det;

  det  = m[0]*(m[4]*m[8] - m[5]*m[7]);
  det += m[1]*(m[5]*m[6] - m[3]*m[8]);
  det += m[2]*(m[3]*m[7] - m[4]*m[7]);

  return det;

}
//
// Return eigen values and eigenvector for a matrix 
//       | a  b  c |
// mat = | 0  d  e |
//       | 0  0  g |
// 
//
int eigenVV(double mat[9],double val[3],double vec[9]){
  
  if( fabs(mat[3]) > 1.E-7 )  return 3;
  if( fabs(mat[6]) > 1.E-7 )  return 6;
  if( fabs(mat[7]) > 1.E-7 )  return 7;

  double a,b,c,d,e,g;

  a = mat[0]; b = mat[1]; c = mat[2];
  d = mat[4]; e = mat[5]; g = mat[8];

  val[0] = a; val[1] = d; val[2] = g;

  vec[0] = (double) 1.;
  vec[1] = (double) 0.;
  vec[2] = (double) 0.;

  vec[3] = b/(d-a);
  vec[4] = (double) 1.;
  vec[5] = (double) 0.;

  vec[6] = (c*d-b*e-c*g)/( (a-g)*(g-d));
  vec[7] = e/(g-d);
  vec[8] = (double) 1.;

  if( fabs(b) < 1.E-7 && fabs(c) < 1.E-7 && fabs(e) < 1.E-7){
    vec[0] = (double) 1.;
    vec[1] = (double) 0.;
    vec[2] = (double) 0.;

    vec[3] = (double) 0.;
    vec[4] = (double) 1.;
    vec[5] = (double) 0.;

    vec[6] = (double) 0.;
    vec[7] = (double) 0.;
    vec[8] = (double) 1.;
  }

  return 0;
}

void getMatT(double *mvec, double *matT,double *matTinv){

  double vecA[3],vecB[3],vecC[3];
  double normA,normB,normC;
  double vol;
  double alpha,beta,gamma;
  double vecTmp[3];

  vecA[0] = mvec[0]; vecA[1] = mvec[1]; vecA[2] = mvec[2];
  vecB[0] = mvec[3]; vecB[1] = mvec[4]; vecB[2] = mvec[5];
  vecC[0] = mvec[6]; vecC[1] = mvec[7]; vecC[2] = mvec[8];

  normA = getNormVec(vecA);
  normB = getNormVec(vecB);
  normC = getNormVec(vecC);

  double tmp;

  tmp   = dotProduct(vecC,vecB);
  tmp  /= (normC*normB);
  alpha = acos(tmp);

  tmp   = dotProduct(vecC,vecA);
  tmp  /= (normC*normA);
  beta  = acos(tmp);

  tmp   = dotProduct(vecA,vecB);
  tmp  /= (normA*normB);
  gamma = acos(tmp);

  crossProduct(vecA,vecB,vecTmp);

  vol = fabs(dotProduct(vecC,vecTmp));



#ifdef DEBUG
  printf(" Norm of vecA =  % 10.6lf\n",normA);
  printf(" Norm of vecB =  % 10.6lf\n",normB);
  printf(" Norm of vecC =  % 10.6lf\n\n",normC);

  printf(" Angle of alpha = % 10.6lf rad  % 6.2lf grad\n",alpha,GRAD(alpha));
  printf(" Angle of beta  = % 10.6lf rad  % 6.2lf grad\n",beta ,GRAD(beta) );
  printf(" Angle of gamma = % 10.6lf rad  % 6.2lf grad\n\n",gamma,GRAD(gamma));

  printf(" Cell volume    = % 10.6lf\n\n",vol);
#endif

  matT[0] = 1.;
  matT[1] = cos(gamma);
  matT[2] = cos(beta);

  matT[3] = 0.;
  matT[4] = sin(gamma);
  matT[5] = cos(alpha)/sin(gamma) - cos(beta)/tan(gamma);
  
  matT[6] = 0.;
  matT[7] = 0.;
  matT[8] = vol/(sin(gamma)*normA*normB*normC);

  tmp = matT[4]*matT[8];

  matTinv[0] = 1.;
  matTinv[1] = -matT[1]/matT[4];
  matTinv[2] = (matT[1]*matT[5]-matT[2]*matT[4])/tmp;

  matTinv[3] = 0.;
  matTinv[4] = 1/matT[4];
  matTinv[5] = -matT[5]/tmp;

  matTinv[6] = 0.;
  matTinv[7] = 0.;
  matTinv[8] = 1/matT[8];

#ifdef DEBUG
  printf("          |% 10.6lf % 10.6lf % 10.6lf|\n",  matT[0],matT[1],matT[2]);
  printf("    T  =  |% 10.6lf % 10.6lf % 10.6lf|\n",  matT[3],matT[4],matT[5]);
  printf("          |% 10.6lf % 10.6lf % 10.6lf|\n\n",matT[6],matT[7],matT[8]);

  printf("    det(T) = % 10.6lf\n\n",detMat(matT));
 
  int info;
  double eval[3],evec[9];

  info = eigenVV(matT,eval,evec);
  if( info != 0 ){
    printf(" Error al usar eigenVV: %3d\n",info);
  }
  
  printf(" Valores y vectores propios\n");
  printf(" Valor : % 10.6lf\n",eval[0]);
  printf(" Vector: % 10.6lf % 10.6lf % 10.6lf\n\n",evec[0],evec[1],evec[2]);
  printf(" Valor : % 10.6lf\n",eval[1]);
  printf(" Vector: % 10.6lf % 10.6lf % 10.6lf\n\n",evec[3],evec[4],evec[5]);
  printf(" Valor : % 10.6lf\n",eval[2]);
  printf(" Vector: % 10.6lf % 10.6lf % 10.6lf\n\n",evec[6],evec[7],evec[8]);
  
  printf("          |% 10.6lf % 10.6lf % 10.6lf|\n",  matTinv[0],matTinv[1],matTinv[2]);
  printf(" T^(-1) = |% 10.6lf % 10.6lf % 10.6lf|\n",  matTinv[3],matTinv[4],matTinv[5]);
  printf("          |% 10.6lf % 10.6lf % 10.6lf|\n\n",matTinv[6],matTinv[7],matTinv[8]);

  printf(" det(T^-1) = % 10.6lf\n",detMat(matTinv));


  info = eigenVV(matTinv,eval,evec);
  if( info != 0 ){
    printf(" Error al usar eigenVV: %3d\n",info);
  }
  
  printf(" Valores y vectores propios\n");
  printf(" Valor : % 10.6lf\n",eval[0]);
  printf(" Vector: % 10.6lf % 10.6lf % 10.6lf\n\n",evec[0],evec[1],evec[2]);
  printf(" Valor : % 10.6lf\n",eval[1]);
  printf(" Vector: % 10.6lf % 10.6lf % 10.6lf\n\n",evec[3],evec[4],evec[5]);
  printf(" Valor : % 10.6lf\n",eval[2]);
  printf(" Vector: % 10.6lf % 10.6lf % 10.6lf\n\n",evec[6],evec[7],evec[8]);

#endif


}


void transform(dataCube cube1,dataCube *cube2,double *coor,double *coor2,double matInv[9]){

   int i,j,natm;
   double vecTmp[3];
   double vecOut[3];

   natm = cube1.natm;
   cube2->natm = natm;
   cube2->npt = cube1.npt;

   
   for(i=0;i<3;i++)
     vecTmp[i] = cube1.min[i]; 
   matVecProduct(vecTmp,matInv,vecOut);
   for(i=0;i<3;i++){
     cube2->min[i] = vecOut[i];
     cube2->pts[i] = cube1.pts[i];
   }


   for(i=0;i<3;i++)
     vecTmp[i] = cube1.mvec[i];

   matVecProduct(vecTmp,matInv,vecOut);
   double htmp;
   htmp = 0.;
   for(i=0;i<3;i++){
     cube2->mvec[i] = vecOut[i];
     htmp += vecOut[i]*vecOut[i];
   }
   cube2->hvec[0] = sqrt(htmp);

   for(i=0;i<3;i++)
     vecTmp[i] = cube1.mvec[3+i];
   matVecProduct(vecTmp ,matInv,vecOut);

   htmp = 0.;
   for(i=0;i<3;i++){
     cube2->mvec[3+i] = vecOut[i];
     htmp += vecOut[i]*vecOut[i];
   }
   cube2->hvec[1] = sqrt(htmp);

   for(i=0;i<3;i++)
     vecTmp[i] = cube1.mvec[6+i];
   matVecProduct(vecTmp,matInv,vecOut);
   htmp = 0.;
   for(i=0;i<3;i++){
     cube2->mvec[6+i] = vecOut[i];
     htmp += vecOut[i]*vecOut[i];
   }
   cube2->hvec[2] = sqrt(htmp);

   
   for( i = 0; i < natm; i++){
     for( j = 0; j < 3; j++)
       vecTmp[j] = coor[3*i+j];

     matVecProduct(vecTmp,matInv,vecOut);

     for( j = 0; j < 3; j++)
       coor2[3*i+j] = vecOut[j];
     

   }

}


/*
void barra(int porc , int n) {
  int i;
  int npor;

  npor = (n*porc/100);

  printf(" [");
  for(i=0;i<n;i++){
    if( i < npor)
      putchar('=');
    else
      putchar(' ');

  }
  printf("] %4d%% Completado...",porc);
  printf("\n");

}

*/
