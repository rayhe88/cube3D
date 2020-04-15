/**
 * @file   transU.c
 * @brief
 * @author Raymundo Hernández-Esparza.
 * @date   August 2018.
 */
#include "transU.h"
#include "lectura.h"

/**
 * @brief
 * @param cube
 * @param *rec
 * @param *matU
 */
void getMatT(dataCube cube, int *rec, double *matT){

  int ntot,val;
  double normA,normB,normC;
  double alpha,beta,gamma,vol,tmp;
  double vecA[3],vecB[3],vecC[3],vecTmp[3];

  ntot = (cube.pts[0] - 1)*(cube.pts[1] - 1)*(cube.pts[2] - 1);
  vecA[0] = cube.mvec[0]; vecA[1] = cube.mvec[1]; vecA[2] = cube.mvec[2];
  vecB[0] = cube.mvec[3]; vecB[1] = cube.mvec[4]; vecB[2] = cube.mvec[5];
  vecC[0] = cube.mvec[6]; vecC[1] = cube.mvec[7]; vecC[2] = cube.mvec[8];

  // Obtenemos las normas de los vectores
  normA = getNormVec(vecA);
  normB = getNormVec(vecB);
  normC = getNormVec(vecC);

  // Obtenemos los ángulos
  tmp   = dotProduct(vecC,vecB);
  tmp  /= (normC*normB);
  alpha = acos(tmp);

  tmp   = dotProduct(vecC,vecA);
  tmp  /= (normC*normA);
  beta  = acos(tmp);

  tmp   = dotProduct(vecA,vecB);
  tmp  /= (normA*normB);
  gamma = acos(tmp);

  // Obtenemos el volumen de una celdilla
  crossProduct(vecA,vecB,vecTmp);

  vol = fabs(dotProduct(vecC,vecTmp));

  // Matriz de transformacion T
  matT[0] = 1.;
  matT[1] = cos(gamma);
  matT[2] = cos(beta);

  matT[3] = 0.;
  matT[4] = sin(gamma);
  matT[5] = cos(alpha)/sin(gamma) - cos(beta)/tan(gamma);
  
  matT[6] = 0.;
  matT[7] = 0.;
  matT[8] = vol/(sin(gamma)*normA*normB*normC);
  // Matriz de transformacion T version 2
  /* 
  matT[0] = vecA[0]/normA;
  matT[1] = vecB[0]/normB;
  matT[2] = vecC[0]/normC;

  matT[3] = vecA[1]/normA;
  matT[4] = vecB[1]/normB;
  matT[5] = vecC[1]/normC;

  matT[6] = vecA[2]/normA;
  matT[7] = vecB[2]/normB;
  matT[8] = vecC[2]/normC;
  */

  printBanner(" Geometrical Info ",stdout);
  printf("      vector A   =  % 10.6lf  % 10.6lf  % 10.6lf\n",vecA[0],vecA[1],vecA[2]);
  printf("      vector B   =  % 10.6lf  % 10.6lf  % 10.6lf\n",vecB[0],vecB[1],vecB[2]);
  printf("      vector C   =  % 10.6lf  % 10.6lf  % 10.6lf\n",vecC[0],vecC[1],vecC[2]);
  printf("\n");
  printf("   Norm of vecA  =  % 10.6lf\n",normA);
  printf("   Norm of vecB  =  % 10.6lf\n",normB);
  printf("   Norm of vecC  =  % 10.6lf\n\n",normC);

  printf("   Angle alpha   = % 10.6lf rad  % 6.2lf deg\n",alpha,GRAD(alpha));
  printf("   Angle beta    = % 10.6lf rad  % 6.2lf deg\n",beta ,GRAD(beta) );
  printf("   Angle gamma   = % 10.6lf rad  % 6.2lf deg\n\n",gamma,GRAD(gamma));

  if( detMat(matT) == 1.0 ){
    printf("    Orthogonal system\n");
    val = 1;
  }else{
    printf("    No-orthogonal system\n");
    val = 0;
  }

  printf("\n   Cell volume   = % 10.6lf (% 10.6lf)\n\n",vol*ntot,vol);


  printf("             |% 10.6lf % 10.6lf % 10.6lf|\n",  matT[0],matT[1],matT[2]);
  printf("       T  =  |% 10.6lf % 10.6lf % 10.6lf|\n",  matT[3],matT[4],matT[5]);
  printf("             |% 10.6lf % 10.6lf % 10.6lf|\n\n",matT[6],matT[7],matT[8]);
  printf("   det(T) = % 10.6lf\n\n",detMat(matT));

  printBar(stdout);

  (*rec) = val;
}

void transformCube (dataCube inp, dataCube *out,
                    int **zatm2, double **coor2,
                    double **field2, double* matU){

  int i;
  double tmp[3];
  double tmvec[9];

  out->natm = inp.natm;
  out->npt  = inp.npt;

  for(i=0;i<3;i++)
    out->pts[i] = inp.pts[i];

  tmp[0] = inp.min[0]; tmp[1] = inp.min[1]; tmp[2] = inp.min[2];
  trans00(tmp,matU);
  out->min[0] = tmp[0]; 
  out->min[1] = tmp[1]; 
  out->min[2] = tmp[2];


  tmp[0] = inp.max[0]; tmp[1] = inp.max[1]; tmp[2] = inp.max[2];
  trans00(tmp,matU);
  out->max[0] = tmp[0];
  out->max[1] = tmp[1]; 
  out->max[2] = tmp[2];

  tmp[0] = inp.mvec[0]; tmp[1] = inp.mvec[1]; tmp[2] = inp.mvec[2];
  trans00(tmp,matU);
  tmvec[0] = tmp[0];  tmvec[1] = tmp[1];  tmvec[2] = tmp[2];

  tmp[0] = inp.mvec[3]; tmp[1] = inp.mvec[4]; tmp[2] = inp.mvec[5];
  trans00(tmp,matU);
  tmvec[3] = tmp[0];  tmvec[4] = tmp[1];  tmvec[5] = tmp[2];

  tmp[0] = inp.mvec[6]; tmp[1] = inp.mvec[7]; tmp[2] = inp.mvec[8];
  trans00(tmp,matU);
  tmvec[6] = tmp[0];  tmvec[7] = tmp[1];  tmvec[8] = tmp[2];

  out->mvec[0] = tmvec[0];  out->mvec[1] = tmvec[1];  out->mvec[2] = tmvec[2];
  out->mvec[3] = tmvec[3];  out->mvec[4] = tmvec[4];  out->mvec[5] = tmvec[5];
  out->mvec[6] = tmvec[6];  out->mvec[7] = tmvec[7];  out->mvec[8] = tmvec[8];

  out->hvec[0] = sqrt( pow(tmvec[0],2) +  pow(tmvec[1],2) +  pow(tmvec[2],2) );
  out->hvec[1] = sqrt( pow(tmvec[3],2) +  pow(tmvec[4],2) +  pow(tmvec[5],2) );
  out->hvec[2] = sqrt( pow(tmvec[6],2) +  pow(tmvec[7],2) +  pow(tmvec[8],2) );
  

  createArrays(zatm2,coor2,field2,inp.natm,inp.npt);

  out->zatm = (*zatm2);
  out->coor = (*coor2);
  out->field = (*field2);
  
  for(i=0;i<inp.natm;i++){
    out -> zatm[i] = inp.zatm[i];
    tmp[0]   = inp.coor[3*i];
    tmp[1]   = inp.coor[3*i+1];
    tmp[2]   = inp.coor[3*i+2];
    trans00(tmp,matU);
    out -> coor[3*i]   = tmp[0];
    out -> coor[3*i+1] = tmp[1];
    out -> coor[3*i+2] = tmp[2];
  }

  for(i=0;i<inp.npt;i++)
    out -> field[i] = inp.field[i];

}


/**
 * @brief
 * @param
 * @param
 */
void getMatInv(double *matT,double *matU){

  matU[0] = 1.;
  matU[1] = -matT[1]/matT[4];
  matU[2] = (matT[1]*matT[5]-matT[2]*matT[4])/(matT[4]*matT[8]);

  matU[3] = 0.;
  matU[4] = 1./matT[4];
  matU[5] = -matT[5]/(matT[4]*matT[8]);

  matU[6] = 0.;
  matU[7] = 0.;
  matU[8] = 1/matT[8];
  
  /*
  printf(" ------------------------------------------------------\n");
  printf("              MATRIZ DE TRANSFORMACION \n");
  printf(" ------------------------------------------------------\n");
  printf("         | % 10.6lf  % 10.6lf  % 10.6lf|\n",matT[0],matT[1],matT[2]);
  printf("   T   = | % 10.6lf  % 10.6lf  % 10.6lf|\n",matT[3],matT[4],matT[5]);
  printf("         | % 10.6lf  % 10.6lf  % 10.6lf|\n",matT[6],matT[7],matT[8]);
  printf(" ------------------------------------------------------\n");
  printf("              MATRIZ DE TRANSFORMACION  INVERSA\n");
  printf(" ------------------------------------------------------\n");
  printf("         | % 10.6lf  % 10.6lf  % 10.6lf|\n",matU[0],matU[1],matU[2]);
  printf("   U   = | % 10.6lf  % 10.6lf  % 10.6lf|\n",matU[3],matU[4],matU[5]);
  printf("         | % 10.6lf  % 10.6lf  % 10.6lf|\n",matU[6],matU[7],matU[8]);
  printf(" ------------------------------------------------------\n");
  */

}

/**
 * @brief
 * @param
 * @param
 *        |  1   u1   u2 |
 *   u  = |  0   u4   u5 |
 *        |  0    0   u8 |
 *
 *   r  = ( rx   ry   rz )
 *
 * return grad
 */

void itrans00(double *vec, const double *u){
  
  double rx,ry,rz;

  rx = vec[0]; ry = vec[1]; rz = vec[2];

  vec[0]  = rx - ry*(u[1]/u[4]);
  vec[0] += rz*((u[1]*u[5] - u[2]*u[4])/(u[4]*u[8]));
  vec[1]  = ry/u[4] - rz*u[5]/(u[4]*u[8]);
  vec[2]  = rz/u[8];


}
/**
 * @brief
 * @param
 * @param
 *        |  1   u1   u2 |
 *   u  = |  0   u4   u5 |
 *        |  0    0   u8 |
 *
 *   r  = ( rx   ry   rz )
 *
 */
void trans00(double *vec, const double *u){
  
  double rx,ry,rz;

  rx = vec[0]; ry = vec[1]; rz = vec[2];

  vec[0]  = rx + ry*u[1] + rz*u[2];
  vec[1]  = ry*u[4] + rz*u[5];
  vec[2]  = rz*u[8];

}
/**
 * @brief
 * @param
 * @param
 *        |  1   u1   u2 |
 *   u  = |  0   u4   u5 |
 *        |  0    0   u8 |
 *
 * grad = ( v1 , v2 , v3 )
 * 
 * return grad.u
 */
void trans01(double *val, const double *u){
  
  double v1,v2,v3;

  v1 = val[1]; v2 = val[2]; v3 = val[3];

  val[1] = v1;
  val[2] = u[1]*v1 + u[4]*v2;
  val[3] = u[2]*v1 + u[5]*v2 + u[8]*v3;

}
/**
 * @brief
 * @param
 * @param
 *        |  1   u1   u2 |
 *   u  = |  0   u4   u5 |
 *        |  0    0   u8 |
 *
 * grad = ( v1 , v2 , v3 )
 *
 *        | v4   v7   v8 |
 *   h  = | v7   v5   v9 |
 *        | v8   v9   v6 |
 *
 */
void trans02(double *val, const double *u){
  
  double v1,v2,v3;

  double v4,v5,v6;
  double v7,v8,v9;
//--------- Transform the gradient
  v1 = val[1]; v2 = val[2]; v3 = val[3];

  val[1] = v1;
  val[2] = u[1]*v1 + u[4]*v2;
  val[3] = u[2]*v1 + u[5]*v2 + u[8]*v3;
//--------- Transform the hessian
  v4 = val[4]; v5 = val[5]; v6 = val[6];
  v7 = val[7]; v8 = val[8]; v9 = val[9];

  val[4]  = v4;
  val[5]  = u[4]*(u[4]*v5 + u[1]*v7) + u[1]*(u[1]*v4 + u[4]*v7);
  val[6]  = u[2]*(u[2]*v4 + u[5]*v7 + u[8]*v8);
  val[6] += u[8]*(u[8]*v6 + u[2]*v8 + u[5]*v9);
  val[6] += u[5]*(u[5]*v5 + u[2]*v7 + u[8]*v9);

  val[7]  = u[1]*v4 + u[4]*v7;
  val[8]  = u[2]*v4 + u[5]*v7 + u[8]*v8;
  val[9]  = u[1]*(u[2]*v4 + u[5]*v7 + u[8]*v8);
  val[9] += u[4]*(u[5]*v5 + u[2]*v7 + u[8]*v9); 
}
