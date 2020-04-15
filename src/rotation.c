#include "rotation.h"
#include "mathTools.h"
#include "utils.h"


// 1 ) getAngles(vecOriginales,angles);

// 2 ) rotacionZYX(vecIn,angles,vecOut);
//
int checkRotation( double *mvec){

  int ret=NOT;
  const double zero = 1.E-8;
  double ay = mvec[1];
  double az = mvec[2];
  double bz = mvec[5];

  printBanner("Checking Axes (rotation)",stdout);

  if( fabs(ay) > zero || fabs(az) > zero || fabs(bz) > zero )
    ret = YES;
  
  return ret;
}

void rotationCube(double *mvec,dataCube* cube){
  int i;
  double angles[3];
  double r[3],rout[3];

  getAngles(mvec,angles);
  printf(" WARNING! The cube files needs to be rotated \n\n");
  printf("   Rotation on the X axis : %8.2lf deg\n",angles[0]*180./M_PI);
  printf("   Rotation on the Y axis : %8.2lf deg\n",angles[1]*180./M_PI);
  printf("   Rotation on the Z axis : %8.2lf deg\n",angles[2]*180./M_PI);
  printBar(stdout);

  r[0] = mvec[0];
  r[1] = mvec[1];
  r[2] = mvec[2];
  transform(r,angles,rout);
  cube->mvec[0] = rout[0];
  cube->mvec[1] = rout[1];
  cube->mvec[2] = rout[2];
  cube->hvec[0] = getNormVec(rout);

  r[0] = mvec[3];
  r[1] = mvec[4];
  r[2] = mvec[5];
  transform(r,angles,rout);
  cube->mvec[3] = rout[0];
  cube->mvec[4] = rout[1];
  cube->mvec[5] = rout[2];
  cube->hvec[1] = getNormVec(rout);

  r[0] = mvec[6];
  r[1] = mvec[7];
  r[2] = mvec[8];
  transform(r,angles,rout);
  cube->mvec[6] = rout[0];
  cube->mvec[7] = rout[1];
  cube->mvec[8] = rout[2];
  cube->hvec[2] = getNormVec(rout);

  r[0] = cube->min[0];
  r[1] = cube->min[1];
  r[2] = cube->min[2];
  transform(r,angles,rout);
  cube->min[0] = rout[0];
  cube->min[1] = rout[1];
  cube->min[2] = rout[2];

  cube->max[0] = cube->min[0] + (cube->pts[0])*(cube->hvec[0]);
  cube->max[1] = cube->min[1] + (cube->pts[1])*(cube->hvec[1]);
  cube->max[2] = cube->min[2] + (cube->pts[2])*(cube->hvec[2]);

  for( i = 0; i < cube->natm ; i++){
    r[0] = cube->coor[3*i];
    r[1] = cube->coor[3*i+1];
    r[2] = cube->coor[3*i+2];
    transform(r,angles,rout);
    cube->coor[3*i]   = rout[0];
    cube->coor[3*i+1] = rout[1];
    cube->coor[3*i+2] = rout[2];
  }


}

void getAngles( double *vec, double *theta){

  int i;
  double tmp1[3],tmp2[3],tmp3[3];
  double axeX[3],axeY[3];

  double axe1[3],axe2[3],axe3[3];

  double axe4[3],axe5[3],axe6[3];
  double axe7[3],axe8[3],axe9[3];
  double axe10[3],axe11[3],axe12[3];

  axe1[0] = vec[0];  axe1[1] = vec[1];  axe1[2] = vec[2];
  axe2[0] = vec[3];  axe2[1] = vec[4];  axe2[2] = vec[5];
  axe3[0] = vec[6];  axe3[1] = vec[7];  axe3[2] = vec[8];

  for(i=0; i<3; i++){
    tmp1[i] = tmp2[i] = tmp3[i] = (double) 0.0;
  }

  axeX[0] = 1.0; axeX[1] = 0.0; axeX[2] = 0.0;
  axeY[0] = 0.0; axeY[1] = 1.0; axeY[2] = 0.0;

  // primer angulo en Z
  tmp1[0] = axe1[0];  tmp1[1] = axe1[1];  tmp1[2] = 0.;

  theta[2] = getTheta(2,tmp1,axeX);
  rotInZ ( axe1,theta[2], axe4);
  rotInZ ( axe2,theta[2], axe5);
  rotInZ ( axe3,theta[2], axe6);

  // segundo angulo en Y
  tmp2[0] = axe4[0];  tmp2[1] = 0.0;  tmp2[2] = axe4[2];

  theta[1] = getTheta(1,tmp2,axeX);
  rotInY ( axe4,theta[1], axe7);
  rotInY ( axe5,theta[1], axe8);
  rotInY ( axe6,theta[1], axe9);

  // tercer angulo en X
  tmp3[0] = 0.0;  tmp3[1] = axe8[1];  tmp3[2] = axe8[2];

  theta[0] = getTheta(0,tmp3,axeY);
  rotInX ( axe7,theta[0], axe10);
  rotInX ( axe8,theta[0], axe11);
  rotInX ( axe9,theta[0], axe12);

}

double getTheta( int i, double *vecV, double *vecU ){
  double vecR[3];
  double normV,normU;
  double prod,thet;

  normV = getNormVec(vecV);
  normU = getNormVec(vecU);


  prod = dotProduct(vecV,vecU);

  prod /= ( normV*normU);

  thet = acos(prod);

  vecR[0] = vecU[2]*vecV[1] - vecU[1]*vecV[2];
  vecR[1] = vecU[0]*vecV[2] - vecU[2]*vecV[0];
  vecR[2] = vecU[1]*vecV[0] - vecU[0]*vecV[1];

  double sign=1.;

  if( vecR[i] >= 0. ) 
    sign = 1.;
  else
    sign = -1.;


  thet *= sign;


  return thet;

}

void rotInX ( double *v, double theta, double *out){
  double cx = cos(theta);
  double sx = sin(theta);

  out[0] =     v[0];
  out[1] =  cx*v[1] - sx*v[2];
  out[2] =  sx*v[1] + cx*v[2];

}
void rotInY ( double *v, double theta, double *out){
  double cy = cos(theta);
  double sy = sin(theta);

  out[0] =  cy*v[0] + sy*v[2];
  out[1] =    v[1];
  out[2] = -sy*v[0] + cy*v[2];

}
void rotInZ ( double *v, double theta, double *out){
  double cz = cos(theta);
  double sz = sin(theta);

  out[0] =  cz*v[0] - sz*v[1];
  out[1] =  sz*v[0] + cz*v[1];
  out[2] =     v[2];

}

void transform (double *v, double *theta, double *out){
 // out = rotx.roty.rotz.v
 // Ayuda a pasar
 //  ( 1  bx  cx )             ( a b c )
 //  ( 0  by  cy ) = transform.( d e f )  
 //  ( 0   0  cz )             ( g h i )

  double cx,sx,cy,sy,cz,sz;

  cx = cos(theta[0]); sx = sin(theta[0]);
  cy = cos(theta[1]); sy = sin(theta[1]);
  cz = cos(theta[2]); sz = sin(theta[2]);

  out[0] = cy*(cz*v[0] - sz*v[1]) + sy*v[2];
  out[1] = cx*(sz*v[0] + cz*v[1]) - sx*(-sy*(cz*v[0]-sz*v[1]) + cy*v[2]);
  out[2] = sx*(sz*v[0] + cz*v[1]) + cx*(-sy*(cz*v[0]-sz*v[1]) + cy*v[2]);
}

void transformInv (double *v, double *theta, double *out){
// out = rotz.roty.rotx.v
 // Ayuda a pasar
 //  ( a  b  c )                ( 1  by  cx )
 //  ( d  e  f ) = transformInv.( 0  bx  cy )  
 //  ( g  h  i )                ( 0   0  cz )

  double cx,sx,cy,sy,cz,sz;

  cx = cos(theta[0]); sx = sin(theta[0]);
  cy = cos(theta[1]); sy = sin(theta[1]);
  cz = cos(theta[2]); sz = sin(theta[2]);

  out[0] = -sz*(cx*v[1] - sx*v[2]) + cz*(cy*v[0] + sy*(sx*v[1] + cx*v[2]));
  out[1] =  cz*(cx*v[1] - sx*v[2]) + sz*(cy*v[0] + sy*(sx*v[1] + cx*v[2]));
  out[2] = -sy*v[0] + cy*(sx*v[1] + cx*v[2]);

}
