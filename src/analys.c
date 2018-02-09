#include "file.h"
#include "array.h"
#include "analys.h"
#include "struct.h"
#include "lectura.h"


int limitsMacro(int pol,int i,int j, int k, dataCube data, double* limits){

  int min,max;

  switch (pol){
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


  limits[0] = data.min[0] + (i+min)*data.hvec[0];
  limits[1] = data.min[0] + (i+max)*data.hvec[0];
  limits[2] = data.min[1] + (j+min)*data.hvec[1];
  limits[3] = data.min[1] + (j+max)*data.hvec[1];
  limits[4] = data.min[2] + (k+min)*data.hvec[2];
  limits[5] = data.min[2] + (k+max)*data.hvec[2];
  
  return 0;
}

int newCube1(dataCube data1, int *zato,
             double *coor, double *field, double *matT,char *name){


  int i,j,k;
  int nu;
  int ip,jp,kp;
  int pts2[3];
  int npt2,mu;
  double min2[3];
  int npx,npy,npz;
  double fun[27];
  double hx,hy,hz;

  FILE *out;
  double *field2;

  strcat(name,"grd.cube");

  openFile(&out,name,"w+");


  hx = data1.hvec[0];
  hy = data1.hvec[1];
  hz = data1.hvec[2];

  min2[0] = data1.min[0] + hx;
  min2[1] = data1.min[1] + hy;
  min2[2] = data1.min[2] + hz;

  npx = data1.pts[0]; npy = data1.pts[1]; npz = data1.pts[2];

  pts2[0] = npx-2;
  pts2[1] = npy-2;
  pts2[2] = npz-2;

  npt2 = (int) pts2[0]*pts2[1]*pts2[2];

  createArrayDou(npt2,&field2,"Campo2");
  
  
  mu = 0;
  for( i = 1 ; i < npx-1 ; i++){
    for( j = 1 ; j < npy-1 ; j++){
      for( k = 1 ; k < npz-1 ; k++){
        nu = 0;
        for(ip = i-1; ip <= i+1; ip++)
          for(jp = j-1; jp <= j+1; jp++)
            for(kp = k-1; kp <= k+1; kp++){
              fun[nu] = field[IDX(ip,jp,kp,npy,npz)];
              nu++;
            }
        
        field2[mu] = fieldNum(hx,hy,hz,fun);
        mu++;
      }
    }
  }


  checkData(data1.natm,pts2,zato,min2,data1.mvec,coor,field2,name,out);

  free(field2); 

  fclose(out);

  return 0;

}

double fieldNum(double hx,double hy, double hz, double *f){

  double grax,gray,graz;
  double num;


  grax = -f[IDX3(0,1,1)] + f[IDX3(2,1,1)];
  gray = -f[IDX3(1,0,1)] + f[IDX3(1,2,1)];
  graz = -f[IDX3(1,1,0)] + f[IDX3(1,1,2)];

  grax /= (2.*hx);
  gray /= (2.*hy);
  graz /= (2.*hz);

  num = grax*grax + gray*gray + graz*graz;

  num = sqrt(num);

#ifdef RED
  double f0,den;
  f0 = f[IDX3(1,1,1)];
  num *= 0.1616204597; 

  den = pow(f0,1.3333333333333333);

  num = num/den;
#endif

  return num;

}


int newCube3(dataCube data1, int *zato,
             double *coor, double *field, double *matT, double *dg,char *name){


  int i,j,k;
  int nu;
  int ip,jp,kp;
  int pts2[3];
  int npt2;
  double min2[3];
  int npx,npy,npz;
  double fun[27];
  double hx,hy,hz;
  int mu;

  FILE *outRed,*outDen,*outData;
  double *red;
  double *den;

  double tmpred,tmpden;
  char denName[120],redName[120];
  

  sprintf(denName,"%sNCIden.cube",name);
  sprintf(redName,"%sNCIred.cube",name);
  strcat(name,".data");


  openFile(&outDen,denName,"w+");
  openFile(&outRed,redName,"w+");
  openFile(&outData,name,"w+");


  hx = data1.hvec[0];
  hy = data1.hvec[1];
  hz = data1.hvec[2];

  min2[0] = data1.min[0] + hx;
  min2[1] = data1.min[1] + hy;
  min2[2] = data1.min[2] + hz;

  npx = data1.pts[0]; npy = data1.pts[1]; npz = data1.pts[2];

  pts2[0] = npx-2;
  pts2[1] = npy-2;
  pts2[2] = npz-2;

  npt2 =(int) pts2[0]*pts2[1]*pts2[2];

  createArrayDou(npt2,&den,"Den");
  createArrayDou(npt2,&red,"Red");

  fprintf(outData,"#  File dat with the data for rho vs grad plot\n");
  fprintf(outData,"#  cube3D project\n");
  fprintf(outData,"#  Cutoff for density*100 : +/- % 10.6lf\n",dg[0]);
  fprintf(outData,"#  Cutoff for reduced grad:     % 10.6lf\n",dg[1]);
  
  
  mu = 0;
  for( i = 1 ; i < npx-1 ; i++){
    for( j = 1 ; j < npy-1 ; j++){
      for( k = 1 ; k < npz-1 ; k++){
        nu = 0;
        for(ip = i-1; ip <= i+1; ip++)
          for(jp = j-1; jp <= j+1; jp++)
            for(kp = k-1; kp <= k+1; kp++){
              fun[nu] = field[IDX(ip,jp,kp,npy,npz)];
              nu++;
            }
          fieldNCI(hx,hy,hz,fun,dg,&tmpred,&tmpden);
          red[mu] = tmpred;
          den[mu] = tmpden;
        if( tmpred <= dg[1] && fabs(tmpden) <= dg[0])
          fprintf(outData," % 10.6lf  % 10.6lf\n",tmpden,tmpred);
        mu++;
      }
    }
  }


  checkData(data1.natm,pts2,zato,min2,data1.mvec,coor,den,name,outDen);
  checkData(data1.natm,pts2,zato,min2,data1.mvec,coor,red,name,outRed);

  fclose(outData);
  fclose(outDen);
  fclose(outRed);

  free(den);
  free(red);

  return 0;

}

int fieldNCI(double hx,double hy, double hz, double *f,double *dg,double *red,double *rho){

  double grax,gray,graz;
  double f0;
  double matA[6];
  double val[3];
  double num,den;


  f0 = f[IDX3(1,1,1)];

  (*red) = 1.E3;
  (*rho) = 0.;

  if( f0 <= dg[0]/100. ){

    grax = -f[IDX3(0,1,1)] + f[IDX3(2,1,1)];
    gray = -f[IDX3(1,0,1)] + f[IDX3(1,2,1)];
    graz = -f[IDX3(1,1,0)] + f[IDX3(1,1,2)];

    grax /= (2.*hx);
    gray /= (2.*hy);
    graz /= (2.*hz);

    num = grax*grax + gray*gray + graz*graz;

    num = sqrt(num);
    num *= 0.1616204597; 

    den = pow(f0,1.3333333333333333);

    num /= den;
    if( num <= dg[1] ){
      matA[0] = (f[IDX3(0,1,1)]- 2.*f[IDX3(1,1,1)] + f[IDX3(2,1,1)])/(hx*hx);
      matA[1] = (f[IDX3(0,0,1)] - f[IDX3(0,2,1)] -f[IDX3(2,0,1)] + f[IDX3(2,2,1)])/(4.*hx*hy);
      matA[2] = (f[IDX3(0,1,0)] - f[IDX3(0,1,2)] -f[IDX3(2,1,0)] + f[IDX3(2,1,2)])/(4.*hx*hz);
      matA[3] = (f[IDX3(1,0,1)] - 2.*f[IDX3(1,1,1)] + f[IDX3(1,2,1)])/(hy*hy);
      matA[4] = (f[IDX3(1,0,0)] - f[IDX3(1,0,2)] -f[IDX3(1,2,0)] + f[IDX3(1,2,2)])/(4.*hy*hz);
      matA[5] = (f[IDX3(1,1,0)] - 2.*f[IDX3(1,1,1)] + f[IDX3(1,1,2)])/(hz*hz);

      valoresPropios3x3(matA,val);

      (*red) = num;
      (*rho) = f0*100.;
      if( val[1] < 0.)
        (*rho) *= -1;
        

    }
 }

  return 0;
}


int valoresPropios3x3(double *matA,double *eigenVal){

  //        | a  b  c |
  //  matA =| b  d  e |
  //        | c  e  f |

  double a, b, c, d, e, g; 
  double aux;
  double p, eig1, eig2, eig3, q, detB, r, phi, pi; 
  int i,j;

  pi = 4.f*atan(1.f);
  a = matA[0];
  b = matA[1];
  c = matA[2];
  d = matA[3];
  e = matA[4];
  g = matA[5];
  p = b*b + c*c + e*e;
  if (p == 0) {
    eig1 = a;
    eig2 = d;
    eig3 = g;
  } else {
    q = (a+d+g)/3.f;
    p = (a - q)*(a - q) + (d - q)*(d - q) + (g - q)*(g - q) + 2.f*p;
    p = sqrt(p/6.f);
    detB = (a-q)*((d-q)*(g-q)-e*e);
    detB = detB - b*(b*(g-q)-e*c);
    detB = detB + c*(b*e-c*(d-q));
    detB = detB/(p*p*p);
    r = detB/2.f;
    if (r <= -1.f)
       phi = pi/3.f;
    else {
      if (r >= 1)
        phi = 0.f;
      else
        phi = acos(r)/3.f;
    }   
    eig1 = q + 2.f*p*cos(phi);
    eig3 = q + 2.f*p*cos(phi + pi * (2.f/3.f));
    eig2 = 3.f*q - eig1 - eig3;
  }   
  eigenVal[0] = eig1;
  eigenVal[1] = eig2;
  eigenVal[2] = eig3;

  for(i=0;i<3;i++)
    for(j=i+1;j<3;j++){
      if( eigenVal[i] > eigenVal[j]){
        aux = eigenVal[i];
        eigenVal[i] = eigenVal[j];
        eigenVal[j] = aux;
      }   
    }   

  return 0;
}

int newCube2(dataCube data1, int *zato,
             double *coor, double *field, double *matT, char *name){


  int i,j,k;
  int nu;
  int ip,jp,kp;
  int pts2[3];
  double min2[3];
  int npx,npy,npz;
  double fun[27];
  double hx,hy,hz;
  double *field2;
  int npt2,mu;

  FILE *out;

  strcat(name,"lap.cube");

  openFile(&out,name,"w+");


  hx = data1.hvec[0];
  hy = data1.hvec[1];
  hz = data1.hvec[2];

  min2[0] = data1.min[0] + hx;
  min2[1] = data1.min[1] + hy;
  min2[2] = data1.min[2] + hz;

  npx = data1.pts[0]; npy = data1.pts[1]; npz = data1.pts[2];

  pts2[0] = npx-2;
  pts2[1] = npy-2;
  pts2[2] = npz-2;

  npt2 = (int)pts2[0]*pts2[1]*pts2[2];

  createArrayDou(npt2,&field2,"Campo2");
  
  
  mu = 0;
  for( i = 1 ; i < npx-1 ; i++){
    for( j = 1 ; j < npy-1 ; j++){
      for( k = 1 ; k < npz-1 ; k++){
        nu = 0;
        for(ip = i-1; ip <= i+1; ip++)
          for(jp = j-1; jp <= j+1; jp++)
            for(kp = k-1; kp <= k+1; kp++){
              fun[nu] = field[IDX(ip,jp,kp,npy,npz)];
              nu++;
            }
        
        field2[mu] = fieldNum2(hx,hy,hz,fun);
        mu++;
      }
    }
  }


  checkData(data1.natm,pts2,zato,min2,data1.mvec,coor,field2,name,out);

  free(field2); 

  fclose(out);

  return 0;

}

double fieldNum2(double hx,double hy, double hz, double *f){

  double lapx,lapy,lapz;

  lapx = (f[IDX3(0,1,1)] - 2.*f[IDX3(1,1,1)] + f[IDX3(2,1,1)])/(hx*hx);
  lapy = (f[IDX3(1,0,1)] - 2.*f[IDX3(1,1,1)] + f[IDX3(1,2,1)])/(hy*hy);
  lapz = (f[IDX3(1,1,0)] - 2.*f[IDX3(1,1,1)] + f[IDX3(1,1,2)])/(hz*hz);
  return lapx+lapy+lapz;

}


int coefficients (double hx,double hy, double hz, double x1,double y1, double z1,
                  double *f, double *c){
  int ix,iy,j;
  double a[27];
  double b[27];
  double x2 = x1*x1;
  double h2x = hx*hx;

  double tmp1,tmp2,tmp3;

  // tmp1 = f000 - 2 f001 + f002
  // tmp2 = f002 - f000
  // tmp3 = f001
  
  for(ix = 0; ix < 3; ix++){
    for(iy = 0; iy < 3; iy++){
   
      tmp1 = f[IDX3(ix,iy,0)] - 2.*f[IDX3(ix,iy,1)] + f[IDX3(ix,iy,2)];
      tmp2 = f[IDX3(ix,iy,2)] -    f[IDX3(ix,iy,0)];
      tmp3 = f[IDX3(ix,iy,1)];

      a[IDX3(ix,iy,2)] = tmp1;
      a[IDX3(ix,iy,1)] = -2.*z1*tmp1 + hz*tmp2;
      a[IDX3(ix,iy,0)] = z1*z1*tmp1 - hz*z1*tmp2 + 2.*hz*hz*tmp3;
    }
  }
  for(j=0;j<27;j++)
    a[j] /= (2.*hz*hz);

  for(ix = 0; ix < 3; ix++){
    j = 9*ix;
    tmp1 = y1*y1 + y1*hy;
    tmp2 = 2.*hy*hy - 2.*y1*y1;
    tmp3 = y1*y1 - y1*hy;
    b[j  ] = tmp1*a[j  ] + tmp2*a[j+3] + tmp3*a[j+6];
    b[j+1] = tmp1*a[j+1] + tmp2*a[j+4] + tmp3*a[j+7];
    b[j+2] = tmp1*a[j+2] + tmp2*a[j+5] + tmp3*a[j+8];
    tmp1 = -2.*y1-hy;
    tmp2 = 4.*y1;
    tmp3 = hy - 2.*y1;
    b[j+3] = tmp1*a[j  ] + tmp2*a[j+3] + tmp3*a[j+6];
    b[j+4] = tmp1*a[j+1] + tmp2*a[j+4] + tmp3*a[j+7];
    b[j+5] = tmp1*a[j+2] + tmp2*a[j+5] + tmp3*a[j+8];
    b[j+6] = a[j  ] - 2.*a[j+3] + a[j+6];
    b[j+7] = a[j+1] - 2.*a[j+4] + a[j+7];
    b[j+8] = a[j+2] - 2.*a[j+5] + a[j+8];
  }
  
  for(j=0;j<27;j++){
    b[j] /= (2.*hy*hy);
    c[j] = 0.;
  }

  c[26] = b[8] - 2.*b[17] + b[26];
  c[25] = b[7] - 2.*b[16] + b[25];
  c[24] = b[6] - 2.*b[15] + b[24];
  c[23] = b[5] - 2.*b[14] + b[23];
  c[22] = b[4] - 2.*b[13] + b[22];
  c[21] = b[3] - 2.*b[12] + b[21];
  c[20] = b[2] - 2.*b[11] + b[20];
  c[19] = b[1] - 2.*b[10] + b[19];
  c[18] = b[0] - 2.*b[ 9] + b[18];

  c[17] = -2.*x1*b[8] - hx*b[8] + 4.*x1*b[17] - 2.*x1*b[26] + hx*b[26];
  c[16] = -2.*x1*b[7] - hx*b[7] + 4.*x1*b[16] - 2.*x1*b[25] + hx*b[25];
  c[15] = -2.*x1*b[6] - hx*b[6] + 4.*x1*b[15] - 2.*x1*b[24] + hx*b[24];
  c[14] = -2.*x1*b[5] - hx*b[5] + 4.*x1*b[14] - 2.*x1*b[23] + hx*b[23];
  c[13] = -2.*x1*b[4] - hx*b[4] + 4.*x1*b[13] - 2.*x1*b[22] + hx*b[22];
  c[12] = -2.*x1*b[3] - hx*b[3] + 4.*x1*b[12] - 2.*x1*b[21] + hx*b[21];
  c[11] = -2.*x1*b[2] - hx*b[2] + 4.*x1*b[11] - 2.*x1*b[20] + hx*b[20];
  c[10] = -2.*x1*b[1] - hx*b[1] + 4.*x1*b[10] - 2.*x1*b[19] + hx*b[19];
  c[ 9] = -2.*x1*b[0] - hx*b[0] + 4.*x1*b[ 9] - 2.*x1*b[18] + hx*b[18];

  c[ 8] = x2*b[8] + x1*hx*b[8] - 2.*x2*b[17] + 2.*h2x*b[17] + x2*b[26] - x1*hx*b[26];
  c[ 7] = x2*b[7] + x1*hx*b[7] - 2.*x2*b[16] + 2.*h2x*b[16] + x2*b[25] - x1*hx*b[25];
  c[ 6] = x2*b[6] + x1*hx*b[6] - 2.*x2*b[15] + 2.*h2x*b[15] + x2*b[24] - x1*hx*b[24];
  c[ 5] = x2*b[5] + x1*hx*b[5] - 2.*x2*b[14] + 2.*h2x*b[14] + x2*b[23] - x1*hx*b[23];
  c[ 4] = x2*b[4] + x1*hx*b[4] - 2.*x2*b[13] + 2.*h2x*b[13] + x2*b[22] - x1*hx*b[22];
  c[ 3] = x2*b[3] + x1*hx*b[3] - 2.*x2*b[12] + 2.*h2x*b[12] + x2*b[21] - x1*hx*b[21];
  c[ 2] = x2*b[2] + x1*hx*b[2] - 2.*x2*b[11] + 2.*h2x*b[11] + x2*b[20] - x1*hx*b[20];
  c[ 1] = x2*b[1] + x1*hx*b[1] - 2.*x2*b[10] + 2.*h2x*b[10] + x2*b[19] - x1*hx*b[19];
  c[ 0] = x2*b[0] + x1*hx*b[0] - 2.*x2*b[ 9] + 2.*h2x*b[ 9] + x2*b[18] - x1*hx*b[18];
  
  for(j=0;j<27;j++)
    c[j] /= (2.*hx*hx);

  return 0;

}

int NumericalCrit01(double x, double y, double z,double *c, double *val){

  double v0,gx,gy,gz;
  double x2,y2,z2;
  x2 = x*x;
  y2 = y*y;
  z2 = z*z;

  v0 = c[0] + c[1]*z + c[2]*z2 + c[3]*y + c[4]*y*z + c[5]*y*z2 + c[6]*y2 +
       c[7]*y2*z + c[8]*y2*z2 + c[9]*x + c[10]*x*z + c[11]*x*z2 + c[12]*x*y +
       c[13]*x*y*z + c[14]*x*y*z2 + c[15]*x*y2 + c[16]*x*y2*z + c[17]*x*y2*z2 +
       c[18]*x2 + c[19]*x2*z + c[20]*x2*z2 + c[21]*x2*y + c[22]*x2*y*z +
       c[23]*x2*y*z2 + c[24]*x2*y2 + c[25]*x2*y2*z + c[26]*x2*y2*z2;

  gx = c[9] + c[10]*z + c[11]*z2 + c[12]*y + c[13]*y*z + c[14]*y*z2 + c[15]*y2 +
       c[16]*y2*z + c[17]*y2*z2 + 2.*c[18]*x + 2.*c[19]*x*z + 2.*c[20]*x*z2 + 
       2.*c[21]*x*y + 2.*c[22]*x*y*z + 2.*c[23]*x*y*z2 + 2.*c[24]*x*y2 +
       2.*c[25]*x*y2*z + 2.*c[26]*x*y2*z2;

  gy = c[3] + c[4]*z + c[5]*z2 + 2.*c[6]*y + 2.*c[7]*y*z + 2.*c[8]*y*z2 +
       c[12]*x + c[13]*x*z + c[14]*x*z2 + 2.*c[15]*x*y + 2.*c[16]*x*y*z +
       2.*c[17]*x*y*z2 + c[21]*x2 + c[22]*x2*z + c[23]*x2*z2 + 2.*c[24]*x2*y +
       2.*c[25]*x2*y*z + 2.*c[26]*x2*y*z2;

  gz = c[1] + 2.*z*c[2] + y*c[4] + 2.*y*z*c[5] + y2*c[7] + 2.*y2*z*c[8] +
       x*c[10] + 2.*x*z*c[11] + x*y*c[13] + 2.*x*y*z*c[14] + 
       x*y2*c[16] + 2.*x*y2*z*c[17] + x2*c[19] + 2.*x2*z*c[20] + x2*y*c[22] +
       2.*x2*y*z*c[23] + x2*y2*c[25] + 2.*x2*y2*z*c[26];

  val[0] = v0;
  val[1] = gx;
  val[2] = gy;
  val[3] = gz;

  return 0;
}

int NumericalCrit02(double x, double y, double z,double *c, double *val){

  double v0,gx,gy,gz;
  double hxx,hxy,hxz,hyy,hyz,hzz;
  double x2,y2,z2;
  x2 = x*x;
  y2 = y*y;
  z2 = z*z;

  v0 = c[0] + c[1]*z + c[2]*z2 + c[3]*y + c[4]*y*z + c[5]*y*z2 + c[6]*y2 +
       c[7]*y2*z + c[8]*y2*z2 + c[9]*x + c[10]*x*z + c[11]*x*z2 + c[12]*x*y +
       c[13]*x*y*z + c[14]*x*y*z2 + c[15]*x*y2 + c[16]*x*y2*z + c[17]*x*y2*z2 +
       c[18]*x2 + c[19]*x2*z + c[20]*x2*z2 + c[21]*x2*y + c[22]*x2*y*z +
       c[23]*x2*y*z2 + c[24]*x2*y2 + c[25]*x2*y2*z + c[26]*x2*y2*z2;

  gx = c[9] + c[10]*z + c[11]*z2 + c[12]*y + c[13]*y*z + c[14]*y*z2 + c[15]*y2 +
       c[16]*y2*z + c[17]*y2*z2 + 2.*c[18]*x + 2.*c[19]*x*z + 2.*c[20]*x*z2 + 
       2.*c[21]*x*y + 2.*c[22]*x*y*z + 2.*c[23]*x*y*z2 + 2.*c[24]*x*y2 +
       2.*c[25]*x*y2*z + 2.*c[26]*x*y2*z2;

  gy = c[3] + c[4]*z + c[5]*z2 + 2.*c[6]*y + 2.*c[7]*y*z + 2.*c[8]*y*z2 +
       c[12]*x + c[13]*x*z + c[14]*x*z2 + 2.*c[15]*x*y + 2.*c[16]*x*y*z +
       2.*c[17]*x*y*z2 + c[21]*x2 + c[22]*x2*z + c[23]*x2*z2 + 2.*c[24]*x2*y +
       2.*c[25]*x2*y*z + 2.*c[26]*x2*y*z2;

  gz = c[1] + 2.*z*c[2] + y*c[4] + 2.*y*z*c[5] + y2*c[7] + 2.*y2*z*c[8] +
       2.*c[25]*x2*y*z + 2.*c[26]*x2*y*z2;

  gz = c[1] + 2.*z*c[2] + y*c[4] + 2.*y*z*c[5] + y2*c[7] + 2.*y2*z*c[8] +
       x*c[10] + 2.*x*z*c[11] + x*y*c[13] + 2.*x*y*z*c[14] + 
       x*y2*c[16] + 2.*x*y2*z*c[17] + x2*c[19] + 2.*x2*z*c[20] + x2*y*c[22] +
       2.*x2*y*z*c[23] + x2*y2*c[25] + 2.*x2*y2*z*c[26];

  hxx = 2.*c[18] + 2.*z*c[19] + 2.*z2*c[20] + 2.*y*c[21] + 2.*y*z*c[22] + 
        2.*y*z2*c[23] + 2.*y2*c[24] + 2.*y2*z*c[25] + 2.*y2*z2*c[26];

  hyy = 2.*c[6] + 2.*z*c[7] + 2.*z2*c[8] + 2.*x*c[15] + 2.*x*z*c[16] + 
        2.*x*z2*c[17] + 2.*x2*c[24] + 2.*x2*z*c[25] + 2.*x2*z2*c[26];
  
  hzz = 2.*c[2] + 2.*y*c[5] + 2.*y2*c[8] + 2.*x*c[11] + 2.*x*y*c[14] +
        2.*x*y2*c[17] + 2.*x2*c[20] + 2.*x2*y*c[23] + 2.*x2*y2*c[26];
  
  hxy = c[12] + z*c[13] + z2*c[14] + 2.*y*c[15] + 2.*y*z*c[16] +
        2.*y*z2*c[17] + 2.*x*c[21] + 2.*x*z*c[22] + 2.*x*z2*c[23] +
        4.*x*y*c[24] + 4.*x*y*z*c[25] + 4.*x*y*z2*c[26];

  hxz = c[10] + 2.*z*c[11] + y*c[13] + 2.*y*z*c[14] + y2*c[16] + 
        2.*y2*z*c[17] + 2.*x*c[19] + 4.*x*z*c[20] + 2.*x*y*c[22] + 
        4.*x*y*z*c[23] + 2.*x*y2*c[25] + 4.*x*y2*z*c[26];

  hyz = c[4] + 2.*z*c[5] + 2.*y*c[7] + 4.*y*z*c[8] + x*c[13] +
        2.*x*z*c[14] + 2.*x*y*c[16] + 4.*x*y*z*c[17] + x2*c[22] + 
        2.*x2*z*c[23] + 2.*x2*y*c[25] + 4.*x2*y*z*c[26];

  val[0] = v0;
  val[1] = gx;
  val[2] = gy;
  val[3] = gz;
  val[4] = hxx;
  val[5] = hxy;
  val[6] = hxz;
  val[7] = hyy;
  val[8] = hyz;
  val[9] = hzz;

  return 0;
}
