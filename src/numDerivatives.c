

//int IDX3(int i, int j, int k){
  //return (i*9 + j*3 + k);
//}

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

      a[IDX3(ix,iy,0)] = tmp1;
      a[IDX3(ix,iy,1)] = -2.*z1*tmp1 + hz*tmp2;
      a[IDX3(ix,iy,2)] = z1*z1*tmp1 - hz*z1*tmp2 + 2.*hz*hz*tmp3;
    }
  }
  for(j=0;j<27;j++)
    a[i] /= (2.*hz*hz);

  for(ix = 0; ix < 3; ixx){
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
    b[j+5] = tmp1*a[j+2] + tmp2*a[j+5] + tmp6*a[j+8];
    b[j+6] = a[j  ] - 2.*a[j+3] + a[j+6];
    b[j+7] = a[j+1] - 2.*a[j+4] + a[j+7];
    b[j+8] = a[j+2] - 2.*a[j+5] + a[j+8];
  }
  
  for(j=0;j<27;j++){
    b[i] /= (2.*hy*hy);
    c[i] = 0.;
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
  c[ 9] = -2.*x1*b[0] - hx*b[1] + 4.*x1*b[ 9] - 2.*x1*b[18] + hx*b[18];

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
    c[i] /= (2.*hx*hx);

  return 0;

}
