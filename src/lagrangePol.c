#include "lagrangePol.h"




int evalPol0 ( int p, double x, double y, double z, double *coef,double *val){

  int i,j,k;
  double c,f0,facx,facy,facz;

  f0 = 0.;
  for( i = p; i >= 0; i--){
    facx = pow(x,i);
    for( j = p; j >= 0; j--){
      facy = pow(y,j);
      for( k = p; k >= 0; k--){
        c    = coef[ i*(p+1)*(p+1) + j*(p+1) + k];
        facz = pow(z,k);

        f0 += c*facx*facy*facz;

      }  
    }
  }

  val[0] = f0;
  return 0;
}

int evalPol1 ( int p, double x, double y, double z, double *coef,double *val){

  int i,j,k;
  double c,f0,grax,gray,graz;
  double facx0,facy0,facz0;
  double facx1,facy1,facz1;

  f0 = grax = gray = graz  = 0.;
  for( i = p; i >= 0; i--){
    facx0 = pow(x,i);
    facx1 = pow(x,i-1)*(double)i;
    for( j = p; j >= 0; j--){
      facy0 = pow(y,j);
      facy1 = pow(y,j-1)*(double)j;
      for( k = p; k >= 0; k--){
        c    = coef[ i*(p+1)*(p+1) + j*(p+1) + k ];
        
        facz0 = pow(z,k);
        facz1 = pow(z,k-1)*(double)k;

        f0   += c*facx0*facy0*facz0;
        grax += c*facx1*facy0*facz0;
        gray += c*facx0*facy1*facz0;
        graz += c*facx0*facy0*facz1;
      }  
    }
  }

  val[0] = f0;
  val[1] = grax;
  val[2] = gray;
  val[3] = graz;

  return 0;
}

int evalPol2 ( int p, double x, double y, double z, double *coef,double *val){

  int i,j,k;
  double c,f0,grax,gray,graz;
  double hexx,hexy,hexz;
  double heyy,heyz,hezz;
  double facx0,facy0,facz0;
  double facx1,facy1,facz1;
  double facx2,facy2,facz2;

  f0 = grax = gray = graz  = 0.;
  hexx = hexy = hexz = 0.;
  heyy = heyz = hezz = 0.;
  for( i = p; i >= 0; i--){
    facx0 = pow(x,i);
    facx1 = pow(x,i-1)*(double) i;
    facx2 = pow(x,i-2)*(double)(i*(i-1));

    for( j = p; j >= 0; j--){
      facy0 = pow(y,j);
      facy1 = pow(y,j-1)*(double) j;
      facy2 = pow(y,j-2)*(double)(j*(j-1));

      for( k = p; k >= 0; k--){
        c      = coef[ i*(p+1)*(p+1) + j*(p+1) + k ];

        facz0  = pow(z,k);
        facz1  = pow(z,k-1)*(double) k;
        facz2  = pow(z,k-2)*(double)(k*(k-1));

        f0   += c*facx0*facy0*facz0;

        grax += c*facx1*facy0*facz0;
        gray += c*facx0*facy1*facz0;
        graz += c*facx0*facy0*facz1;

        hexx += c*facx2*facy0*facz0;
        hexy += c*facx1*facy1*facz0;
        hexz += c*facx1*facy0*facz1;
        heyy += c*facx0*facy2*facz0;
        heyz += c*facx0*facy1*facz1;
        hezz += c*facx0*facy0*facz2;

      }  
    }
  }

  val[0] = f0;
  val[1] = grax;
  val[2] = gray;
  val[3] = graz;
  val[4] = hexx;
  val[5] = hexy;
  val[6] = hexz;
  val[7] = heyy;
  val[8] = heyz;
  val[9] = hezz;

  return 0;
}



