#include "lebedev.h"

void translateAndChange(int num, dataVec *leb, double r0[3], double fac){

    int i;
    for(i=0;i++;i++){
        leb[i].x *= fac;
        leb[i].y *= fac;
        leb[i].z *= fac;
        leb[i].x += r0[0];
        leb[i].y += r0[1];
        leb[i].z += r0[2];
    }
}

int createArrayVec( int num, dataVec **ptr, const char *mess){

  (*ptr) =(dataVec*) malloc (num*sizeof(dataVec));
  if( (*ptr) == NULL ){
    printf(" Error al reservar memoria para [%s]\n",mess);
    exit(EXIT_FAILURE);
  }

  return 0;
}

void genGroup1 ( int *num, dataVec *r, double va, double vb, double v){
  int i = (*num);
  r[i  ].x =  1.; r[i  ].y =  0.; r[i  ].z =  0.;  r[i  ].w = v;
  r[i+1].x = -1.; r[i+1].y =  0.; r[i+1].z =  0.;  r[i+1].w = v;
  r[i+2].x =  0.; r[i+2].y =  1.; r[i+2].z =  0.;  r[i+2].w = v;
  r[i+3].x =  0.; r[i+3].y = -1.; r[i+3].z =  0.;  r[i+3].w = v;
  r[i+4].x =  0.; r[i+4].y =  0.; r[i+4].z =  1.;  r[i+4].w = v;
  r[i+5].x =  0.; r[i+5].y =  0.; r[i+5].z = -1.;  r[i+5].w = v;
  
  (*num) += 6;

}
void genGroup2 ( int *num, dataVec *r, double va, double vb, double v){
  int i = (*num);
  double a = sqrt(0.5);
  r[i   ].x = 0.; r[i   ].y =  a; r[i   ].z =  a;  r[i   ].w = v;
  r[i+ 1].x = 0.; r[i+ 1].y =  a; r[i+ 1].z = -a;  r[i+ 1].w = v;
  r[i+ 2].x = 0.; r[i+ 2].y = -a; r[i+ 2].z =  a;  r[i+ 2].w = v;
  r[i+ 3].x = 0.; r[i+ 3].y = -a; r[i+ 3].z = -a;  r[i+ 3].w = v;
  r[i+ 4].x =  a; r[i+ 4].y = 0.; r[i+ 4].z =  a;  r[i+ 4].w = v;
  r[i+ 5].x =  a; r[i+ 5].y = 0.; r[i+ 5].z = -a;  r[i+ 5].w = v;
  r[i+ 6].x = -a; r[i+ 6].y = 0.; r[i+ 6].z =  a;  r[i+ 6].w = v;
  r[i+ 7].x = -a; r[i+ 7].y = 0.; r[i+ 7].z = -a;  r[i+ 7].w = v;
  r[i+ 8].x =  a; r[i+ 8].y =  a; r[i+ 8].z = 0.;  r[i+ 8].w = v;
  r[i+ 9].x =  a; r[i+ 9].y = -a; r[i+ 9].z = 0.;  r[i+ 9].w = v;
  r[i+10].x = -a; r[i+10].y =  a; r[i+10].z = 0.;  r[i+10].w = v;
  r[i+11].x = -a; r[i+11].y = -a; r[i+11].z = 0.;  r[i+11].w = v;

  (*num) += 12;
};
void genGroup3 ( int *num, dataVec *r, double va, double vb, double v){
  int i = (*num);
  double a = 1./sqrt(3.);
  r[i   ].x =  a; r[i   ].y =  a; r[i   ].z =  a;  r[i   ].w = v;
  r[i+ 1].x =  a; r[i+ 1].y =  a; r[i+ 1].z = -a;  r[i+ 1].w = v;
  r[i+ 2].x =  a; r[i+ 2].y = -a; r[i+ 2].z =  a;  r[i+ 2].w = v;
  r[i+ 3].x =  a; r[i+ 3].y = -a; r[i+ 3].z = -a;  r[i+ 3].w = v;
  r[i+ 4].x = -a; r[i+ 4].y =  a; r[i+ 4].z =  a;  r[i+ 4].w = v;
  r[i+ 5].x = -a; r[i+ 5].y =  a; r[i+ 5].z = -a;  r[i+ 5].w = v;
  r[i+ 6].x = -a; r[i+ 6].y = -a; r[i+ 6].z =  a;  r[i+ 6].w = v;
  r[i+ 7].x = -a; r[i+ 7].y = -a; r[i+ 7].z = -a;  r[i+ 7].w = v;

  (*num) += 8;
};
void genGroup4 ( int *num, dataVec *r, double va, double vb, double v){
  int i = (*num);
  double a = va;
  double b = sqrt( 1. - 2.*va*va);
  r[i   ].x =  a; r[i   ].y =  a; r[i   ].z =  b;  r[i   ].w = v;
  r[i+ 1].x = -a; r[i+ 1].y =  a; r[i+ 1].z =  b;  r[i+ 1].w = v;
  r[i+ 2].x =  a; r[i+ 2].y = -a; r[i+ 2].z =  b;  r[i+ 2].w = v;
  r[i+ 3].x = -a; r[i+ 3].y = -a; r[i+ 3].z =  b;  r[i+ 3].w = v;
  r[i+ 4].x =  a; r[i+ 4].y =  a; r[i+ 4].z = -b;  r[i+ 4].w = v;
  r[i+ 5].x = -a; r[i+ 5].y =  a; r[i+ 5].z = -b;  r[i+ 5].w = v;
  r[i+ 6].x =  a; r[i+ 6].y = -a; r[i+ 6].z = -b;  r[i+ 6].w = v;
  r[i+ 7].x = -a; r[i+ 7].y = -a; r[i+ 7].z = -b;  r[i+ 7].w = v;
  r[i+ 8].x =  a; r[i+ 8].y =  b; r[i+ 8].z =  a;  r[i+ 8].w = v;
  r[i+ 9].x = -a; r[i+ 9].y =  b; r[i+ 9].z =  a;  r[i+ 9].w = v;
  r[i+10].x =  a; r[i+10].y =  b; r[i+10].z = -a;  r[i+10].w = v;
  r[i+11].x = -a; r[i+11].y =  b; r[i+11].z = -a;  r[i+11].w = v;
  r[i+12].x =  a; r[i+12].y = -b; r[i+12].z =  a;  r[i+12].w = v;
  r[i+13].x = -a; r[i+13].y = -b; r[i+13].z =  a;  r[i+13].w = v;
  r[i+14].x =  a; r[i+14].y = -b; r[i+14].z = -a;  r[i+14].w = v;
  r[i+15].x = -a; r[i+15].y = -b; r[i+15].z = -a;  r[i+15].w = v;
  r[i+16].x =  b; r[i+16].y =  a; r[i+16].z =  a;  r[i+16].w = v;
  r[i+17].x =  b; r[i+17].y = -a; r[i+17].z =  a;  r[i+17].w = v;
  r[i+18].x =  b; r[i+18].y =  a; r[i+18].z = -a;  r[i+18].w = v;
  r[i+19].x =  b; r[i+19].y = -a; r[i+19].z = -a;  r[i+19].w = v;
  r[i+20].x = -b; r[i+20].y =  a; r[i+20].z =  a;  r[i+20].w = v;
  r[i+21].x = -b; r[i+21].y = -a; r[i+21].z =  a;  r[i+21].w = v;
  r[i+22].x = -b; r[i+22].y =  a; r[i+22].z = -a;  r[i+22].w = v;
  r[i+23].x = -b; r[i+23].y = -a; r[i+23].z = -a;  r[i+23].w = v;

  (*num) += 24;
};
void genGroup5 ( int *num, dataVec *r, double va, double vb, double v){
  int i = (*num);
  double a = va;
  double b = sqrt( 1. - va*va);
  r[i   ].x =  a; r[i   ].y =  b; r[i   ].z = 0.;  r[i   ].w = v;
  r[i+ 1].x =  a; r[i+ 1].y = -b; r[i+ 1].z = 0.;  r[i+ 1].w = v;
  r[i+ 2].x = -a; r[i+ 2].y =  b; r[i+ 2].z = 0.;  r[i+ 2].w = v;
  r[i+ 3].x = -a; r[i+ 3].y = -b; r[i+ 3].z = 0.;  r[i+ 3].w = v;
  r[i+ 4].x =  b; r[i+ 4].y =  a; r[i+ 4].z = 0.;  r[i+ 4].w = v;
  r[i+ 5].x = -b; r[i+ 5].y =  a; r[i+ 5].z = 0.;  r[i+ 5].w = v;
  r[i+ 6].x =  b; r[i+ 6].y = -a; r[i+ 6].z = 0.;  r[i+ 6].w = v;
  r[i+ 7].x = -b; r[i+ 7].y = -a; r[i+ 7].z = 0.;  r[i+ 7].w = v;
  r[i+ 8].x =  a; r[i+ 8].y = 0.; r[i+ 8].z =  b;  r[i+ 8].w = v;
  r[i+ 9].x =  a; r[i+ 9].y = 0.; r[i+ 9].z = -b;  r[i+ 9].w = v;
  r[i+10].x = -a; r[i+10].y = 0.; r[i+10].z =  b;  r[i+10].w = v;
  r[i+11].x = -a; r[i+11].y = 0.; r[i+11].z = -b;  r[i+11].w = v;
  r[i+12].x =  b; r[i+12].y = 0.; r[i+12].z =  a;  r[i+12].w = v;
  r[i+13].x = -b; r[i+13].y = 0.; r[i+13].z =  a;  r[i+13].w = v;
  r[i+14].x =  b; r[i+14].y = 0.; r[i+14].z = -a;  r[i+14].w = v;
  r[i+15].x = -b; r[i+15].y = 0.; r[i+15].z = -a;  r[i+15].w = v;
  r[i+16].x = 0.; r[i+16].y =  a; r[i+16].z =  b;  r[i+16].w = v;
  r[i+17].x = 0.; r[i+17].y =  a; r[i+17].z = -b;  r[i+17].w = v;
  r[i+18].x = 0.; r[i+18].y = -a; r[i+18].z =  b;  r[i+18].w = v;
  r[i+19].x = 0.; r[i+19].y = -a; r[i+19].z = -b;  r[i+19].w = v;
  r[i+20].x = 0.; r[i+20].y =  b; r[i+20].z =  a;  r[i+20].w = v;
  r[i+21].x = 0.; r[i+21].y = -b; r[i+21].z =  a;  r[i+21].w = v;
  r[i+22].x = 0.; r[i+22].y =  b; r[i+22].z = -a;  r[i+22].w = v;
  r[i+23].x = 0.; r[i+23].y = -b; r[i+23].z = -a;  r[i+23].w = v;

  (*num) += 24;
};
void genGroup6 ( int *num, dataVec *r, double va, double vb, double v){
  int i = (*num);
  double a = va;
  double b = vb;
  double c = sqrt(1. - va*va - vb*vb);

  r[i   ].x =  a; r[i   ].y =  b; r[i   ].z =  c;  r[i   ].w = v;
  r[i+ 1].x =  a; r[i+ 1].y = -b; r[i+ 1].z =  c;  r[i+ 1].w = v;
  r[i+ 2].x = -a; r[i+ 2].y =  b; r[i+ 2].z =  c;  r[i+ 2].w = v;
  r[i+ 3].x = -a; r[i+ 3].y = -b; r[i+ 3].z =  c;  r[i+ 3].w = v;
  r[i+ 4].x =  b; r[i+ 4].y =  a; r[i+ 4].z =  c;  r[i+ 4].w = v;
  r[i+ 5].x = -b; r[i+ 5].y =  a; r[i+ 5].z =  c;  r[i+ 5].w = v;
  r[i+ 6].x =  b; r[i+ 6].y = -a; r[i+ 6].z =  c;  r[i+ 6].w = v;
  r[i+ 7].x = -b; r[i+ 7].y = -a; r[i+ 7].z =  c;  r[i+ 7].w = v;
  r[i+ 8].x =  a; r[i+ 8].y =  b; r[i+ 8].z = -c;  r[i+ 8].w = v;
  r[i+ 9].x =  a; r[i+ 9].y = -b; r[i+ 9].z = -c;  r[i+ 9].w = v;
  r[i+10].x = -a; r[i+10].y =  b; r[i+10].z = -c;  r[i+10].w = v;
  r[i+11].x = -a; r[i+11].y = -b; r[i+11].z = -c;  r[i+11].w = v;
  r[i+12].x =  b; r[i+12].y =  a; r[i+12].z = -c;  r[i+12].w = v;
  r[i+13].x = -b; r[i+13].y =  a; r[i+13].z = -c;  r[i+13].w = v;
  r[i+14].x =  b; r[i+14].y = -a; r[i+14].z = -c;  r[i+14].w = v;
  r[i+15].x = -b; r[i+15].y = -a; r[i+15].z = -c;  r[i+15].w = v;
  r[i+16].x =  a; r[i+16].y =  c; r[i+16].z =  b;  r[i+16].w = v;
  r[i+17].x =  a; r[i+17].y =  c; r[i+17].z = -b;  r[i+17].w = v;
  r[i+18].x = -a; r[i+18].y =  c; r[i+18].z =  b;  r[i+18].w = v;
  r[i+19].x = -a; r[i+19].y =  c; r[i+19].z = -b;  r[i+19].w = v;
  r[i+20].x =  b; r[i+20].y =  c; r[i+20].z =  a;  r[i+20].w = v;
  r[i+21].x = -b; r[i+21].y =  c; r[i+21].z =  a;  r[i+21].w = v;
  r[i+22].x =  b; r[i+22].y =  c; r[i+22].z = -a;  r[i+22].w = v;
  r[i+23].x = -b; r[i+23].y =  c; r[i+23].z = -a;  r[i+23].w = v;
  r[i+24].x =  a; r[i+24].y = -c; r[i+24].z =  b;  r[i+24].w = v;
  r[i+25].x =  a; r[i+25].y = -c; r[i+25].z = -b;  r[i+25].w = v;
  r[i+26].x = -a; r[i+26].y = -c; r[i+26].z =  b;  r[i+26].w = v;
  r[i+27].x = -a; r[i+27].y = -c; r[i+27].z = -b;  r[i+27].w = v;
  r[i+28].x =  b; r[i+28].y = -c; r[i+28].z =  a;  r[i+28].w = v;
  r[i+29].x = -b; r[i+29].y = -c; r[i+29].z =  a;  r[i+29].w = v;
  r[i+30].x =  b; r[i+30].y = -c; r[i+30].z = -a;  r[i+30].w = v;
  r[i+31].x = -b; r[i+31].y = -c; r[i+31].z = -a;  r[i+31].w = v;
  r[i+32].x =  c; r[i+32].y =  a; r[i+32].z =  b;  r[i+32].w = v;
  r[i+33].x =  c; r[i+33].y =  a; r[i+33].z = -b;  r[i+33].w = v;
  r[i+34].x =  c; r[i+34].y = -a; r[i+34].z =  b;  r[i+34].w = v;
  r[i+35].x =  c; r[i+35].y = -a; r[i+35].z = -b;  r[i+35].w = v;
  r[i+36].x =  c; r[i+36].y =  b; r[i+36].z =  a;  r[i+36].w = v;
  r[i+37].x =  c; r[i+37].y = -b; r[i+37].z =  a;  r[i+37].w = v;
  r[i+38].x =  c; r[i+38].y =  b; r[i+38].z = -a;  r[i+38].w = v;
  r[i+39].x =  c; r[i+39].y = -b; r[i+39].z = -a;  r[i+39].w = v;
  r[i+40].x = -c; r[i+40].y =  a; r[i+40].z =  b;  r[i+40].w = v;
  r[i+41].x = -c; r[i+41].y =  a; r[i+41].z = -b;  r[i+41].w = v;
  r[i+42].x = -c; r[i+42].y = -a; r[i+42].z =  b;  r[i+42].w = v;
  r[i+43].x = -c; r[i+43].y = -a; r[i+43].z = -b;  r[i+43].w = v;
  r[i+44].x = -c; r[i+44].y =  b; r[i+44].z =  a;  r[i+44].w = v;
  r[i+45].x = -c; r[i+45].y = -b; r[i+45].z =  a;  r[i+45].w = v;
  r[i+46].x = -c; r[i+46].y =  b; r[i+46].z = -a;  r[i+46].w = v;
  r[i+47].x = -c; r[i+47].y = -b; r[i+47].z = -a;  r[i+47].w = v;

  (*num) += 48;   
};

int Lebedev0006 ( dataVec *r ){
  int num = 0;
  double a=0.,b=0.,v=0.;
  v = 0.1666666666666667;
  genGroup1 ( &num,r,a,b,v);

  return 0;
}

int Lebedev0014 ( dataVec *r ){
  int num = 0;
  double a=0.,b=0.,v=0.;
  v = 0.6666666666666667E-1;
  genGroup1 ( &num,r,a,b,v);
  v = 0.7500000000000000E-1;
  genGroup3 ( &num,r,a,b,v);

  return 0;
}

int Lebedev0026 ( dataVec *r ){
  int num = 0;
  double a=0.,b=0.,v=0.;
  v = 0.4761904761904762E-1; 
  genGroup1 ( &num,r,a,b,v);
  v = 0.3809523809523810E-1;
  genGroup2 ( &num,r,a,b,v);
  v = 0.3214285714285714E-1;
  genGroup3 ( &num,r,a,b,v);

  return 0;
}


int Lebedev0038 ( dataVec *r ){
  int num = 0;
  double a=0.,b=0.,v=0.;
  v = 0.9523809523809524E-2;
  genGroup1 ( &num,r,a,b,v);
  v = 0.3214285714285714E-1;
  genGroup3 ( &num,r,a,b,v);
  a = 0.4597008433809831;
  v = 0.2857142857142857E-1;
  genGroup5 ( &num,r,a,b,v);

  return 0;
}

int Lebedev0050 ( dataVec *r ){
  int num = 0;
  double a=0.,b=0.,v=0.;
  v = 0.1269841269841270E-1;
  genGroup1 ( &num,r,a,b,v);
  v = 0.2257495590828924E-1;
  genGroup2 ( &num,r,a,b,v);
  v = 0.2109375000000000E-1;
  genGroup3 ( &num,r,a,b,v);
  a = 0.3015113445777636;
  v = 0.2017333553791887E-1;
  genGroup4 ( &num,r,a,b,v);

  return 0;
}

int Lebedev0074 ( dataVec *r ){
  int num = 0;
  double a=0.,b=0.,v=0.;
  v = 0.5130671797338464E-3;
  genGroup1 ( &num,r,a,b,v);
  v = 0.1660406956574204E-1;
  genGroup2 ( &num,r,a,b,v);
  v = -0.2958603896103896E-1;
  genGroup3 ( &num,r,a,b,v);
  a = 0.4803844614152614;
  v = 0.2657620708215946E-1;
  genGroup4 ( &num,r,a,b,v);
  a = 0.3207726489807764;
  v = 0.1652217099371571E-1;
  genGroup5 ( &num,r,a,b,v);

  return 0;
}

int Lebedev0086 ( dataVec *r ){
  int num = 0;
  double a=0.,b=0.,v=0.;
  v = 0.1154401154401154E-1;
  genGroup1 ( &num,r,a,b,v);
  v = 0.1194390908585628E-1;
  genGroup3 ( &num,r,a,b,v);
  a = 0.3696028464541502;
  v = 0.1111055571060340E-1;
  genGroup4 ( &num,r,a,b,v);
  a = 0.6943540066026664;
  v = 0.1187650129453714E-1;
  genGroup4 ( &num,r,a,b,v);
  a = 0.3742430390903412;
  v = 0.1181230374690448E-1;
  genGroup5 ( &num,r,a,b,v);

  return 0;
}

int Lebedev0110 ( dataVec *r ){
  int num = 0;
  double a=0.,b=0.,v=0.;
  v = 0.3828270494937162E-2; 
  genGroup1 ( &num,r,a,b,v);
  v = 0.9793737512487512E-2;
  genGroup3 ( &num,r,a,b,v);
  a = 0.1851156353447362;
  v = 0.8211737283191111E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.6904210483822922;
  v = 0.9942814891178103E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.3956894730559419;
  v = 0.9595471336070963E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.4783690288121502;
  v = 0.9694996361663028E-2;
  genGroup5 ( &num,r,a,b,v);

  return 0;
}

int Lebedev0146 ( dataVec *r ){
  int num = 0;
  double a=0.,b=0.,v=0.;
  v = 0.5996313688621381E-3;
  genGroup1 ( &num,r,a,b,v);
  v = 0.7372999718620756E-2;
  genGroup2 ( &num,r,a,b,v);
  v = 0.7210515360144488E-2;
  genGroup3 ( &num,r,a,b,v);
  a = 0.6764410400114264;
  v = 0.7116355493117555E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.4174961227965453;
  v = 0.6753829486314477E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.1574676672039082;
  v = 0.7574394159054034E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.1403553811713183;
  b = 0.4493328323269557;
  v = 0.6991087353303262E-2;
  genGroup6 ( &num,r,a,b,v);

  return 0;
}

int Lebedev0170 ( dataVec *r ){
  int num = 0;
  double a=0.,b=0.,v=0.;
  v = 0.5544842902037365E-2;
  genGroup1 ( &num,r,a,b,v);
  v = 0.6071332770670752E-2;
  genGroup2 ( &num,r,a,b,v);
  v = 0.6383674773515093E-2;
  genGroup3 ( &num,r,a,b,v);
  a = 0.2551252621114134;
  v = 0.5183387587747790E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.6743601460362766;
  v = 0.6317929009813725E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.4318910696719410;
  v = 0.6201670006589077E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.2613931360335988;
  v = 0.5477143385137348E-2;
  genGroup5 ( &num,r,a,b,v);
  a = 0.4990453161796037;
  b = 0.1446630744325115;
  v = 0.5968383987681156E-2;
  genGroup6 ( &num,r,a,b,v);

  return 0;
}

int Lebedev0194 ( dataVec *r ){
  int num = 0;
  double a=0.,b=0.,v=0.;
  v = 0.1782340447244611E-2;
  genGroup1 ( &num,r,a,b,v);
  v = 0.5716905949977102E-2;
  genGroup2 ( &num,r,a,b,v);
  v = 0.5573383178848738E-2;
  genGroup3 ( &num,r,a,b,v);
  a = 0.6712973442695226;
  v = 0.5608704082587997E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.2892465627575439;
  v = 0.5158237711805383E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.4446933178717437;
  v = 0.5518771467273614E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.1299335447650067;
  v = 0.4106777028169394E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.3457702197611283;
  v = 0.5051846064614808E-2;
  genGroup5 ( &num,r,a,b,v);
  a = 0.1590417105383530;
  b = 0.8360360154824589;
  v = 0.5530248916233094E-2;
  genGroup6 ( &num,r,a,b,v);

  return 0;
}


int Lebedev0230 ( dataVec *r ){
  int num = 0;
  double a=0.,b=0.,v=0.;
  v = -0.5522639919727325E-1;
  genGroup1 ( &num,r,a,b,v);
  v = 0.4450274607445226E-2;
  genGroup3 ( &num,r,a,b,v);
  a = 0.4492044687397611;
  v = 0.4496841067921404E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.2520419490210201;
  v = 0.5049153450478750E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.6981906658447242;
  v = 0.3976408018051883E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.6587405243460960;
  v = 0.4401400650381014E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.4038544050097660E-1; 
  v = 0.1724544350544401E-1;
  genGroup4 ( &num,r,a,b,v);
  a = 0.5823842309715585; 
  v = 0.4231083095357343E-2;
  genGroup5 ( &num,r,a,b,v);
  a = 0.3545877390518688;
  v = 0.5198069864064399E-2;
  genGroup5 ( &num,r,a,b,v);
  a = 0.2272181808998187;
  b = 0.4864661535886647;
  v = 0.4695720972568883E-2; 
  genGroup6 ( &num,r,a,b,v);

  return 0;
}

int Lebedev0266 ( dataVec *r ){
  int num = 0;
  double a=0.,b=0.,v=0.;
  v = -0.1313769127326952E-2;
  genGroup1 ( &num,r,a,b,v);
  v = -0.2522728704859336E-2;
  genGroup2 ( &num,r,a,b,v);
  v = 0.4186853881700583E-2;
  genGroup3 ( &num,r,a,b,v);
  a = 0.7039373391585475;
  v = 0.5315167977810885E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.1012526248572414;
  v = 0.4047142377086219E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.4647448726420539;
  v = 0.4112482394406990E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.3277420654971629;
  v = 0.3595584899758782E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.6620338663699974;
  v = 0.4256131351428158E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.8506508083520399;
  v = 0.4229582700647240E-2;
  genGroup5 ( &num,r,a,b,v);
  a = 0.3233484542692899;
  b = 0.1153112011009701;
  v = 0.4080914225780505E-2;
  genGroup6 ( &num,r,a,b,v);
  a = 0.2314790158712601;
  b = 0.5244939240922365;
  v = 0.4071467593830964E-2;
  genGroup6 ( &num,r,a,b,v);

  return 0;
}


int Lebedev0302 ( dataVec *r ){
  int num = 0;
  double a=0.,b=0.,v=0.;
  v = 0.8545911725128148E-3;
  genGroup1 ( &num,r,a,b,v);
  v = 0.3599119285025571E-2;
  genGroup3 ( &num,r,a,b,v);
  a = 0.3515640345570105;
  v = 0.3449788424305883E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.6566329410219612;
  v = 0.3604822601419882E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.4729054132581005;
  v = 0.3576729661743367E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.9618308522614784E-1;
  v = 0.2352101413689164E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.2219645236294178;
  v = 0.3108953122413675E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.7011766416089545;
  v = 0.3650045807677255E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.2644152887060663;
  v = 0.2982344963171804E-2;
  genGroup5 ( &num,r,a,b,v);
  a = 0.5718955891878961;
  v = 0.3600820932216460E-2;
  genGroup5 ( &num,r,a,b,v);
  a = 0.2510034751770465;
  b = 0.8000727494073952;
  v = 0.3571540554273387E-2;
  genGroup6 ( &num,r,a,b,v);
  a = 0.1233548532583327;
  b = 0.4127724083168531;
  v = 0.3392312205006170E-2;
  genGroup6 ( &num,r,a,b,v);

  return 0;
}


int Lebedev0350 ( dataVec *r ){
  int num = 0;
  double a=0.,b=0.,v=0.;
  v = 0.3006796749453936E-2;
  genGroup1 ( &num,r,a,b,v);
  v = 0.3050627745650771E-2;
  genGroup3 ( &num,r,a,b,v);
  a = 0.7068965463912316;
  v = 0.1621104600288991E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.4794682625712025;
  v = 0.3005701484901752E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.1927533154878019;
  v = 0.2990992529653774E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.6930357961327123;
  v = 0.2982170644107595E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.3608302115520091;
  v = 0.2721564237310992E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.6498486161496169;
  v = 0.3033513795811141E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.1932945013230339;
  v = 0.3007949555218533E-2;
  genGroup5 ( &num,r,a,b,v);
  a = 0.3800494919899303;
  v = 0.2881964603055307E-2;
  genGroup5 ( &num,r,a,b,v);
  a = 0.2899558825499574;
  b = 0.7934537856582316;
  v = 0.2958357626535696E-2;
  genGroup6 ( &num,r,a,b,v);
  a = 0.9684121455103957E-1;
  b = 0.8280801506686862;
  v = 0.3036020026407088E-2;
  genGroup6 ( &num,r,a,b,v);
  a = 0.1833434647041659;
  b = 0.9074658265305127;
  v = 0.2832187403926303E-2;
  genGroup6 ( &num,r,a,b,v);

  return 0;
}

int Lebedev0434 ( dataVec *r ){
  int num = 0;
  double a=0.,b=0.,v=0.;
  v = 0.5265897968224436E-3;
  genGroup1 ( &num,r,a,b,v);
  v = 0.2548219972002607E-2;
  genGroup2 ( &num,r,a,b,v);
  v = 0.2512317418927307E-2;
  genGroup3 ( &num,r,a,b,v);
  a = 0.6909346307509111;
  v = 0.2530403801186355E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.1774836054609158;
  v = 0.2014279020918528E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.4914342637784746;
  v = 0.2501725168402936E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.6456664707424256;
  v = 0.2513267174597564E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.2861289010307638;
  v = 0.2302694782227416E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.7568084367178018E-1;
  v = 0.1462495621594614E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.3927259763368002;
  v = 0.2445373437312980E-2;
  genGroup4 ( &num,r,a,b,v);
  a = 0.8818132877794288;
  v = 0.2417442375638981E-2;
  genGroup5 ( &num,r,a,b,v);
  a = 0.9776428111182649;
  v = 0.1910951282179532E-2;
  genGroup5 ( &num,r,a,b,v);
  a = 0.2054823696403044;
  b = 0.8689460322872412;
  v = 0.2416930044324775E-2;
  genGroup6 ( &num,r,a,b,v);
  a = 0.5905157048925271;
  b = 0.7999278543857286;
  v = 0.2512236854563495E-2;
  genGroup6 ( &num,r,a,b,v);
  a = 0.5550152361076807;
  b = 0.7717462626915901;
  v = 0.2496644054553086E-2;
  genGroup6 ( &num,r,a,b,v);
  a = 0.9371809858553722;
  b = 0.3344363145343455;
  v = 0.2236607760437849E-2;
  genGroup6 ( &num,r,a,b,v);

  return 0;
}
