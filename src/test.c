#include "findCrit.h"
#include "lagrange2.h"
#include "struct.h"
#include "jacobi.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

double fun(double x, double y, double z){
    return (x*x*y*z)*exp(-0.2*(x*x+y*y+z*z));
}

int main(int argc, char *argv[]){


    int i,j,k;
    double h = 0.1;
    double x0 = 1.0;
    double y0 = 1.0;
    double z0 = 0.0;
    double xx[4],yy[4],zz[4];
    double funcion[4*4*4];
    double x,y,z;

    double min2,min;
    double f;

    min = 1.E+10;

    for(i=0;i<4;i++){
        x = x0 + i*h;
        xx[i] = x;
        for(j=0;j<4;j++){
            y = y0 + j*h;
            yy[j] = y;
            for(k=0;k<4;k++){
                z = z0 + k*h;
                zz[k] = z;
                f = fun(x,y,z);
                if (f < min )
                    min = f;
                funcion[i*16+j*4+k] = f;

                printf("f(% 6.3lf,% 6.3lf,% 6.3lf) = % 10.8lf\n",x,y,z,fun(x,y,z));
            }
        }
    }


    min2 = min - DELTA;
    printf(" Minimo  de la funcion = % lf \n",min);
    printf(" Minimo2 de la funcion = % lf \n",min2);

//    for(i=0;i<64;i++){
//        funcion[i] = log(funcion[i] - min2);
 //   }

    x = 1.15;
    y = 1.15;
    z = 0.15;

    double matU[9]= { 1.0, 0.0, 0.0,
                      0.0, 1.0, 0.0,
                      0.0, 0.0, 1.0};

    double val[10];

    double cont[10]={0.1338100318,
                     0.1711604842,
                     0.05480393478,
                     0.8840382770,
                     -0.03694664455,
                     -0.1322578355,
                     -0.1600903221,
                     0.07010138092,
                     1.130800266,
                     0.3620713291};
    double cvec1[4]={-1.27283396252762445591911848233269283,
                    -0.65666034932973443449260525121964588,
                    -0.191230213093057630347980553199393638, 
                     0.72953971188588144689078259594010471};
    double cvec2[4]={-0.164109698543837608086731127468782670,
                     0.304766945161377201098435632695759818,
                    -0.95210528895711711964277206247481758, 
                     0.0247513207099462718876177482997535528};
    double cvec3[4]={ 1.10764885897261207574698334829797638,
                    -0.65666034932973443449260525121964588, 
                    -0.191230213093057630347980553199393638, 
                     0.72953971188588144689078259594010471};




    getDerivatives3D(x,y,z,xx,yy,zz,funcion,3,YES,matU,val);

    printf("      f    = % 10.8lf  %10.8lf\n",val[0],cont[0]);
    printf("    df/dx  = % 10.8lf  %10.8lf\n",val[1],cont[1]);
    printf("    df/dy  = % 10.8lf  %10.8lf\n",val[2],cont[2]);
    printf("    df/dz  = % 10.8lf  %10.8lf\n",val[3],cont[3]);
    printf("   d2f/dx2 = % 10.8lf  %10.8lf\n",val[4],cont[4]);
    printf("   d2f/dy2 = % 10.8lf  %10.8lf\n",val[5],cont[5]);
    printf("   d2f/dz2 = % 10.8lf  %10.8lf\n",val[6],cont[6]);
    printf("  d2f/dxdy = % 10.8lf  %10.8lf\n",val[7],cont[7]);
    printf("  d2f/dxdz = % 10.8lf  %10.8lf\n",val[8],cont[8]);
    printf("  d2f/dydz = % 10.8lf  %10.8lf\n",val[9],cont[9]);

    double eval[3],evec[9];
    double matH[9];

    matH[0] = val[4];  matH[1] = val[7];  matH[2] = val[8];
    matH[3] = val[7];  matH[4] = val[5];  matH[5] = val[9];
    matH[6] = val[8];  matH[7] = val[9];  matH[8] = val[6];

    JacobiNxN_2 (matH,eval,evec);

    printf(" 1(%10.8lf) % 10.8lf  % 10.8lf  % 10.8lf\n",
    cvec1[0],cvec1[1],cvec1[2],cvec1[3]);
    printf(" 2(%10.8lf) % 10.8lf  % 10.8lf  % 10.8lf\n",
    cvec2[0],cvec2[1],cvec2[2],cvec2[3]);
    printf(" 3(%10.8lf) % 10.8lf  % 10.8lf  % 10.8lf\n",
    cvec3[0],cvec3[1],cvec3[2],cvec3[3]);


    printf(" 1(%10.8lf) % 10.8lf  % 10.8lf  % 10.8lf\n",
    eval[0],evec[0],evec[1],evec[2]);
    printf(" 2(%10.8lf) % 10.8lf  % 10.8lf  % 10.8lf\n",
    eval[1],evec[3],evec[4],evec[5]);
    printf(" 3(%10.8lf) % 10.8lf  % 10.8lf  % 10.8lf\n",
    eval[2],evec[6],evec[7],evec[8]);

    for(i=0;i<64;i++){
        funcion[i] = log(funcion[i] - min2);
   }
    getDerivatives3DLog(x,y,z,xx,yy,zz,funcion,3,YES,matU,min,val);
    printf("      f    = % 10.8lf  %10.8lf\n",val[0],cont[0]);
    printf("    df/dx  = % 10.8lf  %10.8lf\n",val[1],cont[1]);
    printf("    df/dy  = % 10.8lf  %10.8lf\n",val[2],cont[2]);
    printf("    df/dz  = % 10.8lf  %10.8lf\n",val[3],cont[3]);
    printf("   d2f/dx2 = % 10.8lf  %10.8lf\n",val[4],cont[4]);
    printf("   d2f/dy2 = % 10.8lf  %10.8lf\n",val[5],cont[5]);
    printf("   d2f/dz2 = % 10.8lf  %10.8lf\n",val[6],cont[6]);
    printf("  d2f/dxdy = % 10.8lf  %10.8lf\n",val[7],cont[7]);
    printf("  d2f/dxdz = % 10.8lf  %10.8lf\n",val[8],cont[8]);
    printf("  d2f/dydz = % 10.8lf  %10.8lf\n",val[9],cont[9]);

    return 0;

}
