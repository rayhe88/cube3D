#include "file.h"
#include "array.h"
#include "struct.h"
#include "analys.h"
#include "matvec.h"
#include "lectura.h"
#include "critical.h"
#include "pruebaCoef.h"
#include "lagrangePol.h"

void printDataCube(dataCube cube);

int main(int argc, char* argv[]){

  extern int *zatm;
  extern double *coor;
  extern double *field;

  int pol;
  int task;
  int bound;
  double dg[2];
  char namefld[120];
  char nameout[120];
  char nametmp[120];

  FILE *aux;

  tmpFile(&aux,".cube",nametmp,"w+");


  dataCube data1;
  dataCube auxdata;

  readInput(namefld,nameout,&task,&pol,&bound,dg,argv[1]);

  printf("=====================================================\n");
  printf("     Name field  : %s\n",namefld);
  printf("     Name output : %s\n",nameout);
  printf("     Polynomial  : %3d\n",pol);
  printf("=====================================================\n");

  loadData(&data1,namefld);

  double matT[9],matU[9];
  double *coor2;

  getMatT(data1.mvec,matT,matU);
  createArrayDou(3*data1.natm,&coor2,"CoorT");
  transform(data1,&auxdata,coor,coor2,matU);

  printDataCube(data1);


  fclose(aux);
  remove(nametmp);

  if( task == 1 )
    newCube1(data1,zatm,coor,field,matT,nameout);
  if( task == 2 )
    newCube2(data1,zatm,coor,field,matT,nameout);
  if( task == 3 )
    newCube3(data1,zatm,coor,field,matT,dg,nameout);
  if( task == 4 )
    findCriticalPoints(pol,auxdata,zatm,coor2,field,matT,nameout);
  
  free(coor2);
  if( task == 0) {

    double *f,*c;

    double val1[10];
    double val0[10];
    double hx,hy,hz;
    double x1,y1,z1;
    int i,j,k;
    int ip,jp,kp;
    int nu;
    int npy,npz;
    int sizecf = pow(pol+1,3);

    createArrayDou(sizecf,&f,"funct");
    createArrayDou(sizecf,&c,"coeff");
    npy = data1.pts[1];
    npz = data1.pts[2];

    hx = data1.hvec[0];
    hy = data1.hvec[1];
    hz = data1.hvec[2];

    i = data1.pts[0]/2; j = data1.pts[1]/2; k = data1.pts[2]/2;
    x1  = data1.min[0] + i*hx;
    y1  = data1.min[1] + j*hy;
    z1  = data1.min[2] + k*hz;

    nu = 0;
    for(ip = i-1; ip <= i+1; ip++)
      for(jp = j-1; jp <= j+1; jp++)
        for(kp = k-1; kp <= k+1; kp++){
          f[nu]  = field[IDX(ip,jp,kp,npy,npz)];
          printf(" %4d %4d %4d | %6d ",ip,jp,kp,IDX(ip,jp,kp,npy,npz));
          printf("  f[%3d] = % 20.10lf\n",nu,f[nu]);
          nu++;
        }

  

    val0[0] = f[IDX3(1,1,1)];
    val0[1] = (f[IDX3(2,1,1)] - f[IDX3(0,1,1)])/(2.*hx);
    val0[2] = (f[IDX3(1,2,1)] - f[IDX3(1,0,1)])/(2.*hy);
    val0[3] = (f[IDX3(1,1,2)] - f[IDX3(1,1,0)])/(2.*hz);
    val0[4] = (f[IDX3(0,1,1)]- 2.*f[IDX3(1,1,1)] + f[IDX3(2,1,1)])/(hx*hx);
    val0[5] = (f[IDX3(0,0,1)] - f[IDX3(0,2,1)] -f[IDX3(2,0,1)] + f[IDX3(2,2,1)])/(4.*hx*hy);
    val0[6] = (f[IDX3(0,1,0)] - f[IDX3(0,1,2)] -f[IDX3(2,1,0)] + f[IDX3(2,1,2)])/(4.*hx*hz);
    val0[7] = (f[IDX3(1,0,1)] - 2.*f[IDX3(1,1,1)] + f[IDX3(1,2,1)])/(hy*hy);
    val0[8] = (f[IDX3(1,0,0)] - f[IDX3(1,0,2)] -f[IDX3(1,2,0)] + f[IDX3(1,2,2)])/(4.*hy*hz);
    val0[9] = (f[IDX3(1,1,0)] - 2.*f[IDX3(1,1,1)] + f[IDX3(1,1,2)])/(hz*hz);
  


    //coefficients(hx,hy,hz,x1,y1,z1,f,c);
    getCoeff(pol,data1.hvec,x1,y1,z1,f,c);
    
    //NumericalCrit02(x1,y1,z1,c,val1);
    evalPol2(pol,x1,y1,z1,c,val1);

    printf("                          Estimado0     Estimado Pol\n");
    printf(" Valor densidad :  % 10.6E   % 10.6E  % 10.6E\n",val0[0],val1[0],val0[0]-val1[0]);
    printf(" Valor gradien X:  % 10.6E   % 10.6E  % 10.6E\n",val0[1],val1[1],val0[1]-val1[1]);
    printf(" Valor gradien Y:  % 10.6E   % 10.6E  % 10.6E\n",val0[2],val1[2],val0[2]-val1[2]);
    printf(" Valor gradien Z:  % 10.6E   % 10.6E  % 10.6E\n",val0[3],val1[3],val0[3]-val1[3]);
    printf(" Valor hessia XX:  % 10.6E   % 10.6E  % 10.6E\n",val0[4],val1[4],val0[4]-val1[4]);
    printf(" Valor hessia XY:  % 10.6E   % 10.6E  % 10.6E\n",val0[5],val1[5],val0[5]-val1[5]);
    printf(" Valor hessia XZ:  % 10.6E   % 10.6E  % 10.6E\n",val0[6],val1[6],val0[6]-val1[6]);
    printf(" Valor hessia YY:  % 10.6E   % 10.6E  % 10.6E\n",val0[7],val1[7],val0[7]-val1[7]);
    printf(" Valor hessia YZ:  % 10.6E   % 10.6E  % 10.6E\n",val0[8],val1[8],val0[8]-val1[8]);
    printf(" Valor hessia ZZ:  % 10.6E   % 10.6E  % 10.6E\n",val0[9],val1[9],val0[9]-val1[9]);
  

  
    i = data1.pts[0]/2; j = data1.pts[1]/2; k = data1.pts[2]/2;
    x1  = data1.min[0] + i*hx;
    y1  = data1.min[1] + j*hy;
    z1  = data1.min[2] + k*hz;
    nu = 0;
    for(ip = i-1; ip <= i+1; ip++)
      for(jp = j-1; jp <= j+1; jp++)
        for(kp = k-1; kp <= k+1; kp++){
          f[nu]  = field[IDX(ip,jp,kp,npy,npz)];
          printf(" %4d %4d %4d | %6d ",ip,jp,kp,IDX(ip,jp,kp,npy,npz));
          printf("  f[%3d] = % 20.10lf\n",nu,f[nu]);
          nu++;
        }


    coefficients(hx,hy,hz,x1,y1,z1,f,c);

    double x,y,z,f0;
    for(ip = i-1; ip <= i+1; ip++)
      for(jp = j-1; jp <= j+1; jp++)
        for(kp = k-1; kp <= k+1; kp++){
          x = data1.min[0] + ip*hx;
          y = data1.min[1] + jp*hy;
          z = data1.min[2] + kp*hz;

          //NumericalCrit01(x,y,z,c,val1);
          evalPol1(pol,x1,y1,z1,c,val1);
          f0 = field[IDX(ip,jp,kp,npy,npz)];

          printf(" Coor( % 10.6E, % 10.6E, % 10.6E) = % 10.6E  | % 10.6E | % 10.6E\n",
                 x,y,z,f0,val1[0],fabs(f0-val1[0]));

        
        }
    
    free(c);
    free(f);

  }

   
  //free(coor);
  //free(field);
  //free(zatm);

  exit(EXIT_SUCCESS);
}


void printDataCube(dataCube cube){

  printf(" Number of atoms  : %5d\n",cube.natm);
  printf(" Total points     : %5d\n",cube.npt);
  printf(" Number of points : %5d %5d %5d\n",cube.pts[0],cube.pts[1],cube.pts[2]);
  printf(" Coordinates of r0: % 10.6lf % 10.6lf % 10.6lf\n",cube.min[0],cube.min[1],cube.min[2]);
  printf("   Hx Hy and Hz   : % 10.6lf % 10.6lf % 10.6lf\n",cube.hvec[0],cube.hvec[1],cube.hvec[2]);
  printf("                  | % 10.6lf % 10.6lf % 10.6lf|\n",cube.mvec[0],cube.mvec[1],cube.mvec[2]);
  printf("     Matriz H  =  | % 10.6lf % 10.6lf % 10.6lf|\n",cube.mvec[3],cube.mvec[4],cube.mvec[5]);
  printf("                  | % 10.6lf % 10.6lf % 10.6lf|\n",cube.mvec[6],cube.mvec[7],cube.mvec[8]);
  printf("=====================================================\n");

}
