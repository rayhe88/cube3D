#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lectura.h"
#include "struct.h"
#include "file.h"
#include "array.h"
#include "matvec.h"

int *zatm;
double *coor;
double *field;

int main (int argc, char* argv[]){

  extern int *zatm;
  extern double *coor;
  extern double *field;

  FILE *out,*out2;

  openFile(&out,"check.cube","w+");
  dataCube cube1;

  dataCube cube2;
  double *coor2;


  loadData(&cube1,argv[1]);
  checkData(cube1.natm,cube1.pts,zatm,cube1.min,cube1.mvec,coor,
            field,"check",out);

  fclose(out);

  double matT[9];
  double matTinv[9];
  getMatT(cube1.mvec,matT,matTinv);

  createArrayDou(3*cube1.natm,&coor2,"Coor2t");

  transform(cube1,&cube2,coor,coor2,matTinv);

  openFile(&out2,"Trans.cube","w+");
  checkData(cube2.natm,cube2.pts,zatm,cube2.min,cube2.mvec,coor2,
            field,"modificado",out2);


  fclose(out2);

  free(coor2);
  free(coor);
  free(field);
  free(zatm);


  //barra(100,ND);
  exit(EXIT_SUCCESS);

}
