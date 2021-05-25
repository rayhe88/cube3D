/***************************   ---- MAIN ----   ******************************/
/**
 *  @file main.c
 *  @brief Code to analize numerically scalar fields in quantum
 *         chemistry.
 *
 *  @mainpage Cube3D code v1.2
 *
 *  The main input is cube type files and using numerical interpolation
 *  reconstructs the scalar fields for later analysis.
 *
 *  The fields that it can determine are:
 *    - Gradients.
 *    - Laplacians.
 *
 *  If the field is density it can  determine
 *  the reduced gradient and
 *  the index of non-covalent interactions (NCI).
 *
 *  Additionally it allows to determine the critical points
 *  its characterization and the determination of bond paths.
 *
 *  We are working on incorporating a replicator to
 *  meshes from calculations with periodic conditions
 *  and the determination of the voids.
 *
 *  @authors Raymundo Hernández-Esparza (rayhe88@gmail.com)
 *
 *  @date August 2018
 */

#include "file.h"
#include "array.h"
#include "utils.h"
#include "struct.h"
#include "transU.h"
#include "kernels.h"
#include "lectura.h"

double getDenInCube2(int i, int j, int k, int n1, int n2, double *field);
void chargeOfSystem(dataCube cube);

int main(int argc, char* argv[]){


  char namefld[120];
  char nameout[120];
  char nametmp[120];

  int rec;
  int *zatm1;
  int rotate;

  double matT[9];
  double matU[9];
  double *coor1,*field1;

  FILE *aux;

  dataCube cube;
  dataRun  parameters;

  checkCommandLine(argc, argv);

  tmpFile(&aux,".c3dInp",nametmp,"w+");

  readInput(namefld,nameout,&parameters,argv[1]);

  printHead(argv[1],namefld,nameout);

  printRunning(parameters);

  loadData(&cube,&zatm1,&coor1,&field1,&rotate,namefld);

  printCubeRot(rotate,nameout,cube);

  printInfoCube(cube);

  fclose(aux);

  remove(nametmp);

  getMatT(cube,&rec,matT);

  parameters.orth = rec;

  getMatInv(matT,matU);
/*
  dataCube cube2;
  int *zatm2;
  double *coor2,*field2;
  FILE *out2;
  transformCube(cube,&cube2,&zatm2,&coor2,&field2,matU);
  openFile(&out2,"ortogonal.cube","w+");
  printCube("Ortogonal",cube2,out2);
  fclose(out2);
  unloadData(&cube2,&zatm2,&coor2,&field2);
 */

  //comienza la ejecución real

  //printTapas(cube);
  selectExec(cube,parameters,matU,nameout);

  unloadData(&cube,&zatm1,&coor1,&field1);

  printTime("  End  time: ");

  exit(EXIT_SUCCESS);
}
