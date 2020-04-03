#include "utils.h"

void printSizeCube(dataCube cube){

  int memory;
  char ascmem[20];

  // Al cargar el cube file se cargan
  // 5 valores enteros por default
  // y 15 valores tipo double
  memory  =  5*sizeof (int);
  memory += 15*sizeof (double);
  // Se cargan dependiendo del tamagno del cube:
  //   natm (   int) para los numeros atómicos,
  // 3*natm (double) para las coordenadas,
  //   npt  (double) para el campo escalar.
  memory +=   cube.natm*sizeof(int);
  memory += 3*cube.natm*sizeof(double);
  memory +=   cube.npt *sizeof(double);

  sizeUnits(memory,ascmem);
  
  printBar(stdout);
  printBanner(" Memory Info ",stdout);
  printf(" Size of  int         : %10ld (B)\n",sizeof(int));
  printf(" Size of double       : %10ld (B)\n",sizeof(double));
  printf(" Information loaded from cube file : %-20s\n",ascmem);

}

void fieldMinMax(dataCube cube, double *min, double *max){
  unsigned int i;
  double val,fmin,fmax;
  fmin = 1.E10;
  fmax = -fmin;
  
  for(i=0;i<cube.npt;i++){
    val = cube.field[i];
    if( val < fmin )
      fmin = val;
    if( val > fmax )
      fmax = val;
  }

  (*min) = fmin;
  (*max) = fmax;

}

void cpyDataCube(dataCube cIn,dataCube *cOut){

  int i;
  cOut -> natm = cIn.natm;
  cOut -> npt  = cIn.npt;
  cOut -> zatm = cIn.zatm;
  cOut -> coor = cIn.coor;
  cOut -> field = cIn.field;

  for(i=0;i<3;i++){
    cOut -> pts[i]  = cIn.pts[i];
    cOut -> min[i]  = cIn.min[i];
    cOut -> hvec[i] = cIn.hvec[i];
  }

  for(i=0;i<9;i++)
    cOut -> mvec[i] = cIn.mvec[i];

}

void printBar   (FILE *out){
  int i;
  for(i=0;i<NCHAR;i++)
    fprintf(out,"=");
  fprintf(out,"\n");
}

void printBanner (char *mess, FILE *out){

  int i;
  int nd,nz;
  int nc; 
  int nco2;

  char c = ' ';

  nc = strlen(mess);

  if((nc-4) > NCHAR )
    fprintf(out,"%s\n",mess);

  if( nc%2 == 0 ) 
    nco2 = nc/2;
  else{
    nco2 = (nc+1)/2;
    fprintf(out,"%c",c);
  }

  nz = (NCHAR/2) - nco2 -2; 
  nd = (NCHAR/2) - nco2 -2; 

  for(i=0;i<nz;i++)
    fprintf(out,"%c",c);
  fprintf(out,"  %s  ",mess);
  for(i=0;i<nd;i++)
    fprintf(out,"%c",c);
  fprintf(out,"\n");

  printBar(out);
}

void printBar82   (FILE *out){
  int i;
  for(i=0;i<82;i++)
    fprintf(out,"=");
  fprintf(out,"\n");
}

void printBanner82 (char *mess, FILE *out){

  int i;
  int nd,nz;
  int nc; 
  int nco2;

  char c = ' ';

  nc = strlen(mess);

  if((nc-4) > 82 )
    fprintf(out,"%s\n",mess);

  if( nc%2 == 0 ) 
    nco2 = nc/2;
  else{
    nco2 = (nc+1)/2;
    fprintf(out,"%c",c);
  }

  nz = (41) - nco2 -2; 
  nd = (41) - nco2 -2; 

  for(i=0;i<nz;i++)
    fprintf(out,"%c",c);
  fprintf(out,"  %s  ",mess);
  for(i=0;i<nd;i++)
    fprintf(out,"%c",c);
  fprintf(out,"\n");

  printBar82(out);
}


void printTime(char *text){
  time_t timeNow;
  struct tm *timeinfo;
  time (&timeNow);
  timeinfo = localtime(&timeNow);
  printf(" %s :%45s",text,asctime(timeinfo));
  printBar(stdout);
}

void printHead (char *namei, char *namef, char *nameout){

  printBar(stdout);
  printf(" Compilation date :   %-18s  %18s\n",__DATE__,__TIME__);
  printBar(stdout);
  printf("               ____      _          _____ ____             \n");
  printf("              / ___|   _| |__   ___|___ /|  _ \\           \n");
  printf("             | |  | | | | '_ \\ / _ \\ |_ \\| | | |        \n");
  printf("             | |__| |_| | |_) |  __/___) | |_| |           \n");
  printf("              \\____\\__,_|_.__/ \\___|____/|____/         \n");
  printf("                                                           \n");
  printBar(stdout);
  printTime(" Start time: ");
  printBanner(" Input Info ",stdout);
  printf(" Input's name         : %s\n",namei);
  printf(" Field file's name    : %s\n",namef);
  printf(" Output's name        : %s\n",nameout);
  printBar(stdout);
}

void printRunning (dataRun param){
  int task = param.task;
  char ascTask[6];
  char ascBol[2][4]={"NOT","YES"};

  switch( task ){
    case RED: strcpy(ascTask,"  RED"); break;
    case GRA: strcpy(ascTask," GRAD"); break;
    case LAP: strcpy(ascTask,"  LAP"); break;
    case KIN: strcpy(ascTask,"  KIN"); break;
    case VIR: strcpy(ascTask,"  VIR"); break;
    case NCI: strcpy(ascTask,"  NCI"); break;
    case CRI: strcpy(ascTask," CRIT"); break;
    case VOI: strcpy(ascTask,"VOIDS"); break;
    case REP: strcpy(ascTask,"REPLI"); break;
  }

  printBanner(" Running Info ",stdout);
  printf(" Task to be performed : %10s\n",ascTask);
  if ( task == NCI ){
    printf(" Cutoff for density   : % 10.6lf\n",param.la2);
    printf(" Cutoff for reduced   : % 10.6lf\n",param.rgd);
  }

  if( task >= RED && task <= CRI ){
    printf(" Polynomial degree    : %10d\n",param.pol);
    printf(" Periodic boundary    : %10s\n",ascBol[param.pbc]);
  }

  if( task == REP )
    printf(" Copies in X, Y and Z : %10d %10d %10d\n",
           param.rep[0],param.rep[1],param.rep[2]);

  printBar(stdout);
}

void printInfoCube(dataCube c){
  double min,max;
  fieldMinMax(c, &min,&max);

  printBanner(" Cube Info ",stdout);
  printf(" Number of atoms      : % 10d\n",c.natm);
  printf(" Total field points   : % 10d\n",c.npt);
  printf(" Points in X Y Z axes : % 10d % 10d % 10d\n",c.pts[0],c.pts[1],c.pts[2]);
  printf(" Coordinates of r0    : % 10.6lf % 10.6lf % 10.6lf\n\n",c.min[0],c.min[1],c.min[2]);
  printf("                      | % 10.6lf % 10.6lf % 10.6lf|\n"  ,c.mvec[0],c.mvec[1],c.mvec[2]);
  printf("        H Matrix  =   | % 10.6lf % 10.6lf % 10.6lf|\n"  ,c.mvec[3],c.mvec[4],c.mvec[5]);
  printf("                      | % 10.6lf % 10.6lf % 10.6lf|\n\n",c.mvec[6],c.mvec[7],c.mvec[8]);

  printf("        H vector  =  (  % 10.6lf % 10.6lf % 10.6lf )\n\n",c.hvec[0],c.hvec[1],c.hvec[2]);

  
  //SkewParameters(c.mvec);

  printf(" Minimum and maximun  : % 10.6lf  % 10.6lf\n",min,max);
  printf(" Field range          : % 10.6lf\n",max-min);

  printSizeCube(c);

  printBar(stdout);

}

void printCubeRot(int rotate,const char *name, dataCube cube ){

  char nameOut[128]; 
  FILE *out;

  if( rotate == YES ){

    sprintf(nameOut,"%sRotate.cube",name);
    openFile(&out,nameOut,"w+");
    printCube("Cube in the format of fraccional coordinates",cube,out);
    fclose(out);

    printf("  File %s was generated\n",nameOut);
    printBar(stdout);
  }

}


void printCube(char *text, dataCube c,FILE *out){
  int i,j;

  fprintf(out," File generated by Cube3D code\n");
  fprintf(out," %s\n",text);

  fprintf(out," %5d % 10.6lf % 10.6lf % 10.6lf\n",c.natm  ,c.min[0] ,c.min[1] ,c.min[2] );
  fprintf(out," %5d % 10.6lf % 10.6lf % 10.6lf\n",c.pts[0],c.mvec[0],c.mvec[1],c.mvec[2]);
  fprintf(out," %5d % 10.6lf % 10.6lf % 10.6lf\n",c.pts[1],c.mvec[3],c.mvec[4],c.mvec[5]);
  fprintf(out," %5d % 10.6lf % 10.6lf % 10.6lf\n",c.pts[2],c.mvec[6],c.mvec[7],c.mvec[8]);

  for(i=0;i<c.natm;i++)
    fprintf(out," %5d % 10.6lf % 10.6lf % 10.6lf % 10.6lf\n",c.zatm[i],(double) c.zatm[i],
           c.coor[3*i],c.coor[3*i+1],c.coor[3*i+2]);

   j = 0;
   for(i=0;i<c.npt;i++){
     j++;
     fprintf(out," % 10.6E ",c.field[i]);
     if( j == 6 || i == c.npt - 1){
       fprintf(out,"\n");
       j = 0;
     }
   }


}

int sizeUnits(int memory, char *out){
 
  double tmp = (double) memory;
  if( tmp <=  1024.) {
    sprintf(out," %6.2lf (B)",tmp);
    return 0;
  }
  
  tmp /= 1024.;
  if( tmp <= 1024.){
    sprintf(out," %6.2lf (kiB)",tmp);
    return 0;
  }

  tmp /= 1024.;
  if( tmp <= 1024.){
    sprintf(out," %6.2lf (MiB)",tmp);
    return 0;
  }

  tmp /= 1024.;
  sprintf(out," %6.2lf (GiB)",tmp);
  return 0;
}


void printTapas(dataCube cube ){
  
  int i,j,k;
  int npx,npy,npz;
  double x0,y0,z0;
  double xf,yf,zf;
  double hx,hy,hz;
  double x,y,z;
  double t1,t2;
  int n1,n2;
  FILE *out;
  double err1,err2,err3;

  err1 = err2 = err3 = (double) 0.;

  npx = cube.pts[0]; hx = cube.hvec[0];
  npy = cube.pts[1]; hy = cube.hvec[1];
  npz = cube.pts[2]; hz = cube.hvec[2];

  n1 = npy*npz;
  n2 = npz;

  x0 = cube.min[0];
  xf = x0 + (npx-1)*hx;

  y0 = cube.min[1]; 
  yf = y0 + (npy-1)*hy;

  z0 = cube.min[2];
  zf = z0 + (npz-1)*hz;

  printBanner("TAPA YZ",stdout);
  openFile(&out,"TapaYZ.dat","w+");
  fprintf(out,"#Tapa YZ a X constante\n");
  fprintf(out,"#   Xi = % 10.6lf\n",x0);
  fprintf(out,"#   Xf = % 10.6lf\n",xf);
  fprintf(out,"#    Y     Z     Xi    Xf    Diff\n");

  for(j=0;j<npy;j++){
    for(k=0;k<npz;k++){
      y = y0 + j*hy;
      z = z0 + k*hz;
      i = 0;
      t1 = cube.field[i*n1+j*n2+k];
      i = npx-1;
      t2 = cube.field[i*n1+j*n2+k];
      fprintf(out," % 10.6lf % 10.6lf % 10.6lf % 10.6lf % 10.6E\n",y,z,t1,t2,t2-t1);
      err1 += fabs(t2-t1);
    }
    fprintf(out,"\n");
  }
  fclose(out);

  
  printBanner("TAPA XZ",stdout);
  openFile(&out,"TapaXZ.dat","w+");
  fprintf(out,"#Tapa XZ a Y constante\n");
  fprintf(out,"#   Yi = % 10.6lf\n",y0);
  fprintf(out,"#   Yf = % 10.6lf\n",yf);
  fprintf(out,"#    X     Z     Yi    Yf    Diff\n");
  for(i=0;i<npx;i++){
    for(k=0;k<npz;k++){
      x = x0 + i*hx;
      z = z0 + k*hz;
      j=0;
      t1 = cube.field[i*n1+j*n2+k];
      j= npy-1;
      t2 = cube.field[i*n1+j*n2+k];
      fprintf(out," % 10.6lf % 10.6lf % 10.6lf % 10.6lf % 10.6E\n",x,y,t1,t2,t2-t1);
      err2 += fabs(t2-t1);
    }
    fprintf(out,"\n");
  }
  fclose(out);


  printBanner("TAPA XY",stdout);
  openFile(&out,"TapaXY.dat","w+");
  fprintf(out,"#Tapa XY a Z constante\n");
  fprintf(out,"#   Zi = % 10.6lf\n",z0);
  fprintf(out,"#   Zf = % 10.6lf\n",zf);
  fprintf(out,"#    X     Y     Zi    Zf    Diff\n");
  for(i=0;i<npx;i++){
    for(j=0;j<npy;j++){
      x = x0 + i*hx;
      y = y0 + j*hy;
      k=0;
      t1 = cube.field[i*n1+j*n2+k];
      k=npz-1;
      t2 = cube.field[i*n1+j*n2+k];
      fprintf(out," % 10.6lf % 10.6lf % 10.6lf % 10.6lf % 10.6E\n",x,z,t1,t2,t2-t1);
      err3 += fabs(t2-t1);
    }
    fprintf(out,"\n");
  }
  fclose(out);


  printf(" Error en tapas\n");
  printf(" Tapa XY : % 10.6lf\n",err1);
  printf(" Tapa XZ : % 10.6lf\n",err2);
  printf(" Tapa YZ : % 10.6lf\n",err3);

}

void printDataNCI ( int npt, dataRun param, double *fred, double
                    *frho, char *name){

  int i;
  char nameOut[120];
  double tmp;
  FILE *out;


  strcpy(nameOut,name);
  strcat(nameOut,".dat");

  openFile(&out,nameOut,"w+");

  fprintf(out,"########################################################\n");
  fprintf(out,"#   File 100*lambda_2*rho(r) vs s(r) create by Cube3D   \n");
  fprintf(out,"# Cutoff in RGD : %10.8lf\n",param.rgd);
  fprintf(out,"# Cutoff in Rho : %10.8lf\n",param.la2);
  fprintf(out,"########################################################\n");
  fprintf(out,"#        Rho (r)          RGD (r) \n");
  fprintf(out,"########################################################\n");
  for(i=0;i<npt;i++){
    tmp = fabs(frho[i])/100.;
    if( tmp < param.la2 && fred[i] < param.rgd )
      fprintf(out,"  % 15.8lf  % 15.8lf\n",frho[i],fred[i]);
  }
  fprintf(out,"########################################################\n");

  fclose(out);
}
void SkewParameters(double *mvec){
  int i,j;
  double norm;
  double sqrtG;
  double cov[3][3];
  double ctr[3][3];

  double tmp[3];

  printf(" Skew Coordinates !! \n");
  printBar(stdout);

  cov[0][0] = mvec[0];  cov[0][1] = mvec[1];  cov[0][2] = mvec[2];
  cov[1][0] = mvec[3];  cov[1][1] = mvec[4];  cov[1][2] = mvec[5];
  cov[2][0] = mvec[6];  cov[2][1] = mvec[7];  cov[2][2] = mvec[8];

  norm = getNormVec(cov[0]);
  for(i=0;i<3;i++)
    cov[0][i] /= norm;

  norm = getNormVec(cov[1]);
  for(i=0;i<3;i++)
    cov[1][i] /= norm;

  norm = getNormVec(cov[2]);
  for(i=0;i<3;i++)
    cov[2][i] /= norm;

  crossProduct(cov[1],cov[2],tmp);
  sqrtG = dotProduct(cov[0],tmp);


  printf("                         % 10.6lfi % 10.6lfj % 10.6lfk\n"  ,cov[0][0],cov[0][1],cov[0][2]);
  printf("     Covariant Basis =   % 10.6lfi % 10.6lfj % 10.6lfk\n"  ,cov[1][0],cov[1][1],cov[1][2]);
  printf("                         % 10.6lfi % 10.6lfj % 10.6lfk\n\n",cov[2][0],cov[2][1],cov[2][2]);

  printf("           g11       =   % 10.6lf \n",dotProduct(cov[0],cov[0]));
  printf("           g22       =   % 10.6lf \n",dotProduct(cov[1],cov[1]));
  printf("           g33       =   % 10.6lf \n",dotProduct(cov[2],cov[2]));

  printf("           g12       =   % 10.6lf \n",dotProduct(cov[0],cov[1]));
  printf("           g13       =   % 10.6lf \n",dotProduct(cov[0],cov[2]));
  printf("           g23       =   % 10.6lf \n",dotProduct(cov[1],cov[2]));
  printf("        sqrt(g)      =   % 10.6lf \n",sqrtG);

  crossProduct(cov[1],cov[2],ctr[0]);
  crossProduct(cov[2],cov[0],ctr[1]);
  crossProduct(cov[0],cov[1],ctr[2]);

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      ctr[i][j] /= sqrtG;

  printf("                         % 10.6lfi % 10.6lfj % 10.6lfk\n"  ,ctr[0][0],ctr[0][1],ctr[0][2]);
  printf(" Contravariant Basis =   % 10.6lfi % 10.6lfj % 10.6lfk\n"  ,ctr[1][0],ctr[1][1],ctr[1][2]);
  printf("                         % 10.6lfi % 10.6lfj % 10.6lfk\n\n",ctr[2][0],ctr[2][1],ctr[2][2]);

  printBar(stdout);
}


void checkBoundaryCond(dataCube cube, int *check){
  
  int i,j,k;
  int npx,npy,npz;
  double t1,t2;
  int n1,n2;
  double err1,err2,err3;

  err1 = err2 = err3 = (double) 0.;

  npx = cube.pts[0]; 
  npy = cube.pts[1]; 
  npz = cube.pts[2]; 

  n1 = npy*npz;
  n2 = npz;

  for(j=0;j<npy;j++){
    for(k=0;k<npz;k++){
      i = 0;
      t1 = cube.field[i*n1+j*n2+k];
      i = npx-1;
      t2 = cube.field[i*n1+j*n2+k];
      err1 += fabs(t2-t1);
    }
  }
  
  for(i=0;i<npx;i++){
    for(k=0;k<npz;k++){
      j=0;
      t1 = cube.field[i*n1+j*n2+k];
      j= npy-1;
      t2 = cube.field[i*n1+j*n2+k];
      err2 += fabs(t2-t1);
    }
  }


  for(i=0;i<npx;i++){
    for(j=0;j<npy;j++){
      k=0;
      t1 = cube.field[i*n1+j*n2+k];
      k=npz-1;
      t2 = cube.field[i*n1+j*n2+k];
      err3 += fabs(t2-t1);
    }
  }


  if( err1 > 1.E-7 || err2 > 1.E-7 || err3 > 1.E-7){
    // the point n-1 != 0
    (*check) = 0;
  }else{
    // the point n-1 == 0
    (*check) = 1;
  }

}

void checkCommandLine(int argc, char *argv[]){

    int i;
    char name[64];
    
 //   if( argc == 1 )
 //       exit(EXIT_FAILURE);

    for( i=0; i<argc; i++){
        if( argv[i][0] == '-' && argv[i][1] == 'i'){
            printf("%s A input file [input.in] will be created!%s\n",TGBI,TRST);
            if ( argv[i+1] != NULL )
                strcpy(name,argv[i+1]);
            else{
            printf("%s Change the name 'foobar.cube' in the file [input.in]%s\n",TRBI,TRST);
                strcpy(name,"foobar.cube");
            }
            createInput (name);
            exit(EXIT_SUCCESS);
        }            

    }
    

}

void getNameOut(char *nameInp, char *nameOut){
    char c;
    int i,nchar;
    nchar = strlen(nameInp);

    strcpy(nameOut,nameInp);
    
    for(i=0;i<nchar;i++){
        c = nameInp[i];
        if( c == '.')
          nameOut[i] = '\0';
    }
}

void createInput (char* nameInp ){

    char nameOut[64];

    FILE *fin;
    openFile(&fin,"input.in","w+");
    
    getNameOut(nameInp,nameOut);

    fprintf(fin,"#Name of input file\n");
    fprintf(fin,">> INPUT\n");
    fprintf(fin,"%s\n",nameInp);
    fprintf(fin,"#Name of output file\n");
    fprintf(fin,">> OUTPUT\n");
    fprintf(fin,"%s\n",nameOut);
    fprintf(fin,"# Task to evaluate\n");
    fprintf(fin,"# REDU  Reduced gradient\n");
    fprintf(fin,"# GRAD  Gradient module\n");
    fprintf(fin,"# LAP   Laplacian\n");
    fprintf(fin,"# KIN   Kinetic energy (if the cube is density)\n");
    fprintf(fin,"#       The expression is approximate. It is taken from\n");
    fprintf(fin,"#       Yu.A. Abramov, Acta Cryst. A53 1997 264–272.\n");
    fprintf(fin,"#\n");
    fprintf(fin,"# VIR   virial... (lap/4 - 2 kin)\n");
    fprintf(fin,"# NCI   Index of non-covalent interactions\n");
    fprintf(fin,"# CRIT  to search critical points \n");
    fprintf(fin,"# VOID  to determine the percentage of the void\n");
    fprintf(fin,"# REP   Replicates nx ny and nz the original cube\n");
    fprintf(fin,">> TASK\n");
    fprintf(fin,"CRIT\n");
    fprintf(fin,"# It is the polynomial used in the interpolation\n");
    fprintf(fin,"# 2 and 3 are suitable for critical points, 4 takes a long time\n");
    fprintf(fin,">> POLY\n");
    fprintf(fin,"3\n");
    fprintf(fin,"# If the system is periodic mark yes\n");
    fprintf(fin,">> PERIODIC\n");
    fprintf(fin,"yes\n");
    fprintf(fin,"# Extra properties for  NCI. The default values are\n");
    fprintf(fin,"# den 0.05 and grad  2.0\n");
    fprintf(fin,"#>> NCI PROP\n");
    fprintf(fin,"#0.1 1.5\n");
    fprintf(fin,"# Void properties the default value is 0.0003 u.a.\n");
    fprintf(fin,"# it is only loaded if VOID was written in task\n");
    fprintf(fin,"#>> VOI PROP\n");
    fprintf(fin,"#0.03\n");
    fprintf(fin,"# Properties to replicate the cube\n");
    fprintf(fin,"# only replicates if REP was TASK\n");
    fprintf(fin,"# needs 3 integer values nx ny and nz\n");
    fprintf(fin,"#>> REP PROP\n");
    fprintf(fin,"#1 2 3\n");
    fprintf(fin,"# If this line is added it is optional the field is calculated\n");
    fprintf(fin,"# given by task (grad, red,  lap, kin, vir) on one line\n");
    fprintf(fin,"# given by at1 and at2 (integers)\n");
    fprintf(fin,"# or by the plane formed by at1, at2 and at3 (integers)\n");
    fprintf(fin,"#>> GEOM\n");
    fprintf(fin,"#LINE at1 at2\n");
    fprintf(fin,"# or\n");
    fprintf(fin,"#>> GEOM\n");
    fprintf(fin,"#PLANE at1 at2 at3\n");

    fclose(fin);

}
