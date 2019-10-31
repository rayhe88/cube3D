/**
 * @file lectura.c
 * @brief
 * @author Raymundo Hernández-Esparza.
 * @date   August 2018.
 */

#include "file.h"
#include "array.h"
#include "struct.h"
#include "lectura.h"
#include "rotation.h"
#include "cubeIndex.h"

int unloadData( dataCube *cube, int **atmnum,double **coor, double **field){
  cube->zatm = NULL;
  cube->coor = NULL;
  cube->field = NULL;
  free(*atmnum);
  free(*coor);
  free(*field);

  return 0;
}

/**
 * @brief Esta función carga la información del archivo cube
 *        en una estructura llamada dataCube.
 *
 *        Esta función llama a readData1 y readData2 para
 *        almacenar los datos y rear los arrays necesarios
 *        para la estructura dataCube.
 *
 * @param *data Información numérica del cubo.
 * @param **zatm Puntero al puntero de números atómicos.
 * @param **coor Puntero al puntero de coordenadas.
 * @param **field Puntero al puntero del campo.
 * @param *name Nombre del archivo cub/cube
 */
int loadData(dataCube *data, int **zatm, double **coor, double **field, int *rot,const char *name){

  int i,n,npt;
  int points[3];
  double min[3];
  double vec[9];
  FILE *inp;

  openFile(&inp,name,"r");

  // Primer paso en la lectura. En esta parte cargamos
  // lo que es invariante del sistema
  readData1(&n,points,min,vec,inp);

  npt =(int) points[0]*points[1]*points[2];

  createArrays(zatm,coor,field,n,npt);

  // Segundo paso en la lectura. En esta parte cargamos
  // lo que depende del sistema (números atómicos,
  // coordenadas y el valor del campo)
  readData2((*zatm),(*coor),(*field),n,npt,inp);

  fclose(inp);
  
  // Asignamos la lectura que se realizo de manera temporal
  // a la estructura dataCube.
  data->natm = n;
  for(i=0;i<3;i++){
    data->min[i] = min[i];
    data->pts[i] = points[i];
  }
  for(i=0;i<9;i++){
    data->mvec[i] = vec[i];
  }

  data->npt = npt;

  data->hvec[0] = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
  data->hvec[1] = sqrt(vec[3]*vec[3] + vec[4]*vec[4] + vec[5]*vec[5]);
  data->hvec[2] = sqrt(vec[6]*vec[6] + vec[7]*vec[7] + vec[8]*vec[8]);

  data->max[0] = min[0] + points[0]*(data->hvec[0]) ;
  data->max[1] = min[1] + points[1]*(data->hvec[1]) ;
  data->max[2] = min[2] + points[2]*(data->hvec[2]) ;


  data->zatm = (*zatm);
  data->coor = (*coor);
  data->field = (*field);

  if( checkRotation(data->mvec) == YES ){
    rotationCube(data->mvec,data);
    (*rot) = YES;
  }

  return 0;
}

/** 
 * @brief 
 * @param
 * @param
 * @param
 * @param
 * @param
 */
int readData1(int *nnuc, int *points, double *min, double *vec, FILE *inp){

  char buffer[128];
  int n,p0,p1,p2;
  double x0,y0,z0;
  double v1,v2,v3,v4,v5,v6,v7,v8,v9;

  fgets(buffer,128,inp); // Comentario 1
  fgets(buffer,128,inp); // Comentario 2

  fscanf(inp,"%d %lf %lf %lf",&n,&x0,&y0,&z0);
  fscanf(inp,"%d %lf %lf %lf",&p0,&v1,&v2,&v3);
  fscanf(inp,"%d %lf %lf %lf",&p1,&v4,&v5,&v6);
  fscanf(inp,"%d %lf %lf %lf",&p2,&v7,&v8,&v9);

  (*nnuc) = n;

  points[0] = p0;
  points[1] = p1;
  points[2] = p2;

  min[0] = x0;  min[1] = y0;  min[2] = z0;

  vec[0] = v1;  vec[1] = v2;  vec[2] = v3;
  vec[3] = v4;  vec[4] = v5;  vec[5] = v6;
  vec[6] = v7;  vec[7] = v8;  vec[8] = v9;

  
  rewind(inp);

  return 0;
}

/** 
 * @brief 
 * @param
 * @param
 * @param
 * @param
 * @param
 * @param
 */
int readData2( int *natom, double *coor,double *field,int n, int npt, FILE *inp){

  char buffer[128];

  int i,j;
  int zatom;
  double tmp;
  double x,y,z;

  fgets(buffer,128,inp);
  fgets(buffer,128,inp);
  fgets(buffer,128,inp);
  fgets(buffer,128,inp);
  fgets(buffer,128,inp);
  fgets(buffer,128,inp);

  for( i = 0; i < n ; i++ ){
    fscanf(inp,"%d %lf %lf %lf %lf",&zatom,&tmp,&x,&y,&z);
    natom[i] = zatom;
    coor[3*i  ] = x;
    coor[3*i+1] = y;
    coor[3*i+2] = z;
  }

  for( j = 0; j < npt ; j++ ){
    fscanf(inp,"%lf",&tmp);
    field[j] = tmp;
  }

  return 0;
}

/** 
 * @brief 
 * @param
 * @param
 * @param
 * @param
 * @param
 */
int createArrays(int **natom, double **coor, double **field, int n, int npt){
  if( n < 1 || npt <= 0 ){
    printf("There is an error in the arrays' size\n");
    exit(EXIT_FAILURE);
  }


  *natom = (int*) malloc( n*sizeof(int));
  if( *natom == NULL ){
    printf(" Failed to allocate memory for atomic number\n");
    exit(EXIT_FAILURE);
  }

  *coor = (double*) malloc(3*n*sizeof(double));
  if( *coor == NULL ){
    printf(" Failed to allocate memory for coordinates\n");
    exit(EXIT_FAILURE);
  }
  
  *field = (double*) malloc(npt*sizeof(double));
  if( *field == NULL ){
    printf(" Failed to allocate memory for field\n");
    exit(EXIT_FAILURE);
  }

  return 0;
}

/**
 * @brief  This function reads the input file to load 
 *         the general data.
 * @param  namefld is the name of the cube file.
 * @param  nameout is the gneral name for output files.
 * @param  task is the value correspondent to job.
 * @param  poly is the degree of the polynomial.
 * @param  bound is a variable to activate the periodic conditios.
 * @param  dg is an array to store the prop of NCI.
 * @param  nameinp  is the name of the input file.
 * @return The correspondent value for task, properties and name of input file.
 */
int readInput (char *namefld, char *nameout, dataRun *param, char *nameinp){

  int i;
  int tmptask;
  int tmppoly;
  int tmpper=0;
  int rpx,rpy,rpz;
  int bin[10];
  int izq,der,size;
  int at1,at2,at3;
  int gt=0;
  double tmpden,tmpgra;
  double tmpvac;
  FILE *inp;
  FILE *tmp;

  char buffer[120],bufper[10],bufgeom[10];
  char nametmp[120];
  char buffertmp[120];
  char nametask[120];
  char tmp1[120],tmp2[120];

  const char flag_0[] = ">> INPUT";
  const char flag_1[] = ">> OUTPUT";
  const char flag_2[] = ">> TASK";
  const char flag_3[] = ">> POLY";
  const char flag_4[] = ">> NCI PROP";
  const char flag_5[] = ">> PERIODIC";
  const char flag_6[] = ">> REP PROP";
  const char flag_7[] = ">> VOI PROP";
  const char flag_8[] = ">> GEOM";

  for(i=0;i<10;i++)
    bin[i] = 0;
 
  openFile(&inp,nameinp,"r");

  tmpFile(&tmp,".c3dLec",nametmp,"w+");

  while( !feof(inp) ){
    fgets(buffer,100,inp);

    if( buffer[0] != '#' ){
      i = 0;
      while( i < 100 && !feof(inp) && buffer[i] != '\n'){
        fprintf(tmp,"%c",buffer[i]);
        i++;
      }
      fprintf(tmp,"\n");
    }

  }
 
  rewind(inp);
  rewind(tmp);

  while( !feof(tmp)){
    fgets(buffer,100,tmp);
    if( !strncmp(buffer,flag_0,7) ){ // INPUT
      fscanf(tmp,"%s",tmp1);
      strcpy(namefld,tmp1);
      bin[0] = 1;
    }
    if( !strncmp(buffer,flag_1,7) ){ // OUTPUT
      fscanf(tmp,"%s",tmp2);
      strcpy(nameout,tmp2);
      bin[1] = 1;
    }
    if( !strncmp(buffer,flag_2,7) ){ // TASK
      fscanf(tmp,"%s",nametask);
      bin[2] = 1;
    }
    if( !strncmp(buffer,flag_3,7) ){ // POLY
      fscanf(tmp,"%d",&tmppoly);
      bin[3] = 1;
    }
    if( !strncmp(buffer,flag_4,7) ){ // NCI PROP
      fscanf(tmp,"%lf %lf",&tmpden,&tmpgra);
      bin[4] = 1;
    }
    if( !strncmp(buffer,flag_5,7) ){ // PERIODIC
      fscanf(tmp,"%s",bufper);
      bin[5] = 1;
    }
    if( !strncmp(buffer,flag_6,7) ){ // REP PROP
      fscanf(tmp,"%d %d %d",&rpx,&rpy,&rpz);
      if( rpx <= 0 ){
        rpx = 1;
        printf("  [ERROR] n_rep_x must be greater than zero\n");
      }
      if( rpy <= 0 ){ 
        rpy = 1;
        printf("  [ERROR] n_rep_y must be greater than zero\n");
      }
      if( rpz <= 0 ){ 
        rpz = 1;
        printf("  [ERROR] n_rep_z must be greater than zero\n");
      }
      bin[6] = 1;
    }
    if( !strncmp(buffer,flag_7,7) ){ // VOI PROP
      fscanf(tmp,"%lf",&tmpvac);
      bin[7] = 1;
    }
    if( !strncmp(buffer,flag_8,7) ){ // GEOM
      at1 = at2 = at3 = -1;

      fgets(buffertmp,120,tmp);
      sscanf(buffertmp,"%s",bufgeom);

      for(i=0;i<strlen(bufgeom);i++)
        bufgeom[i] = toupper(bufgeom[i]);

      if( !strncmp(bufgeom,"LINE",4) ){
        gt = LIN;
        sscanf(buffertmp,"%s %d %d",bufgeom,&at1,&at2);
      }
      if( !strncmp(bufgeom,"PLANE",4) ){
        gt = PLA;
        sscanf(buffertmp,"%s %d %d %d",bufgeom,&at1,&at2,&at3);
      }
      bin[8] = 1;
    }
  }


  for(i=0;i<strlen(bufper);i++)
    bufper[i] = toupper(bufper[i]);

  for(i=0;i<strlen(nametask);i++)
    nametask[i] = toupper(nametask[i]);


  tmptask = getTask(nametask);


  if( !strncmp(bufper,"YES",2))
    tmpper=YES;

  if( tmptask == 4 ){
    if( tmppoly < 1 || tmppoly > 6)
      tmppoly = 3;
  }
 

    


// Load defaults
//
// Default for NCI prop
  if( bin[4] == 0 ){
    tmpden  = 0.05;
    tmpgra  = 2.0;
  }
// Default for periodic boundary
  if( bin[5] == 0 ){
    tmpper = NOT; 
  }
// Default for replicate
  if( bin[6] == 0 ){  
    rpx = rpy = rpz = 1;
  }  
// Default for void
  if( bin[7] == 0 ){
    tmpvac = 0.0003;
  }



  size = getLeftRight(tmppoly,&izq,&der);

  param -> pol     = tmppoly;
  param -> task    = tmptask;
  param -> pbc     = tmpper;
  param -> la2     = tmpden;
  param -> rgd     = tmpgra;
  param -> vac     = tmpvac;
  param -> izq     = izq;
  param -> der     = der;
  param -> size    = size;

  param -> rep[0]  = rpx;
  param -> rep[1]  = rpy;
  param -> rep[2]  = rpz;

  param -> geotask = gt;

  param -> geom[0] = at1;
  param -> geom[1] = at2;
  param -> geom[2] = at3;

  fclose(tmp);
  fclose(inp);
  remove(nametmp);

  if(checkInput(bin) == -1)
    exit(EXIT_FAILURE);

  
  return 0;
}

/**
 * @brief  This function gets the macro for the task.
 * @param  name for task read.
 * @return The correspondent value for task.
 */
int getTask(char* task){
  int valor = -1;

  if( !strncmp(task,"REDU",3) )
    valor = RED;
  if( !strncmp(task,"GRAD",3) )
    valor = GRA;
  if( !strncmp(task,"LAPL",3) )
    valor = LAP;
  if( !strncmp(task,"KIN",3) )
    valor = KIN;
  if( !strncmp(task,"VIR",3) )
    valor = VIR;
  if( !strncmp(task,"NCI",3) )
    valor = NCI;
  if( !strncmp(task,"CRIT",3) )
    valor = CRI;
  if( !strncmp(task,"VOID",3) )
    valor = VOI;
  if( !strncmp(task,"REPI",3) )
    valor = REP;
  
  if( valor == -1 ){
    printf("  [ERROR] There is an error in >> TASK\n");
    exit(EXIT_FAILURE);
  }


  return valor;
}

/** 
 * @brief 
 * @param
 * @param
 * @param
 * @param
 * @param
 */
int checkInput(int bin[8]){
  
  int i,error=0;
  int ret;

  for(i=0;i<4;i++)
    error += bin[i];
  
  ret = 0;

  if( error != 4) {
    ret = -1;
    printf(" There is an error in the\n");
    if( bin[0] == 0 )
      printf(" >> INPUT section\n");
    if( bin[1] == 0 )
      printf(" >> OUTPUT section\n");
    if( bin[2] == 0 )
      printf(" >> TASK section \n");
    if( bin[3] == 0 )
      printf(" >> POLY section \n");
   }

   return ret;
}
