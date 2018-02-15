#include "file.h"
#include "array.h"
#include "struct.h"
#include "lectura.h"

int unloadData( int **atmnum,double **coor, double **field){
  free(*atmnum);
  free(*coor);
  free(*field);

  return 0;
}

int loadData(dataCube *data, const char *name){

  extern int *zatm;
  extern double *coor;
  extern double *field;

  int i,n,npt;
  int points[3];
  double min[3];
  double vec[9];
  FILE *inp;

  openFile(&inp,name,"r");

  readData1(&n,points,min,vec,inp);
  npt =(int) points[0]*points[1]*points[2];
  createArrays(&zatm,&coor,&field,n,npt);
  readData2(zatm,coor,field,n,npt,inp);

  fclose(inp);

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

  return 0;
}


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


int readData2( int *natom, double *coor,double *field,int n, int npt, FILE *inp){

  char buffer[128];

  int i,zatom;
  int j;
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
int checkData(int nnuc, int *pts, int *zatm,double *min,double *vec,double *coor, double *field,const char *name,FILE *out){
   

   int i,k;
   int j,npt = (int) pts[0]*pts[1]*pts[2];

   fprintf(out,"This file was reading: [%s]\n",name);
   fprintf(out,"This file was reading: [%s]\n",name);
   fprintf(out," %2d % 10.6lf % 10.6lf % 10.6lf\n",nnuc,min[0],min[1],min[2]);
   fprintf(out," %2d % 10.6lf % 10.6lf % 10.6lf\n",pts[0],vec[0],vec[1],vec[2]);
   fprintf(out," %2d % 10.6lf % 10.6lf % 10.6lf\n",pts[1],vec[3],vec[4],vec[5]);
   fprintf(out," %2d % 10.6lf % 10.6lf % 10.6lf\n",pts[2],vec[6],vec[7],vec[8]);
   for(i=0;i<nnuc;i++)
     fprintf(out," %2d % 10.6lf % 10.6lf % 10.6lf % 10.6lf\n",zatm[i],(double)zatm[i],
                                                         coor[3*i  ],
                                                         coor[3*i+1],
                                                         coor[3*i+2]);
   for(j=0,k=0;j<npt;j++){
     fprintf(out," % 10.6E",field[j]);
     k++;
     if( k == 6){
       fprintf(out,"\n");
       k = 0;
     }
     

   }

   return 0;
}

int readInput (char *namefld, char *nameout, 
               int *task, int *poly,int *bound, double *dg, char *nameinp){

  int i;
  int tmptask;
  int tmppoly;
  int tmpper=0;
  int bin[8];
  double tmpden;
  double tmpgra;
  FILE *inp;
  FILE *tmp;

  char buffer[120],bufper[10];
  char nametmp[120];
  char nametask[120];
  char tmp1[120],tmp2[120];

  const char flag_0[] = ">> INPUT";
  const char flag_1[] = ">> OUTPUT";
  const char flag_2[] = ">> TASK";
  const char flag_3[] = ">> POLY";
  const char flag_4[] = ">> NCI PROP";
  const char flag_5[] = ">> PERIODIC";

  for(i=0;i<8;i++)
    bin[i] = 0;
 
  openFile(&inp,nameinp,"r");

  tmpFile(&tmp,".",nametmp,"w+");

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
    if( !strncmp(buffer,flag_0,7) ){
      fscanf(tmp,"%s",&tmp1);
      strcpy(namefld,tmp1);
      bin[0] = 1;
    }
    if( !strncmp(buffer,flag_1,7) ){
      fscanf(tmp,"%s",&tmp2);
      strcpy(nameout,tmp2);
      bin[1] = 1;
    }
    if( !strncmp(buffer,flag_2,7) ){
      fscanf(tmp,"%s",&nametask);
      bin[2] = 1;
    }
    if( !strncmp(buffer,flag_3,7) ){
      fscanf(tmp,"%d",&tmppoly);
      bin[3] = 1;
    }

    if( !strncmp(buffer,flag_4,7) ){
      fscanf(tmp,"%lf %lf",&tmpden,&tmpgra);
      bin[4] = 1;
    }

    if( !strncmp(buffer,flag_5,7) ){
      fscanf(tmp,"%s",&bufper);
      bin[5] = 1;
    }

  }


  for(i=0;i<strlen(bufper);i++)
    bufper[i] = toupper(bufper[i]);

  for(i=0;i<strlen(nametask);i++)
    nametask[i] = toupper(nametask[i]);

  if( !strncmp(nametask,"TEST",2) )
    tmptask = 0 ;

  if( !strncmp(nametask,"GRAD",2) )
    tmptask = 1;

  if( !strncmp(nametask,"LAP",2) )
    tmptask = 2;

  if( !strncmp(nametask,"NCI",2) )
    tmptask = 3;

  if( !strncmp(nametask,"CRIT",2) )
    tmptask = 4;

  if( !strncmp(bufper,"YES",2))
    tmpper=1;

  if( tmptask == 4 ){
    if( tmppoly < 1 || tmppoly > 10)
      tmppoly = 2;
  }

  if( bin[4] == 0 ){
    tmpden = 5.0;
    tmpgra = 2.0;
  }

  if( bin[5] == 0 ){
    tmpper = 0;
  }

  
  //printf(" name field  : %s\n",namefld);
  //printf(" name output : %s\n",nameout);
  //printf("  task name  : %s\n",nametask);
  //printf("     task    : %3d\n",tmptask);
  //printf("     poly    : %3d\n",tmppoly);

  (*poly) = tmppoly;

  (*task) = tmptask;

  (*bound) = tmpper;

  dg[0] = tmpden;
  dg[1] = tmpgra;

  fclose(tmp);
  fclose(inp);
  remove(nametmp);

  if(checkInput(bin) == -1)
    exit(EXIT_FAILURE);

  
  return 0;
}

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
