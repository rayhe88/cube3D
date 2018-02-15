#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "struct.h"



#ifndef _LECTURA_H_
 #define _LECTURA_H_

 int readInput(char*,char*,int*,int*,int*,double*,char*);
 
 int unloadData (int**,double**,double**);
 int loadData (dataCube*,const char *);

 int readData1(int*,int*,double*,double*,FILE*);
 int readData2(int*,double*,double*,int,int,FILE*);
 int createArrays(int**,double**,double**,int,int);

 int checkInput(int*);


 int checkData (int,int*,int*,double*,double*,double*,double*,const char *,FILE*);

#endif
