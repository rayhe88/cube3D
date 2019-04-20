/**
 * @file   lectura.h
 * @brief 
 * @author Raymundo Hernández-Esparza.
 * @date   August 2018.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "struct.h"

#ifndef _LECTURA_H_
 #define _LECTURA_H_

 int readInput(char*,char*,dataRun*,char*);
 
 int unloadData (dataCube*,int**,double**,double**);
 int loadData (dataCube*,int**,double**,double**,int*,const char *);

 int readData1(int*,int*,double*,double*,FILE*);
 int readData2(int*,double*,double*,int,int,FILE*);
 int createArrays(int**,double**,double**,int,int);

 int checkInput(int*);

 int getTask(char*);
 int checkData (int,int*,int*,double*,double*,double*,double*,const char *,FILE*);

#endif
