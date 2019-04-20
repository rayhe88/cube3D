/**
 * @file   file.h
 * @brief 
 * @author Raymundo Hernández-Esparza.
 * @date   August 2018.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>



#ifndef _FILE_H_
 #define _FILE_H_


 int openFile(FILE **, const char*,const char*);

 int tmpFile(FILE **,char*,char*,const char*);

#endif
