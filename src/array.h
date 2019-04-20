/**
 *  @file   array.h
 *  @brief  Defining functions to allocate dynamic memory.
 *
 *  This set of functions allocates dynamic memory for 
 *  integer, float, double and long types. For use in the
 *  code.
 *
 *  @author Raymundo Hernández-Esparza.
 *  @date   August 2017.
 */
#include <stdio.h>
#include <stdlib.h>

#ifndef _ARRAY_H_
 #define _ARRAY_H_
 
 
 int createArrayInt (int, int    **, const char*);
 int createArrayFlo (int, float  **, const char*);
 int createArrayDou (int, double **, const char*);
 int createArrayLong(int, long int**,const char*);

#endif
