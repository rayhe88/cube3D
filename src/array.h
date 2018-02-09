/**
  @file array.h
  @brief Definition of functions that allocate dynamic memory

  This set of functions allocate dynamic memory for data type
  integer, float and double. For the use en all the code.

 */
#include <stdio.h>
#include <stdlib.h>

#ifndef _ARRAY_H_
 #define _ARRAY_H_
 
 
 int createArrayInt (int, int    **, char const*);
 int createArrayFlo (int, float  **, char const*);
 int createArrayDou (int, double **, char const*);
 int createArrayLong(int, long int**,char const*);

#endif
