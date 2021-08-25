/**
 * @file tableP.h
 * @brief
 * @author Raymundo Hernández-Esparza.
 * @date   August 2018.
 */
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef LSTRING
#define LSTRING 10
#endif

#ifndef _TABLE_P_H_
#define _TABLE_P_H_

#define TPNMAX 119

int clearString(char *);
int printSymbol(int, FILE *);
int printProperties(int, FILE *);
int getAtomicNumber(const char *);
int getAtomicSymbol(int, int, char *);

float getChemProperty(int, const char *);

#endif
