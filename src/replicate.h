#include "struct.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef _REPLICATE_H_
#define _REPLICATE_H_

int replicate(dataCube, int[], const char *);

int replicateCoor(int, dataCube, int *, double *, int *);

int getIndexOld(int, int, int[], int[]);

#endif
