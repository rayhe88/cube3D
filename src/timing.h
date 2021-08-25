#include <sys/time.h>
#include <time.h>

#ifndef _TIMING_H_
#define _TIMING_H_
typedef struct timeval timeType;

double timing(timeType *tf, timeType *ti);

#endif
