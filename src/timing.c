#include "timing.h"

double timing ( timeType *tf, timeType *ti ){
  double value = (double) tf -> tv_sec 
               + (double) tf -> tv_usec/1000000.
               - (double) ti -> tv_sec
               - (double) ti -> tv_usec/1000000.;
  return 1000.*value;
}
