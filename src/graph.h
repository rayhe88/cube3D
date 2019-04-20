#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include  "struct.h"

#ifndef _GRAPH_H_
 #define _GRAPH_H_

 #define GNUPLOT_PATH "/usr/bin/gnuplot-qt -persist"
 #define GNUPLOT_PATH_WX "/usr/bin/gnuplot-wx -persist &> /dev/null"

 #define NPX 621
 #define NPY 384

 #define _2D 102 
 #define _3D 103

 #define DEFAULT "#CFD8DC"
 void InitPlot  (FILE*,char*,char*,int,int,int,int);
 void InitPlotWx(FILE*,char*,char*,int,int,int,int);
 
 void plotNCIdat(char*,dataRun );

 void getNCIgnu (char*,dataRun );


#endif
