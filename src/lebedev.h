#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef _LEBEDEV_H
 #define _LEBEDEV_H_
 typedef struct{
   double x;
   double y;
   double z;
   double w;
 } dataVec;

 void translateAndChange( int, dataVec*, double[],double);

 int createArrayVec(int, dataVec**,const char*);

 void genGroup1 ( int *num, dataVec *r, double va, double vb, double v);
 void genGroup2 ( int *num, dataVec *r, double va, double vb, double v);
 void genGroup3 ( int *num, dataVec *r, double va, double vb, double v);
 void genGroup4 ( int *num, dataVec *r, double va, double vb, double v);
 void genGroup5 ( int *num, dataVec *r, double va, double vb, double v);
 void genGroup6 ( int *num, dataVec *r, double va, double vb, double v);

 int Lebedev0006 ( dataVec *r );

 int Lebedev0014 ( dataVec *r );

 int Lebedev0026 ( dataVec *r );
 
 int Lebedev0038 ( dataVec *r );

 int Lebedev0050 ( dataVec *r );

 int Lebedev0074 ( dataVec *r );

 int Lebedev0086 ( dataVec *r );

 int Lebedev0110 ( dataVec *r );

 int Lebedev0146 ( dataVec *r );

 int Lebedev0170 ( dataVec *r );

 int Lebedev0194 ( dataVec *r );

 int Lebedev0230 ( dataVec *r );

 int Lebedev0266 ( dataVec *r );

 int Lebedev0302 ( dataVec *r );

 int Lebedev0350 ( dataVec *r );

 int Lebedev0434 ( dataVec *r );

#endif
