#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "array.h"

void merge( int l, int m, int r, double *array,int *iarray){
  
  int i,j,k;
  int n1 = m - l + 1;
  int n2 = r - m;

  //double larray[n1],rarray[n2];
  //int    lintarr[n1],rintarr[n2];

  double *larray,*rarray;
  int    *lintarr,*rintarr;

  createArrayInt(n1,&lintarr,"integer left");
  createArrayInt(n2,&rintarr,"integer right");

  createArrayDou(n1,&larray,"double left");
  createArrayDou(n2,&rarray,"double right");

  for( i = 0; i < n1; i++){
    larray[i]  = array[ l + i ];
    lintarr[i] = iarray[ l + i ];
  }

  for( i = 0; i < n2; i++){
    rarray[i]  = array[ m + 1 + i];
    rintarr[i] = iarray[ m + 1 + i];
  }

  i = 0;
  j = 0;
  k = l;

  while( i < n1 && j < n2 ){
    if( larray[i] >=  rarray[j] ){
      array[k] = larray[i];
      iarray[k] = lintarr[i];
      i++;
    }else{
      array[k] = rarray[j];
      iarray[k] = rintarr[j];
      j++;
    }
    k++;
  }
  
  while( i < n1 ){
    array[k] = larray[i];
    iarray[k] = lintarr[i];
    i++;
    k++;
  }

  while( j < n2 ){
    array[k] = rarray[j];
    iarray[k] = rintarr[j];
    j++;
    k++;

  }

  free(lintarr);
  free(rintarr);

  free(larray);
  free(rarray);
}

void mergeSort( int l, int r, double *array,int *iarray){
  if( l < r ){
    int m = l + (r-l)/2;

    mergeSort(l,m,array,iarray);
    mergeSort(m+1,r,array,iarray);

    merge(l,m,r,array,iarray);

  }
}

void printArray(int n, double *array, int *iarray){

  int i;
  for(i=0;i<n;i++)
    printf(" % 12d   % 12.8lf\n",iarray[i],array[i]);
  
}

double genRandom(double min, double max){
  return drand48()*(max-min) + min;
}

void fillArrays(int n, int *iarray, double *darray){
  int i;

  for(i=0;i<n;i++){
    iarray[i] = i;
    darray[i] = genRandom(0,100);
  }

}
/*
int main (int argc, char *argv[]){

  int np;
  int *iarray;
  double *darray;
  
  if( argc > 1 )
    np = atoi(argv[1]);
  else
    np = 10;

  iarray = malloc( np*sizeof(int) );
  darray = malloc( np*sizeof(double) );

  if( iarray == NULL || darray == NULL ){
    printf(" There are some errors in memory\n");
    exit(EXIT_FAILURE);
  }

  fillArrays(np,iarray,darray);

  printf(" Array original\n");
  printArray(np, darray,iarray);

  mergeSort(0,np-1,darray,iarray);

  printf(" Array Sorted\n");
  printArray(np, darray,iarray);

  return 0;
}
*/
