#include "array.h"


int createArrayInt ( int n, int **ptr, char const *mess){

  if( n < 1 ){
    printf(" There is an error in the size for [%s]\n",mess);
    exit(EXIT_FAILURE);
  }
  
  *ptr = (int*) malloc( n*sizeof(int));
  if( *ptr == NULL ){
    printf(" Failed to allocate memory: [%s]\n",mess);
    exit(EXIT_FAILURE);
   }

   return 0;
}

int createArrayFlo ( int n, float **ptr, char const *mess){

  if( n < 1 ){
    printf(" There is an error in the size for [%s]\n",mess);
    exit(EXIT_FAILURE);
  }
  
  *ptr = (float*) malloc( n*sizeof(float));
  if( *ptr == NULL ){
    printf(" Failed to allocate memory: [%s]\n",mess);
    exit(EXIT_FAILURE);
   }

   return 0;
}
  
int createArrayDou ( int n, double **ptr, char const *mess){

  if( n < 1 ){
    printf(" There is an error in the size for [%s]\n",mess);
    exit(EXIT_FAILURE);
  }
  
  *ptr = (double*) malloc( n*sizeof(double));
  if( *ptr == NULL ){
    printf(" Failed to allocate memory: [%s]\n",mess);
    exit(EXIT_FAILURE);
   }

   return 0;
}

int createArrayLong( int n, long int **ptr, char const *mess){
  if( n < 1 ){
    printf(" There is an error in the size for [%s]\n",mess);
    exit(EXIT_FAILURE);
  }
  
  *ptr = (long int*) malloc( n*sizeof(long int));
  if( *ptr == NULL ){
    printf(" Failed to allocate memory: [%s]\n",mess);
    exit(EXIT_FAILURE);
   }

   return 0;
}
