/**
 * @file   array.c
 * @brief  Implementation of functions to allocate dynamic memory.
 *
 * This set of functios allocates dynamic memory for
 * integer, float, double and long types. For use in the
 * code.
 *
 * @author Raymundo Hernández-Esparza.
 * @date   August 2017.
 */
#include "array.h"


/**
 *  @brief Allocate memory for int.
 *  @param n size to assign.
 *  @param **ptr pointer to memory array.
 *  @param *mess  exit message for possible erors.
 */
int createArrayInt ( int n, int **ptr, const char *mess){

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

/**
 *  @brief Allocate memory for float.
 *  @param n  size to assign.
 *  @param **ptr pointer to memory array.
 *  @param *mess  exit message for possible erors.
 */
int createArrayFlo ( int n, float **ptr, const char *mess){

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
  
/**
 *  @brief Allocate memory for double.
 *  @param n size to assign.
 *  @param **ptr pointer to memory array.
 *  @param *mess  exit message for possible erors.
 */
int createArrayDou ( int n, double **ptr, const char *mess){

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

/**
 *  @brief Allocate memory for long int.
 *  @param n size to assign.
 *  @param **ptr pointer to memory array.
 *  @param *mess  exit message for possible erors.
 */
int createArrayLong( int n, long int **ptr, const char *mess){
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
