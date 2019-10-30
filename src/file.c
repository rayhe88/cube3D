/**
 * @file  file.c
 * @brief 
 * @author Raymundo Hernández-Esparza.
 * @date   August 2018.
 */

#include "file.h"

/**
 * @brief  Esta función abre y verifica que se abrio
 *         de forma correcta el archivo.
 * @param  **file doble puntero al archivo.
 * @param  *name parámetro const char con el nombre del archivo.
 * @param  *type parámetro const char con el tipo de archivo.
 */
int openFile( FILE **file,const char *name, const char *type ){

  (*file) = fopen(name,type);

  if( (*file) == NULL ){
    printf(" Fail to open [%s]\n",name);
    exit(EXIT_FAILURE);
  }

  return 0;
}

/**
 * @brief  Esta función abre un archivo temporal, cuyo nombre depende
 *         del momento de la ejecución.
 * @param  **file double puntero al archivo
 * @param  *prefix parámetro que guarda el prefijo del archivo a abrir.
 * @param  *name   es el nombre que regresará la función.
 * @param  *type   es el tipo de archivo que se abrirá.
 */
int tmpFile ( FILE **file, char *prefix, char *name, const char *type ){
  char path[] = "/tmp/";
  char timeascii[120];
  time_t t;

  t = time(NULL);
  sprintf(timeascii,"-%d-%ld",getpid(),t);
  strcpy(name,path);
  strcat(name,prefix);
  strcat(name,timeascii);
  strcat(name,".tmp");

  (*file) = fopen(name,type);

  
  if( (*file) == NULL ){
    printf(" Fail to open [%s]\n",name);
    exit(EXIT_FAILURE);
  }
#ifdef DEBUG
  printf(" Temporary file's name [%s]\n",name);
#endif

  return 0;
}
