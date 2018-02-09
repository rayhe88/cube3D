#include "file.h"


int openFile( FILE **file,const char *name, const char *type ){

  (*file) = fopen(name,type);

  if( (*file) == NULL ){
    printf(" Fail to open [%s]\n",name);
    exit(EXIT_FAILURE);
  }

  return 0;
}


int tmpFile ( FILE **file, char *prefix, char *name, const char *type ){
  
  char timeascii[120];
  time_t t;

  t = time(NULL);
  sprintf(timeascii,"%d",t);
  strcpy(name,prefix);
  strcat(name,timeascii);
  strcat(name,".tmp");

  (*file) = fopen(name,type);

  
  if( (*file) == NULL ){
    printf(" Fail to open [%s]\n",name);
    exit(EXIT_FAILURE);
  }

  return 0;
}
