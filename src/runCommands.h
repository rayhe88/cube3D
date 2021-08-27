/**
 *  @file   loadRunCommands.h
 *  @brief  header file for loading run commands.
 *
 *  @author Raymundo Hernández-Esparza.
 *  @date   August 2021.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <unistd.h>
#include <pwd.h>
#include <sys/types.h>

#include "file.h"
#include "struct.h"

#ifndef _LOAD_RUN_COMMANDS_H_
#define _LOAD_RUN_COMMANDS_H_

int checkFile(dataRC *config);
int createFileRC(char* namerc, dataRC config);
int readFileRC(char *path, dataRC *config);
int checkRunCommands(dataRC *config);
int defaultRunCommands(dataRC *config);
#endif
