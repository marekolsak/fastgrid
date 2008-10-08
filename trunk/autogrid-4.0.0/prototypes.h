/*

 $Id: prototypes.h,v 1.8 2007/05/03 20:46:06 garrett Exp $

 AutoGrid 

 Copyright (C) 1989-2007,  Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson, 
 All Rights Reserved.

 AutoGrid is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 */

#pragma once
#ifndef _WIN32
    #include <sys/times.h>
#endif
#include <sys/types.h>
#include <cstdio>
#include "parameters.h"
#include "structs.h"
#include "autocomm.h"

// implemented in banner.cpp
void banner(double version_num, FILE *logFile);

// implemented in process_program_parameters.cpp
int process_program_parameters(int argc, char **argv, FILE *&GPF, FILE *&logFile, char *&programname, char (&grid_param_fn)[MAX_CHARS], int &debug, int &oldpdbq);

// implemented in atom_parameter_manager.cpp
ParameterEntry *apm_find(const char key[]);
void apm_enter(const char key[], ParameterEntry value);

// implemented in check_size.cpp
int check_size(int nelements, char axischar, const char *programname, FILE *logFile);

// implemented in gpfparser.cpp
int gpfparser(char line[LINE_LEN]);

// implemented in parsetypes.cpp
int parsetypes(char * line, char *words[], int maxwords);

// implemented in prHMSfixed.cpp
void prHMSfixed(float t, FILE *logFile);

// implemented in printdate.cpp
void printdate(FILE *fp, int flag);

// implemented in printhms.cpp
void printhms(float t, FILE *logFile);

// implemented in read_parameter_library.cpp
void read_parameter_library(char FN_parameter_library[MAX_CHARS], int outlev, const char *programname, int debug, FILE *logFile, Linear_FE_Model &AD4);
void setup_parameter_library(int outlev, const char *programname, int debug, FILE *logFile, Linear_FE_Model &AD4);

// implemented in strindex.cpp
int strindex(char s[], char t[]);

// implemented in timesys.cpp
void timesys(Clock duration, struct tms *start, struct tms *end, Real idct, FILE *logFile);

// implemented in timesyshms.cpp
void timesyshms(Clock duration, struct tms *start, struct tms *end, Real idct, FILE *logFile);
