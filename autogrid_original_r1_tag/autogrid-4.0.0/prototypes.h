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

#ifndef _WIN32
#include <sys/times.h>
#endif

#include <sys/types.h>
#include <stdio.h>
#include "parameters.h"


void	banner( double version_num );
int	    setflags( int argc, char **argv );
ParameterEntry * apm_find( const char key[] );
void    apm_enter( const char key[], ParameterEntry value );
int	    check_size( int nelements, char axischar );
int	    gpfparser( char line[LINE_LEN] );
int	    main( int argc, char **argv );
int	    parsetypes(char * line, char *words[], int maxwords);
void	prHMSfixed( float t );
void	printdate( FILE *fp, int flag );
void	printhms( float t );
int	    strindex( char s[], char t[] );
void	timesys( Clock duration, struct tms *start, struct tms *end );
void	timesyshms( Clock duration, struct tms *start, struct tms *end );

/* EOF */
