/*

 $Id: printdate.cpp,v 1.7 2007/05/03 20:46:06 garrett Exp $

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

#include <stdio.h>
#include <sys/types.h>

#ifndef _WIN32
#   include <sys/time.h>
#else
#   include "times.h"
#endif

#ifdef HAVE_CONFIG_H
#   include <config.h>
#endif

#include "printdate.h"

void printdate( FILE *fp, int flag )
{
    time_t tn; /* tn = "time_now" */
    char *StringTimeDate;
    struct tm *ts;

    tn = time( &tn );

    ts = localtime( &tn );
    
    if (flag==1) {
        fprintf(fp, "%d:%02d %02d\" %s, %02d/%02d/%4d\n", 
        ( (ts->tm_hour >  12) ? (ts->tm_hour-12) : ts->tm_hour ), ts->tm_min, ts->tm_sec,
        ( (ts->tm_hour >= 12) ? "p.m." : "a.m." ),
        (ts->tm_mon + 1), ts->tm_mday, 1900+ts->tm_year );
    } else if (flag==2) {
          StringTimeDate = ctime( &tn );
          fprintf(fp, "%s", StringTimeDate);
    } else {
        fprintf(fp, "%d:%02d %02d\" %s\n", 
        ( (ts->tm_hour >  12) ? (ts->tm_hour-12) : ts->tm_hour ), ts->tm_min, ts->tm_sec,
        ( (ts->tm_hour >= 12) ? "pm" : "am" ) );
    }
}
/* EOF */
